# ============================================
# GUTS vs WORKER — within each GENUS, per SITE
# Macrotermes & Nitiditermes
# NMDS (facet Site x Genus) + PERMANOVA(PoolType) per panel
# ============================================

library(tidyverse)
library(vegan)
library(phyloseq)
library(ggrepel)

FEATURE_TSV <- "feature-table.tsv"
TAXON_TSV   <- "taxonomy.tsv"
TREE_NWK    <- "tree.nwk"
META_TSV    <- "sample-metadata.tsv"

# ============================================
# 1) Read feature table (BIOM->TSV)
# ============================================
otu_raw <- read.delim(
  FEATURE_TSV,
  skip = 1,
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

asv_col <- names(otu_raw)[1]     # usually "#OTU ID"
rownames(otu_raw) <- trimws(otu_raw[[asv_col]])
otu_raw[[asv_col]] <- NULL

otu_raw[] <- lapply(otu_raw, function(x) as.numeric(as.character(x)))
otu_raw[is.na(otu_raw)] <- 0

cat("OTU dims (taxa x samples): ", dim(otu_raw), "\n")
cat("OTU taxa head: ", paste(head(rownames(otu_raw), 5), collapse=" | "), "\n")
cat("OTU sample head: ", paste(head(colnames(otu_raw), 5), collapse=" | "), "\n")

# ============================================
# 2) Read taxonomy.tsv (robust column detection)
# ============================================
tax <- read.delim(TAXON_TSV, sep="\t", header=TRUE, check.names=FALSE)

feature_col <- if ("Feature.ID" %in% names(tax)) "Feature.ID" else
  if ("Feature ID" %in% names(tax)) "Feature ID" else
    names(tax)[1]

taxon_col <- if ("Taxon" %in% names(tax)) "Taxon" else names(tax)[2]

parse_q2_tax <- function(x){
  ranks <- c("Kingdom","Phylum","Class","Order","Family","Genus","Species")
  out <- rep(NA_character_, length(ranks))
  parts <- strsplit(x, ";\\s*")[[1]]
  parts <- gsub("^D_[0-6]__", "", parts)
  out[seq_len(min(length(parts), 7))] <- parts[seq_len(min(length(parts), 7))]
  out
}

tax_mat <- t(sapply(tax[[taxon_col]], parse_q2_tax))
rownames(tax_mat) <- trimws(tax[[feature_col]])
tax_mat <- as.matrix(tax_mat)

cat("TAX taxa head: ", paste(head(rownames(tax_mat), 5), collapse=" | "), "\n")

# ============================================
# 3) Align taxa between OTU and taxonomy
# ============================================
common_taxa <- intersect(rownames(otu_raw), rownames(tax_mat))
cat("Common taxa: ", length(common_taxa), "\n")

if(length(common_taxa) == 0){
  stop("Common taxa = 0. OTU and taxonomy do not match (or wrong file).")
}

otu2 <- otu_raw[common_taxa, , drop=FALSE]
tax2 <- tax_mat[common_taxa, , drop=FALSE]

OTU <- otu_table(as.matrix(otu2), taxa_are_rows = TRUE)
TAX <- tax_table(as.matrix(tax2))

# ============================================
# 4) Read metadata + create PoolType (guts/worker) + filter genera
# ============================================
meta <- read.delim(META_TSV, sep="\t", header=TRUE, comment.char="", check.names=FALSE)
meta <- meta[meta$`sample-id` != "#q2:types", ]

# PoolType : guts vs worker (case-insensitive)
meta$PoolType <- ifelse(grepl("guts", tolower(meta$`Nmb-indiv-Pools`)), "guts", "worker")
meta$PoolType <- factor(meta$PoolType, levels = c("guts","worker"))

genera_keep <- c("Macrotermes", "Nitiditermes")

meta2 <- meta %>%
  filter(genus %in% genera_keep) %>%
  mutate(site = factor(site),
         genus = factor(genus))

rownames(meta2) <- meta2$`sample-id`

# ============================================
# 5) Align samples between OTU and metadata
# ============================================
common_samples <- intersect(colnames(otu2), rownames(meta2))
cat("Common samples (all pools, selected genera): ", length(common_samples), "\n")

if(length(common_samples) < 6){
  stop("Too few samples after filtering. Check sample IDs and metadata.")
}

otu_s  <- otu2[, common_samples, drop=FALSE]
meta_s <- meta2[common_samples, , drop=FALSE]

OTU_s  <- otu_table(as.matrix(otu_s), taxa_are_rows = TRUE)
META_s <- sample_data(meta_s)

tree <- read_tree(TREE_NWK)

ps <- phyloseq(OTU_s, TAX, META_s, tree)
ps

# ============================================
# 6) Transform + Bray distance (Hellinger)
# ============================================
ps_rel <- transform_sample_counts(ps, function(x) x / sum(x))

otu_mat <- t(as(otu_table(ps_rel), "matrix"))  # samples x taxa
otu_hel <- decostand(otu_mat, method="hellinger")

df <- as(sample_data(ps_rel), "data.frame")

dist_bray <- vegdist(otu_hel, method="bray")

# ============================================
# 7) PERMANOVA within each panel (Site x Genus): PoolType effect
# ============================================
panels <- expand.grid(
  site  = levels(df$site),
  genus = levels(df$genus),
  stringsAsFactors = FALSE
)

panel_stats <- purrr::pmap_dfr(panels, function(site, genus){
  
  idx <- rownames(df)[df$site == site & df$genus == genus]
  df_p  <- df[idx, , drop=FALSE]
  otu_p <- otu_hel[idx, , drop=FALSE]
  
  # conditions minimales : 2 groupes (guts/worker) et >= 4 échantillons
  if(nrow(df_p) < 4 || length(unique(df_p$PoolType)) < 2){
    return(tibble(site=site, genus=genus, n=nrow(df_p), R2=NA_real_, p=NA_real_))
  }
  
  dist_p <- vegdist(otu_p, method="bray")
  a <- adonis2(dist_p ~ PoolType, data=df_p)
  
  tibble(
    site = site,
    genus = genus,
    n = nrow(df_p),
    R2 = a$R2[1],
    p  = a$`Pr(>F)`[1]
  )
})

print(panel_stats)

# Ajout d'étoiles de significativité
panel_stats <- panel_stats %>%
  mutate(sig = case_when(
    is.na(p) ~ "NA",
    p < 0.001 ~ "***",
    p < 0.01  ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ "ns"
  ))

# ============================================
# 8) NMDS plot (facet Site x Genus) + labels (p/R2)
# ============================================
set.seed(42)
ord <- metaMDS(dist_bray, k=2, trymax=200)

site_scores <- as.data.frame(vegan::scores(ord, display = "sites"))
site_scores$sample <- rownames(site_scores)

if (!all(c("NMDS1","NMDS2") %in% colnames(site_scores))) {
  colnames(site_scores)[1:2] <- c("NMDS1","NMDS2")
}

scores_df <- site_scores %>%
  left_join(df %>% rownames_to_column("sample"), by="sample")

# Positionnement du label dans chaque panneau (coin haut-gauche)
label_pos <- scores_df %>%
  group_by(site, genus) %>%
  summarise(
    x = min(NMDS1, na.rm=TRUE),
    y = max(NMDS2, na.rm=TRUE),
    .groups="drop"
  )

lab_df <- panel_stats %>%
  mutate(label = paste0("PERMANOVA(PoolType)\nR2=", round(R2, 3),
                        "\np=", signif(p, 3), "  ", sig)) %>%
  left_join(label_pos, by=c("site","genus"))

p <- ggplot(scores_df, aes(NMDS1, NMDS2, color = PoolType, shape = PoolType)) +
  geom_point(size = 4) +
  
  # Labels ID des échantillons
  geom_text_repel(
    aes(label = sample),
    size = 3,
    max.overlaps = 100
  ) +
  
  facet_grid(site ~ genus, scales="free") +
  theme_bw() +
  ggtitle("Guts vs Worker within genus — per site (Bantaco vs Saraya)") +
  
  # Ajout des p-values
  geom_text(
    data = lab_df,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0, vjust = 1,
    size = 3.3
  )

print(p)

# Export
ggsave("NMDS_guts_vs_worker_by_site_and_genus.pdf", p, width=11, height=6)
write.csv(panel_stats, "PERMANOVA_guts_vs_worker_by_site_and_genus.csv", row.names = FALSE)
