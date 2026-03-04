# ============================================
# WORKERS ONLY — Bantaco vs Saraya by Genus
# Macrotermes & Nitiditermes
# NMDS (facets) + PERMANOVA(site) per genus
# ============================================

library(tidyverse)
library(vegan)
library(phyloseq)

FEATURE_TSV <- "feature-table.tsv"
TAXON_TSV   <- "taxonomy.tsv"
TREE_NWK    <- "tree.nwk"
META_TSV    <- "sample-metadata.tsv"

# ============================================
# 1) Read feature table (BIOM->TSV)
# ============================================
otu_raw <- read.delim(
  FEATURE_TSV,
  skip = 1,                 # remove "# Constructed from biom file"
  sep = "\t",
  header = TRUE,
  check.names = FALSE
)

asv_col <- names(otu_raw)[1]      # usually "#OTU ID"
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
# 4) Read metadata + keep WORKERS only + keep 2 genera
# ============================================
meta <- read.delim(META_TSV, sep="\t", header=TRUE, comment.char="", check.names=FALSE)
meta <- meta[meta$`sample-id` != "#q2:types", ]

# Worker/Guts from "Nmb-indiv-Pools" values like "30guts" / "40worker" / "30Worker"
meta$Type <- ifelse(grepl("worker", tolower(meta$`Nmb-indiv-Pools`)), "Worker", "Guts")

genera_keep <- c("Macrotermes", "Nitiditermes")

meta_w <- meta %>%
  filter(Type == "Worker", genus %in% genera_keep) %>%
  mutate(site = factor(site),
         genus = factor(genus))

rownames(meta_w) <- meta_w$`sample-id`

# ============================================
# 5) Align samples between OTU and metadata
# ============================================
common_samples <- intersect(colnames(otu2), rownames(meta_w))
cat("Common samples (workers, selected genera): ", length(common_samples), "\n")

if(length(common_samples) < 4){
  stop("Too few samples after filtering. Check sample IDs and metadata.")
}

otu_w  <- otu2[, common_samples, drop=FALSE]
meta_w <- meta_w[common_samples, , drop=FALSE]

OTU_w  <- otu_table(as.matrix(otu_w), taxa_are_rows = TRUE)
META_w <- sample_data(meta_w)

# Tree (optional for this analysis, but we include it)
tree <- read_tree(TREE_NWK)

ps_w <- phyloseq(OTU_w, TAX, META_w, tree)
ps_w

# ============================================
# 6) Transform + Bray distance (Hellinger recommended)
# ============================================
ps_rel <- transform_sample_counts(ps_w, function(x) x / sum(x))

otu_mat <- t(as(otu_table(ps_rel), "matrix"))           # samples x taxa
otu_hel <- decostand(otu_mat, method="hellinger")

df <- as(sample_data(ps_rel), "data.frame")

dist_bray <- vegdist(otu_hel, method="bray")

# ============================================
# 7) PERMANOVA(site) per genus
# ============================================
perm_by_genus <- lapply(genera_keep, function(g){
  
  idx <- rownames(df)[df$genus == g]
  df_g  <- df[idx, , drop=FALSE]
  otu_g <- otu_hel[idx, , drop=FALSE]
  
  # Must have both sites and at least 3 samples
  if(length(unique(df_g$site)) < 2 || nrow(df_g) < 3){
    return(tibble(genus=g, n=nrow(df_g), R2=NA_real_, p=NA_real_))
  }
  
  dist_g <- vegdist(otu_g, method="bray")
  a <- adonis2(dist_g ~ site, data=df_g)
  
  tibble(
    genus = g,
    n     = nrow(df_g),
    R2    = a$R2[1],
    p     = a$`Pr(>F)`[1]
  )
}) %>% bind_rows()

print(perm_by_genus)

# ============================================
# 8) One NMDS plot (facetted by genus) + p-values
# ============================================

set.seed(42)
ord <- metaMDS(dist_bray, k=2, trymax=200)

# Extraire les coordonnées des sites (robuste)
site_scores <- as.data.frame(vegan::scores(ord, display = "sites"))
site_scores$sample <- rownames(site_scores)

# Renommer en NMDS1/NMDS2 si besoin
if (!all(c("NMDS1","NMDS2") %in% colnames(site_scores))) {
  colnames(site_scores)[1:2] <- c("NMDS1","NMDS2")
}

# Joindre avec le metadata df (sites/genus)
scores_df <- site_scores %>%
  left_join(df %>% rownames_to_column("sample"), by = "sample")

# Labels p-values/R2 par genus
lab_df <- perm_by_genus %>%
  mutate(label = paste0("PERMANOVA(site)\nR2=", round(R2, 3), "\np=", signif(p, 3))) %>%
  left_join(scores_df %>%
              group_by(genus) %>%
              summarise(x = min(NMDS1, na.rm=TRUE),
                        y = max(NMDS2, na.rm=TRUE),
                        .groups="drop"),
            by = "genus")

p <- ggplot(scores_df, aes(NMDS1, NMDS2, color = site)) +
  geom_point(size = 4) +
  facet_wrap(~ genus, scales = "free") +
  theme_bw() +
  ggtitle("Workers only — Bantaco vs Saraya (Macrotermes & Nitiditermes)") +
  geom_text(data = lab_df,
            aes(x = x, y = y, label = label),
            inherit.aes = FALSE, hjust = 0, vjust = 1, size = 3.5)

print(p)

# Optionnel : export
ggsave("NMDS_workers_site_by_genus.pdf", p, width=10, height=5)
