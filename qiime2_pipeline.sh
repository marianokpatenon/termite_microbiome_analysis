#!/usr/bin/env bash
# Utilise bash depuis l'environnement courant

set -euo pipefail
# -e : stoppe le script dès qu'une commande échoue
# -u : erreur si une variable non définie est utilisée
# -o pipefail : si une commande d'un pipe échoue, le pipe entier échoue

PROJECT_DIR="/home/pkpat/qiime2_project"
# Dossier projet = 1er argument du script, sinon le dossier courant

MANIFEST="${PROJECT_DIR}/Cleandata/manifest.tsv"
# Chemin vers le fichier manifest qui liste les FASTQ (nom échantillon + chemins fichiers)

METADATA="${PROJECT_DIR}/metadata/sample-metadata.tsv"
# Chemin vers le fichier de métadonnées (variables expérimentales par échantillon)

CLASSIFIER="${PROJECT_DIR}/classifiers/silva-16s-classifier.qza"
# Classifieur taxonomique QIIME2 (SILVA/GG2) au format .qza

TRIM_LEFT_F=0
TRIM_LEFT_R=0
# Nombre de bases à couper au début (5' primer/adapters) pour forward et reverse
# 0 signifie "ne rien couper au début"

TRUNC_LEN_F=290
TRUNC_LEN_R=290
# Longueur de troncature (coupure en fin de reads) pour forward et reverse
# Doit être choisie en regardant demux summarize (qualité)

SAMPLING_DEPTH=10000
# Profondeur de rarefaction utilisée pour les métriques de diversité
# À choisir selon la distribution des reads par échantillon (table-asv.qzv)

THREADS=10
# Threads = 0 => QIIME2 utilise tous les cœurs disponibles (souvent)

OUT="${PROJECT_DIR}/results"
# Dossier de sortie

mkdir -p "${OUT}"
# Crée le dossier results/ s'il n'existe pas

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path "${MANIFEST}" \
  --output-path "${OUT}/demux-paired.qza" \
  --input-format PairedEndFastqManifestPhred33V2
# Importe les FASTQ (paired-end) décrits par le manifest
# Produit un artefact QIIME2 .qza contenant les séquences + qualités

qiime demux summarize \
  --i-data "${OUT}/demux-paired.qza" \
  --o-visualization "${OUT}/demux-paired.qzv"
# Génère un rapport .qzv (qualité par position, nb reads, etc.)
# Sert à choisir TRUNC_LEN_* et éventuellement TRIM_LEFT_*

qiime dada2 denoise-paired \
  --i-demultiplexed-seqs "${OUT}/demux-paired.qza" \
  --p-trim-left-f "${TRIM_LEFT_F}" --p-trim-left-r "${TRIM_LEFT_R}" \
  --p-trunc-len-f "${TRUNC_LEN_F}" --p-trunc-len-r "${TRUNC_LEN_R}" \
  --p-n-threads "${THREADS}" \
  --o-table "${OUT}/table-asv.qza" \
  --o-representative-sequences "${OUT}/rep-seqs.qza" \
  --o-denoising-stats "${OUT}/denoising-stats.qza"
# Étape clé : DADA2
# - filtration/denoising/correction d'erreurs + merge paired-end + suppression chimères
# Sorties :
# - table-asv.qza : table d’abondance des ASV par échantillon
# - rep-seqs.qza : séquences des ASV
# - denoising-stats.qza : stats (reads filtrés, merged, chimères, etc.)

qiime feature-table summarize \
  --i-table "${OUT}/table-asv.qza" \
  --m-sample-metadata-file "${METADATA}" \
  --o-visualization "${OUT}/table-asv.qzv"
# Résume la table ASV : nb reads/échantillon, nb features, etc.
# Très utile pour choisir SAMPLING_DEPTH

qiime feature-classifier classify-sklearn \
  --i-classifier "${CLASSIFIER}" \
  --i-reads "${OUT}/rep-seqs.qza" \
  --p-n-jobs "${THREADS}" \
  --o-classification "${OUT}/taxonomy.qza"
# Assigne une taxonomie à chaque ASV (rep-seqs) via le classifieur (SILVA)
# Sortie : taxonomy.qza (table ASV -> taxon)

qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "${OUT}/rep-seqs.qza" \
  --o-rooted-tree "${OUT}/rooted-tree.qza"
# Aligne les séquences ASV, construit un arbre (FastTree) et le root
# Sortie : rooted-tree.qza (nécessaire pour UniFrac, Faith PD, core-metrics-phylogenetic)

qiime diversity core-metrics-phylogenetic \
  --i-phylogeny "${OUT}/rooted-tree.qza" \
  --i-table "${OUT}/table-asv.qza" \
  --p-sampling-depth "${SAMPLING_DEPTH}" \
  --m-metadata-file "${METADATA}" \
  --output-dir "${OUT}/core-metrics-results"
# Calcule les métriques de diversité alpha/beta + PCoA + Emperor
# IMPORTANT : rarefaction à SAMPLING_DEPTH

echo "[DONE] ${OUT}"
# Message de fin

