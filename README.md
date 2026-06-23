# termite_microbiome_analysis
# Termite Microbiome Analysis using QIIME2 and R

## Overview

This repository contains the bioinformatics and statistical workflows used to analyze the gut microbiota of termites. The analyses were performed using QIIME2 for sequence processing and taxonomic assignment, followed by R for diversity analyses, visualization, and statistical testing.

The project focuses on understanding the structure and diversity of termite-associated microbial communities across termite genera, sampling sites, and sample types.

---

## Repository Structure

```text
.
├── qiime2_pipeline.sh
├── analysis_guts_worker.R
├── analysis_workers_site.R
├── README.md
└── results/
```

### Files Description

| File                      | Description                                                             |
| ------------------------- | ----------------------------------------------------------------------- |
| `qiime2_pipeline.sh`      | QIIME2 workflow from raw sequences to ASV table and taxonomy assignment |
| `analysis_guts_worker.R`  | Statistical analyses comparing gut and worker microbiomes               |
| `analysis_workers_site.R` | Statistical analyses comparing worker microbiomes across sampling sites |
| `README.md`               | Project documentation                                                   |

---

## Biological Questions

This study aims to answer the following questions:

1. Does microbiome composition differ between termite genera?
2. Does microbiome diversity vary among sampling sites?
3. Are gut microbiomes different from whole-worker microbiomes?
4. Which microbial taxa contribute most to observed differences?

---

## Sequencing Data Processing

### QIIME2 Workflow

Raw paired-end Illumina reads were processed using QIIME2.

Main steps:

1. Import sequencing data
2. Quality assessment
3. Denoising and ASV inference
4. Chimera removal
5. Taxonomic assignment
6. Diversity analyses
7. Export of feature tables and taxonomic profiles

### Main Outputs

* Feature table (ASV abundance matrix)
* Representative sequences
* Taxonomic assignments
* Alpha diversity metrics
* Beta diversity distance matrices

---

## Alpha Diversity Analyses

The following metrics were evaluated:

* Shannon diversity index
* Simpson diversity index
* Observed ASVs
* Faith's Phylogenetic Diversity (Faith-PD)

Statistical comparisons were performed using:

* Kruskal-Wallis tests
* Pairwise Wilcoxon tests when appropriate

---

## Beta Diversity Analyses

Community composition differences were investigated using:

* Bray-Curtis distance
* Jaccard distance
* Weighted UniFrac
* Unweighted UniFrac

Ordination methods:

* Principal Coordinates Analysis (PCoA)
* Non-metric Multidimensional Scaling (NMDS)

Statistical tests:

* PERMANOVA
* PERMDISP

---

## Study Design

### Termite Genera

* Macrotermes
* Nitiditermes

### Sampling Sites

* Bantaco
* Saraya

### Sample Types

* Gut samples
* Worker samples

---

## R Packages

```r
library(phyloseq)
library(vegan)
library(ggplot2)
library(dplyr)
library(tidyr)
library(picante)
library(ape)
library(microbiome)
```

---

## Example Workflow

### QIIME2

```bash
bash qiime2_pipeline.sh
```

### Gut vs Worker Analysis

```r
source("analysis_guts_worker.R")
```

### Site Comparison Analysis

```r
source("analysis_workers_site.R")
```

---

## Main Outputs

### Diversity Results

* Alpha diversity boxplots
* Beta diversity ordinations
* Statistical comparison tables

### Community Structure

* Relative abundance barplots
* Heatmaps
* Taxonomic summaries

### Statistical Outputs

* PERMANOVA results
* Pairwise comparisons
* Diversity significance tests

---

## Reproducibility

All analyses were performed using:

* QIIME2
* R (≥ 4.3)
* Linux environment

The workflow is fully reproducible provided that the raw sequencing data and metadata files are available.

---

## Citation

If you use this workflow, please cite:

KPATENON MJ et al. Termite gut microbiome analyses using QIIME2 and R.

---

## Author

Mariano Joly KPATENON

PhD in Genetics and Genomics

University of Montpellier

Population Genetics • Evolutionary Biology • Microbial Ecology • Bioinformatics
