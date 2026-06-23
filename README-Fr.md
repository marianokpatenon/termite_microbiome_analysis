# Analyse du microbiote intestinal des termites avec QIIME2 et R

## Présentation du projet

Ce dépôt contient l'ensemble des scripts bioinformatiques et statistiques utilisés pour analyser les communautés microbiennes associées aux termites. Les données de séquençage sont traitées à l'aide de **QIIME2**, puis analysées sous **R** afin d'étudier la diversité, la structure et la composition taxonomique du microbiote.

L'objectif principal est de caractériser les variations du microbiote en fonction :

* du genre de termite ;
* du site d'échantillonnage ;
* du type d'échantillon (tube digestif ou individu entier).

---

## Structure du dépôt

```text
.
├── qiime2_pipeline.sh
├── analysis_guts_worker.R
├── analysis_workers_site.R
├── README.md
└── results/
```

### Description des fichiers

| Fichier                   | Description                                                                                                                                     |
| ------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------- |
| `qiime2_pipeline.sh`      | Pipeline QIIME2 permettant de traiter les données brutes de séquençage jusqu'à la génération des tables d'ASV et des assignations taxonomiques. |
| `analysis_guts_worker.R`  | Analyses statistiques comparant le microbiote des tubes digestifs et celui des individus entiers.                                               |
| `analysis_workers_site.R` | Analyses statistiques comparant les microbiotes entre les différents sites d'échantillonnage.                                                   |
| `README.md`               | Documentation du projet.                                                                                                                        |

---

## Contexte scientifique

Les termites jouent un rôle majeur dans le fonctionnement des écosystèmes tropicaux grâce à leur capacité à dégrader la matière organique. Cette capacité repose en grande partie sur les communautés microbiennes présentes dans leur tube digestif.

L'étude du microbiote intestinal des termites permet de mieux comprendre :

* les interactions hôte-microorganismes ;
* les mécanismes de dégradation de la matière végétale ;
* l'écologie microbienne associée aux insectes sociaux ;
* les facteurs environnementaux influençant la structure des communautés microbiennes.

---

## Questions biologiques étudiées

Ce travail vise notamment à répondre aux questions suivantes :

1. Le microbiote diffère-t-il entre les genres de termites étudiés ?
2. Les communautés microbiennes varient-elles selon les sites d'échantillonnage ?
3. Existe-t-il des différences entre le microbiote des tubes digestifs et celui des individus entiers ?
4. Quels groupes microbiens sont responsables des différences observées ?

---

## Traitement bioinformatique des séquences

### Pipeline QIIME2

Les données de séquençage Illumina sont traitées selon les étapes suivantes :

1. Importation des données brutes.
2. Contrôle qualité des séquences.
3. Dénombrement des ASV (Amplicon Sequence Variants).
4. Élimination des chimères.
5. Construction de la table d'abondance.
6. Attribution taxonomique.
7. Calcul des indices de diversité.
8. Exportation des résultats vers R.

### Sorties principales

* Table d'ASV.
* Séquences représentatives.
* Assignations taxonomiques.
* Matrices de diversité alpha.
* Matrices de diversité bêta.

---

## Analyses de diversité alpha

Les indices suivants sont calculés :

* Richesse observée (Observed ASVs)
* Indice de Shannon
* Indice de Simpson
* Faith's Phylogenetic Diversity (Faith-PD)

Les comparaisons statistiques sont réalisées à l'aide :

* du test de Kruskal-Wallis ;
* de comparaisons post-hoc lorsque nécessaire.

---

## Analyses de diversité bêta

Les différences de composition microbienne sont étudiées à l'aide de :

### Distances écologiques

* Bray-Curtis
* Jaccard
* UniFrac pondérée
* UniFrac non pondérée

### Méthodes d'ordination

* PCoA (Principal Coordinates Analysis)
* NMDS (Non-metric Multidimensional Scaling)

### Tests statistiques

* PERMANOVA
* PERMDISP

---

## Plan d'échantillonnage

### Genres étudiés

* Macrotermes
* Nitiditermes

### Sites d'échantillonnage

* Bantaco
* Saraya

### Types d'échantillons

* Tubes digestifs (Gut)
* Individus ouvriers (Worker)

---

## Packages R utilisés

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

## Exécution des analyses

### Traitement des données de séquençage

```bash
bash qiime2_pipeline.sh
```

### Comparaison Gut vs Worker

```r
source("analysis_guts_worker.R")
```

### Comparaison des sites d'échantillonnage

```r
source("analysis_workers_site.R")
```

---

## Résultats produits

Les scripts génèrent notamment :

### Diversité alpha

* Boxplots
* Résultats statistiques
* Comparaisons entre groupes

### Diversité bêta

* Graphiques PCoA
* Graphiques NMDS
* Résultats PERMANOVA

### Composition taxonomique

* Barplots d'abondance relative
* Tableaux taxonomiques
* Visualisation des principaux taxons

---

## Reproductibilité

L'ensemble du pipeline est reproductible à partir :

* des données brutes de séquençage ;
* du fichier de métadonnées ;
* des scripts fournis dans ce dépôt.

Environnement recommandé :

* Linux Ubuntu
* QIIME2
* R ≥ 4.3

---

## Perspectives

Ce pipeline peut être facilement adapté à d'autres études de microbiote d'insectes ou d'organismes non modèles. Il constitue une base robuste pour l'exploration de la diversité microbienne, l'étude de la structuration des communautés et l'identification de taxons d'intérêt écologique.

---

## Auteur

**Mariano Joly KPATENON**

Docteur en Génétique et Génomique
Université de Montpellier

Spécialités :

* Génétique des populations
* Genome-Wide Association Studies
* Génomique évolutive
* Bioinformatique
* Analyses statistiques de données omiques

---

**Mots-clés :** Microbiote, Termites, QIIME2, Phyloseq, Diversité alpha, Diversité bêta, PERMANOVA, Écologie microbienne, Bioinformatique.
