# SIBO 16S rRNA Microbiome Analysis Pipeline

This repository contains the reproducible bioinformatics and statistical analysis pipeline developed for a Master's Thesis investigating differences in the duodenal microbiota between SIBO-positive and SIBO-negative samples using 16S rRNA gene sequencing data.

The analysis uses publicly available sequencing data from **BioProject PRJNA1273190** and implements an R-based workflow for read processing, amplicon sequence variant (ASV) inference, taxonomic assignment, microbial diversity analysis, and differential abundance testing.

## Study cohort

The final cohort comprises **123 samples**:

| SIBO status | Hypothyroidism negative | Hypothyroidism positive | Total |
|---|---:|---:|---:|
| SIBO negative | 40 | 29 | 69 |
| SIBO positive | 40 | 14 | 54 |
| **Total** | **80** | **43** | **123** |

Hypothyroidism was included as a covariate in the adjusted statistical analyses.

## Analysis workflow

The main R script performs the following steps:

1. Import and validation of paired-end FASTQ files.
2. Inspection of sequencing quality profiles.
3. Read filtering and truncation with DADA2.
4. Estimation of forward and reverse error models.
5. Dereplication and ASV inference.
6. Merging of forward and reverse reads.
7. Construction of the sequence table.
8. Chimera removal.
9. Taxonomic assignment using SILVA.
10. Integration of sample metadata and construction of a `phyloseq` object.
11. Alpha diversity analysis.
12. Beta diversity analysis using Bray-Curtis dissimilarities and principal coordinates analysis (PCoA).
13. Assessment of multivariate dispersion homogeneity using PERMDISP.
14. Unadjusted and hypothyroidism-adjusted PERMANOVA.
15. Descriptive genus-level relative abundance analysis.
16. Differential abundance analysis using ANCOM-BC2, including prevalence filtering and adjustment for hypothyroidism.
17. Automated export of tables, figures, and reproducibility information.

## Repository structure

The pipeline expects the following directory structure:

```text
TFM-SIBO-16S/
├── README.md
├── .gitignore
├── scripts/
│   └── 01_SIBO_pipeline.R
├── data/
│   ├── fastq/
│   │   ├── SIBO_POS/
│   │   ├── SIBO_NEG/
│   │   └── CONTROL/
│   ├── metadata/
│   │   └── SraRunTable.csv
│   └── silva/
│       ├── silva_nr99_v138.2_toGenus_trainset.fa.gz
│       └── silva_v138.2_assignSpecies.fa.gz
└── results/
    ├── figures/
    ├── tables/
    ├── processed/
    └── reproducibility/
```

Output directories are created automatically by the pipeline when they do not already exist.

## Input data

### Raw sequencing data

Raw FASTQ files are **not included in this repository** because of their size. The sequencing data originate from the publicly available **BioProject PRJNA1273190**.

The pipeline expects paired-end FASTQ files following this naming convention:

```text
<sample>_1.fastq
<sample>_2.fastq
```

Files must be placed in the corresponding `SIBO_POS`, `SIBO_NEG`, or `CONTROL` directory according to the study cohort used in the analysis.

### Metadata

Sample metadata are expected at:

```text
data/metadata/SraRunTable.csv
```

The metadata file must contain the sample identifiers and the clinical variables required to define SIBO status and hypothyroidism status.

### Taxonomic reference database

Taxonomic assignment was performed using **SILVA v138.2**. The pipeline expects the following reference files:

```text
data/silva/silva_nr99_v138.2_toGenus_trainset.fa.gz
data/silva/silva_v138.2_assignSpecies.fa.gz
```

These reference files are not distributed with this repository and must be obtained separately.

## Running the pipeline

The final analysis was performed using **R 4.6.0**.

Before running the pipeline, set the working directory to the root directory of the repository. The script defines the project directory using:

```r
project_dir <- normalizePath(getwd(), winslash = "/", mustWork = FALSE)
```

Once the required input files have been placed in the expected directories, the complete analysis can be run with:

```r
source("scripts/01_SIBO_pipeline.R")
```

Running the pipeline from a **clean R session** is recommended.

## Main software dependencies

The main packages used in the final analysis include:

- DADA2 1.40.0
- phyloseq 1.56.0
- vegan 2.7-5
- ANCOMBC 2.14.0
- microbiome 1.34.0
- ggplot2 4.0.3

The complete computational environment and installed package versions are automatically exported to:

```text
results/reproducibility/sessionInfo.txt
results/reproducibility/package_versions.csv
```

## Main analysis parameters

DADA2 filtering was performed using:

```r
truncLen = c(280, 220)
maxN = 0
maxEE = c(2, 2)
truncQ = 2
rm.phix = TRUE
```

Beta diversity analyses were based on **Bray-Curtis dissimilarities**. PERMDISP and PERMANOVA analyses used **999 permutations**.

For ANCOM-BC2, taxa were filtered using a minimum prevalence threshold of **10%** (`prv_cut = 0.10`). The effect of SIBO status was first evaluated in an unadjusted model and subsequently in a model including hypothyroidism as a covariate.

## Output

The pipeline automatically exports the main analysis outputs to the `results/` directory, including:

- DADA2 read-tracking tables;
- alpha diversity results;
- Bray-Curtis PCoA;
- PERMDISP and PERMANOVA results;
- genus-level relative abundance tables;
- unadjusted and adjusted ANCOM-BC2 results;
- figures used for interpretation and reporting;
- complete R session information and package versions.

## Reproducibility

The script provided in this repository corresponds to the pipeline used to generate the final results reported in the Master's Thesis.

Raw sequencing data are not redistributed through GitHub. Instead, the original public BioProject accession is provided together with the expected project structure and analysis code required to reproduce the workflow.

For exact software and package versions used in the final analysis, see `results/reproducibility/sessionInfo.txt`.

## Author

**Javier Sánchez Collado**

Master's Thesis, 2026
