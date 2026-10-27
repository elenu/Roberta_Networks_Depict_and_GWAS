# Roberta Networks, DEPICT and GWAS

This repository contains R scripts and R Markdown analyses associated with the
[scientific paper](https://www.sciencedirect.com/science/article/pii/S2212877820300351).
The analyses include expression and trajectory visualizations of POMC neuron
RiboTag/RNA-seq data, examination of iron-metabolism genes, and human–mouse
homology and GWAS analyses.

## Repository contents

| File | Description |
|---|---|
| `GWAS_homology_mmu_hsa.R` | Matches mouse and human homologs with obesity GWAS data and prepares gene/SNP lists for downstream analysis. |
| `Heatmap.R` and `Heatmap_new data .R` | Heatmap visualizations. |

## Data

Input datasets that are not included in this repository and paths need to be adapted for each user. 

In particular, the R Markdown files use an Excel workbook containing expression
data and metadata. The GWAS/homology analysis uses gene-homology and GWAS
association data. See the input-loading sections of each script for the
expected filenames, sheets, and formats.

## Requirements

These analyses are implemented in **R** and rely on CRAN and Bioconductor packages.
Because scripts are run independently, package needs vary by file.

### R version

- Recommended: **R >= 4.1** (R 4.2+ preferred)

### Core packages used across scripts

- `readxl` (Excel input)
- `dplyr`, `tidyr`, `tibble` (data wrangling)
- `stringr` (string handling)
- `ggplot2` (visualization)
- `pheatmap` and/or `ComplexHeatmap` (heatmaps; depends on script)
- `openxlsx` or `writexl` (table export; if enabled in script)

### Single-cell / trajectory analysis packages (used in trajectory-focused scripts)

- `SingleCellExperiment`
- `scater`
- `scran`
- `slingshot`
- `tradeSeq`

### Optional / script-specific packages

Some scripts may additionally require packages for:
- annotation and gene ID conversion (e.g., `biomaRt`, `AnnotationDbi`, organism-specific `org.*.eg.db`)
- GWAS summary table handling and joins
- plotting themes / label management (e.g., `cowplot`, `patchwork`, `ggrepel`)

Please check the `library(...)` calls at the top of each script for the exact set needed for that analysis.

## Requirements

### Install packages

CRAN packages:

```r
install.packages(c(
  "readxl", "dplyr", "tidyr", "tibble", "stringr",
  "ggplot2", "pheatmap", "openxlsx", "writexl"
))
```

Bioconductor packages:

```r
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c(
  "SingleCellExperiment", "scater", "scran", "slingshot", "tradeSeq"
))
```

### Reproducibility note

This repository currently does **not** include a pinned environment (`renv.lock` or container).
For reproducible results, we recommend recording `sessionInfo()` output for each run and, if possible, creating a local `renv` environment.

## Publication

Please cite the associated paper when using these analyses:

Haddad-Tóvolli R, Altirriba J, Obri A, Sánchez EE, Chivite I, Milà-Guasch M, Ramírez S, Gómez-Valadés AG, Pozo M, Burguet J, Velloso LA, Claret M. Pro-opiomelanocortin (POMC) neuron translatome signatures underlying obesogenic gestational malprogramming in mice. Mol Metab. 2020 Jun;36:100963. doi: 10.1016/j.molmet.2020.02.006. Epub 2020 Feb 15. PMID: 32283518; PMCID: PMC7152705.

## Notes

These scripts document analyses associated with the publication; they are not
packaged as a standalone R pipeline. Results depend on the input data,
software versions, and any local path or parameter changes.
