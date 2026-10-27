# Roberta Networks, DEPICT and GWAS

This repository contains R scripts and R Markdown analyses associated with the
[scientific paper](https://www.sciencedirect.com/science/article/pii/S2212877820300351).
The analyses include expression and trajectory visualizations of POMC neuron
RiboTag/RNA-seq data, examination of iron-metabolism genes, and human–mouse
homology and GWAS analyses.

## Repository contents

| File | Description |
|---|---|
| `Trajectories_P0_vs_P21.Rmd` | Pseudotime analysis comparing P0 and P21 samples. |
| `Trajectories_Ribotag_P0_P21_mice.Rmd` | Trajectory analysis of mouse RiboTag data. |
| `Plots_62_genes_Fe_metabolism_P0_P21_mice.Rmd` | Expression plots for 62 iron-metabolism-related genes. |
| `GWAS_homology_mmu_hsa.R` | Matches mouse and human homologs with obesity GWAS data and prepares gene/SNP lists for downstream analysis. |
| `Heatmap.R` and `Heatmap_new data .R` | Heatmap visualizations. |
| `Alluvial try.R` | Alluvial plot code (not published). |

## Data

Input datasets that are not included in this repository and paths need to be adapted for each user. 

In particular, the R Markdown files use an Excel workbook containing expression
data and metadata. The GWAS/homology analysis uses gene-homology and GWAS
association data. See the input-loading sections of each script for the
expected filenames, sheets, and formats.

## Requirements

R and the packages used by the individual scripts are required. Depending on
which analysis you run, these include packages such as `readxl`, `dplyr`,
`ggplot2`, `SingleCellExperiment`, `scater`, `scran`, `slingshot`, and
`tradeSeq`.

The repository does not currently provide a locked package environment, so
package versions may affect reproducibility.

## Publication

Please cite the associated paper when using these analyses:

Haddad-Tóvolli R, Altirriba J, Obri A, Sánchez EE, Chivite I, Milà-Guasch M, Ramírez S, Gómez-Valadés AG, Pozo M, Burguet J, Velloso LA, Claret M. Pro-opiomelanocortin (POMC) neuron translatome signatures underlying obesogenic gestational malprogramming in mice. Mol Metab. 2020 Jun;36:100963. doi: 10.1016/j.molmet.2020.02.006. Epub 2020 Feb 15. PMID: 32283518; PMCID: PMC7152705.

## Notes

These scripts document analyses associated with the publication; they are not
packaged as a standalone R pipeline. Results depend on the input data,
software versions, and any local path or parameter changes.
