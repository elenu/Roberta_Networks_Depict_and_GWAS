# 1. Set the working directory ----
setwd("~/IDIBAPS")

# 2. Upload packages of functions ----
library(dplyr)
library(tidyr)
library(readxl)

# 3.1 Upload raw data ----
library(readr)
mmu_hsa <- read_delim("C:/Users/eyre/Documents/IDIBAPS/Roberta/GWAS and Dropseq paper and analysis/depict_BI/depict_first_try/data/Complete List of Human and Mouse Homologs with phenotype annotations.csv", 
                      ";", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
library(readxl)
new_chow <- read_excel("C:/Users/eyre/Documents/IDIBAPS/Roberta/Programming_FDR_Adapted_Elena.xlsx", 
                       sheet = "New cd0vcd21 fdr005 fc15 exp20", 
                       skip = 1)
new_hfd <- read_excel("C:/Users/eyre/Documents/IDIBAPS/Roberta/Programming_FDR_Adapted_Elena.xlsx", 
                      sheet = "New hfd0hfd21 fdr005 fc15 Exp20", 
                      skip = 1)

library(readr)
All_ancestries_SNP_gwas_mc_merge_nogc_tbl_uniq <- read_delim("C:/Users/eyre/Documents/IDIBAPS/Roberta/GWAS and Dropseq paper and analysis/depict_BI/depict_first_try/data/All_ancestries_SNP_gwas_mc_merge_nogc.tbl.uniq.gz", 
                                                             "\t", escape_double = FALSE, trim_ws = TRUE)

gwas_catalog_v1_0_associations_e96_r2019_11_21_tsv <- read_delim("C:/Users/eyre/Documents/IDIBAPS/Roberta/GWAS and Dropseq paper and analysis/depict_BI/depict_first_try/data/gwas_catalog_v1.0.2-associations_e96_r2019-11-21.tsv", 
                                                                      "\t", escape_double = FALSE, trim_ws = TRUE)
# 3.2 Pre-process data ----
mmu_hsa_genes <- mmu_hsa[,c(1,5)]
names(mmu_hsa_genes) <- c('Genes', 'gene...1')
sum(is.na(mmu_hsa_genes)) # 0
dim(mmu_hsa_genes) # 18756     2

sum(is.na(new_chow)) # 326665
new_chow <- new_chow[-which(is.na(new_chow$gene...1)),]
sum(is.na(new_chow)) # 0
dim(new_chow) # 3342   58

sum(is.na(new_hfd)) # 880507
new_hfd <- new_hfd[-which(is.na(new_hfd$gene...1)),]
sum(is.na(new_hfd)) # 0
dim(new_hfd) # 4673    58

gwas <- gwas_catalog_v1_0_associations_e96_r2019_11_21_tsv
#m1 <- grep("diabetes", gwas$`DISEASE/TRAIT`,value = TRUE) # 2287
#m2 <- grep("Diabetes", gwas$`DISEASE/TRAIT`, value = TRUE) # 15
m1 <- grep("obesity", gwas$`DISEASE/TRAIT`,value = TRUE) # 48
m2 <- grep("Obesity", gwas$`DISEASE/TRAIT`, value = TRUE) # 1175
gwas_1 <- gwas[gwas$`DISEASE/TRAIT` %in% m1, ] 
gwas_2 <- gwas[gwas$`DISEASE/TRAIT` %in% m2, ] 
gwas <- rbind(gwas_1,gwas_2)
gwas <- gwas[,c(15,22,28)]
names(gwas) <- c('Genes', 'SNPs','P-value')
dim(gwas) # 1223      3

# 4. Intersect Roberta's list & list of mmu-hsa homologs & GWAS catalog ----
chow_mmu_hsa <- left_join(new_chow, mmu_hsa_genes, 'gene...1')
sum(is.na(chow_mmu_hsa)) # 103
chow_mmu_hsa <- chow_mmu_hsa[-which(is.na(chow_mmu_hsa$Genes)),]
sum(is.na(chow_mmu_hsa)) # 0
chow_mmu_hsa <- left_join(chow_mmu_hsa, gwas, 'Genes')
#chow_mmu_hsa <- left_join(gwas, chow_mmu_hsa, 'Genes')
sum(is.na(chow_mmu_hsa)) # 6318
sum(is.na(chow_mmu_hsa$Genes)) # 0
sum(is.na(chow_mmu_hsa$SNPs)) # 3159
chow_mmu_hsa <- chow_mmu_hsa[-which(is.na(chow_mmu_hsa$SNPs)),]
sum(is.na(chow_mmu_hsa)) # 0
#chow_mmu_hsa <- chow_mmu_hsa[new_chow$gene...1 %in% chow_mmu_hsa$gene...1,]
#sum(is.na(chow_mmu_hsa)) # 0
dim(chow_mmu_hsa) # 127   61
snps_chow <- unique(chow_mmu_hsa$SNPs)
length(snps_chow) # 5668
write.table(snps_chow, "chow_in_hsa_snps_list_obesity.txt",sep="\t",dec = ",", row.names = F,col.names = F)
# Remember to remove qotes from SNPs before submitting them to Depict!!!
write.table(chow_mmu_hsa$Genes, "chow_in_hsa_homologs_list_obesity.txt",sep="\t",dec = ",", row.names = F,col.names = T)

hfd_mmu_hsa <- left_join(new_hfd, mmu_hsa_genes, 'gene...1')
sum(is.na(hfd_mmu_hsa)) # 153
hfd_mmu_hsa <- hfd_mmu_hsa[-which(is.na(hfd_mmu_hsa$Genes)),]
sum(is.na(hfd_mmu_hsa)) # 0
hfd_mmu_hsa <- left_join(hfd_mmu_hsa, gwas, 'Genes')
sum(is.na(hfd_mmu_hsa)) # 8860
sum(is.na(hfd_mmu_hsa$Genes)) # 0
sum(is.na(hfd_mmu_hsa$SNPs)) # 4430
hfd_mmu_hsa <- hfd_mmu_hsa[-which(is.na(hfd_mmu_hsa$SNPs)),]
sum(is.na(hfd_mmu_hsa)) # 0
dim(hfd_mmu_hsa) # 144   61
snps_hfd <- unique(hfd_mmu_hsa$SNPs)
length(snps_hfd) # 115
write.table(snps_hfd, "hfd_in_hsa_snps_list_obesity.txt",sep="\t",dec = ",", row.names = F,col.names = F)
# Remember to remove quotes from SNPs before submitting them to Depict!!!
write.table(hfd_mmu_hsa$Genes, "hfd_in_hsa_homologs_list_obesity.txt",sep="\t",dec = ",", row.names = F,col.names = T)

# 5. Check results and save them as csv ----
# Chow-Obesity

