# 1. Set the working directory ----
setwd("~/Desktop/Elena/Trabajo/IDIBAPS/Roberta")

# 2. Upload packages of functions ----
library(dplyr)
library(tidyr)
library(readxl)

# 3.1 Upload raw data ----
library(readxl)
all_probes <- read_excel("Programming_FDR_Adapted_Elena.xlsx", skip = 1)

cd_new <- read_excel("Programming_FDR_Adapted_Elena.xlsx", skip = 1,
                                            sheet = "New cd0vcd21 fdr005 fc15 exp20")
sum(is.na(cd_new)) # 326665
cd_new <- cd_new[-which(is.na(cd_new$gene...1)==TRUE),]
sum(is.na(cd_new)) # 0
dim(cd_new) # 3342   58

hfd_new <- read_excel("Programming_FDR_Adapted_Elena.xlsx", skip = 1,
                      sheet = "New hfd0hfd21 fdr005 fc15 Exp20")
sum(is.na(hfd_new)) # 880507
hfd_new <- hfd_new[-which(is.na(hfd_new$gene...1)==TRUE),]
sum(is.na(hfd_new)) # 0
dim(hfd_new) # 4673   58

probes_l <- c('Tbx3', 'Nkx2-1', 'Ascl1', 'Isl1', 'Sox2', 'Shh', 'Rax', 'Nhlh2', #'Mash1',
             'Notch1', 'Notch2', 'Notch3', 'Notch4','Rbpj', 'Hes1', 'Neurog2')
probes_s <- c('Tbx3', 'Nkx2-1', 'Isl1', 'Ascl1')
#paper <- c('Pomc','Agrp','Aldh1l1','Nkx2-1','Isl1','Mash1','Tbx3', 'Th')

fc2_s <- all_probes[all_probes$gene...1 %in% probes_s,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
fc2_l <- all_probes[all_probes$gene...1 %in% probes_l,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
#cd_s <- cd_new[cd_new$gene...1 %in% probes_s,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
#hfd_s <- hfd_new[hfd_new$gene...1 %in% probes_s,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
#cd_l <- cd_new[cd_new$gene...1 %in% probes_l,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
#hfd_l <- hfd_new[hfd_new$gene...1 %in% probes_l,c(1,51:52)] ## SAME list than Short_list_of_genes_of_POMC_identity :)

names(fc2_s) <- c('Genes','FC_CD0vsCD21','FC_HFD0vsHFD21')
names(fc2_l) <- c('Genes','FC_CD0vsCD21','FC_HFD0vsHFD21')

# 4. Visualisations ----
library(gplots)
library(RColorBrewer)
library(pheatmap)
library(ggrepel)
library(ggpubr)

matriu <- fc2_l %>%
  group_by_all()%>%
  as.matrix()
dim(matriu) # 40    8
matriu <- as.matrix(apply(matriu[,-1],2,as.numeric))
#matriu <- t(matriu)
#colnames(matriu) <- fc2_l$Genes
rownames(matriu) <- fc2_l$Genes
#matriu <- matriu[c(2,7),]
#matriu <- matriu[c(2,3:5),]
dim(matriu) # 4 40

library(pheatmap)
pdf(file = "Robertas HFD0 vs HFD21 long list of Pomc identity_vertical.pdf", 
    #pdf(file = "Robertas HFD signif genes in common with Dropseq pomc population pools.pdf", 
    width = 10, height = 15)
my_palette <- colorRampPalette(c("dodgerblue4", "white", "firebrick"))(n = 1000)
paletteLength <- 1000
myBreaks <- c(seq(min(matriu), 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(max(matriu)/paletteLength, max(matriu), length.out=floor(paletteLength/2)))

pheatmap(matriu,  
         cellwidth = 14, cellheight = 22, width =10, height = 10, 
         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
         #cellwidth = 22, cellheight = 0.5, width =10, height = 10, 
         #fontsize = 10,fontsize_row = 0.05, fontsize_col = 10, 
         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
         legend = TRUE, border_color = 'gray', drop_levels = TRUE, breaks=myBreaks,
         main = "Overlapping with DropSeq_HFD. \n FCs. Not scaled.", 
         scale = "none",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
         cluster_rows= T, cluster_cols = F) 
#pheatmap(matriu,  
#         cellwidth = 14, cellheight = 22, width =10, height = 10, 
#         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
#         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
#         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
#         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
#         legend = TRUE, border_color = NA, drop_levels = TRUE, 
#         main = "Overlapping with DropSeq_HFD. \n Means of Pomc's FCs. Scaled by rows (differences between genes).", 
#         scale = "row",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
#         cluster_rows= F, cluster_cols = T) 
#pheatmap(matriu,  
#         cellwidth = 14, cellheight = 22, width =10, height = 10, 
#         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
#         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
#         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
#         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
#         legend = TRUE, border_color = NA, drop_levels = TRUE, 
#         main = "Overlapping with DropSeq_HFD. \n Means of Pomc's FCs. Scaled by columns (differences between treatments).", 
#         scale = "column",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
#         cluster_rows= F, cluster_cols = T) 
dev.off()

matriu <- fc2_s %>%
  group_by_all()%>%
  as.matrix()
dim(matriu) # 40    8
matriu <- as.matrix(apply(matriu[,-1],2,as.numeric))
matriu <- t(matriu)
colnames(matriu) <- fc2_s$Genes
#matriu <- matriu[c(2,7),]
#matriu <- matriu[c(2,3:5),]
dim(matriu) # 4 40

library(pheatmap)
pdf(file = "Robertas HFD0 vs HFD21 short list of Pomc identity.pdf", 
    #pdf(file = "Robertas HFD signif genes in common with Dropseq pomc population pools.pdf", 
    width = 15, height = 10)
my_palette <- colorRampPalette(c("dodgerblue4", "white", "firebrick"))(n = 1000)
#paletteLength <- 1000
#myBreaks <- c(seq(min(matriu), 0, length.out=ceiling(paletteLength/2) + 1), 
#              seq(max(matriu)/paletteLength, max(matriu), length.out=floor(paletteLength/2)))
pheatmap(matriu,  
         cellwidth = 14, cellheight = 22, width =10, height = 10, 
         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
         #cellwidth = 22, cellheight = 0.5, width =10, height = 10, 
         #fontsize = 10,fontsize_row = 0.05, fontsize_col = 10, 
         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
         legend = TRUE, border_color = 'gray', drop_levels = TRUE, #breaks=myBreaks,
         main = "Overlapping with DropSeq_HFD. \n FCs. Not scaled.", 
         scale = "none",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
         cluster_rows= F, cluster_cols = T) 
#pheatmap(matriu,  
#         cellwidth = 14, cellheight = 22, width =10, height = 10, 
#         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
#         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
#         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
#         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
#         legend = TRUE, border_color = NA, drop_levels = TRUE, 
#         main = "Overlapping with DropSeq_HFD. \n Means of Pomc's FCs. Scaled by rows (differences between genes).", 
#         scale = "row",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
#         cluster_rows= F, cluster_cols = T) 
#pheatmap(matriu,  
#         cellwidth = 14, cellheight = 22, width =10, height = 10, 
#         fontsize = 10,fontsize_row = 10, fontsize_col = 10, 
#         color = my_palette, #annotation_colors = ann_colors,#annotation = annotation_col,
#         show_colnames = TRUE,labels_col = colnames(matriu),#annotation_col = annotation_col[,-3], 
#         show_rownames = TRUE, labels_row = rownames(matriu), #annotation_row = gex2.means$Probes...1, 
#         legend = TRUE, border_color = NA, drop_levels = TRUE, 
#         main = "Overlapping with DropSeq_HFD. \n Means of Pomc's FCs. Scaled by columns (differences between treatments).", 
#         scale = "column",angle_col = "45",treeheight_col = 15,treeheight_row = 15, 
#         cluster_rows= F, cluster_cols = T) 
dev.off()

