# 1. Set the working directory ----
setwd("~/Desktop/Elena/Trabajo/IDIBAPS/Roberta")

# 2. Upload packages of functions ----
library(dplyr)
library(tidyr)
library(readxl)

# 3.1 Upload raw data ----

#Short_list_of_genes_of_POMC_identity <- read_excel("Short list of genes of POMC identity.xlsx")
#fc_s <- Short_list_of_genes_of_POMC_identity
#names(fc_s) <- c('Probes','FC_CD0vsCD21','FC_HFD0vsHFD21')
#rownames(fc_s) <-fc_s$Probes
#dim(fc_s) # 4  3

#all_probes <- read_excel("RNAseq_POMCRibotag_P0vsP21_HFD.xlsx", sheet = "All data", skip = 1)
all_probes <- read_excel("Programming_FDR_Adapted_Elena.xlsx", skip = 1)

#drop_seq_fc <- read_excel("Supl_Table_High_Fat_Diet_41593_2017_BFnn4495_MOESM74_ESM (1).xlsx")
drop_seq_fc <- read_excel("Supl_Table_High_Fat_Diet_41593_2017_BFnn4495_MOESM74_ESM (1).xlsx", 
                          sheet = "fold-change with FDR<25")

# 3.2 Pre-process data ----

drop_seq_fc <- drop_seq_fc[,c(1,15,16,21)]
names(drop_seq_fc) <- c('Gene','n14.Pomc/Ttr','n15.Pomc/Anxa2','n21.Pomc/Glipr1')
drop_seq_fc[,-1] <- apply(drop_seq_fc[,-1],1,as.numeric)
drop_seq_fc2 <- drop_seq_fc[-which(drop_seq_fc$`n14.Pomc/Ttr` == 0 & drop_seq_fc$`n15.Pomc/Anxa2` == 0 & drop_seq_fc$`n21.Pomc/Glipr1` == 0),]

#write.table(drop_seq_fc2, "Dropseq_significant_list.txt",sep="\t",dec = ".", row.names = F,col.names = T)

#probes_l <- c('Tbx3', 'Nkx2-1', 'Ascl1', 'Isl1', 'Sox2', 'Shh', 'Rax', 'Nhlh2', #'Mash1',
#              'Notch1', 'Notch2', 'Notch3', 'Notch4','Rbpj', 'Hes1', 'Neurog2')
#probes_s <- c('Tbx3', 'Nkx2-1', 'Isl1', 'Ascl1')
#fc2_s <- all_probes[all_probes$gene...1 %in% probes_s,c(1,53:56)] ## SAME list than Short_list_of_genes_of_POMC_identity :)
#names(fc2_s) <- c('Probes','FC_CD0vsCD21','FC_HFD0vsHFD21','FC_CD0vsHFD0','FC_CD21vsHFD21')

Venn_diagram_Dropseq_with_HFD <- read_excel("Venn diagram Dropseq with CD and HFD.xlsx", 
                                            sheet = "FCandFRD Dropseq vs P21 fdr025")
drop_seq <- Venn_diagram_Dropseq_with_HFD$'36 common elements in "HFD_fdr025" and "Drop-seq":'
sum(is.na(drop_seq)) # 1810
drop_seq <- drop_seq[-which(is.na(drop_seq)==TRUE)]
sum(is.na(drop_seq)) # 0
length(drop_seq) # 36

fc2_l <- all_probes[all_probes$gene...1 %in% drop_seq,c(1,51:52)]
names(fc2_l) <- c('Gene','FC_CD0vsCD21','FC_HFD0vsHFD21')
fc2_l <- left_join(fc2_l, drop_seq_fc, 'Gene')
sum(is.na(fc2_l)) # 0
fc2_l$pomc_means <- rowMeans(fc2_l[,4:6])
fc2_l$pomc_pool <- apply(fc2_l[,4:6],1,sum)
#fc2_l$pomc_median <- apply(fc2_l[,6:8],1,median)
#fc2_l[,-1]<- apply(fc2_l[,-1],1,as.numeric)
#write.table(fc2_l, "Robertas and dropseq_significant_list.txt",sep="\t",dec = ",", row.names = F,col.names = T)

#rownames(all_probes) <- all_probes$gene...1
#gex <- all_probes[,1:18]
#names(gex)<- c('Probes', 'CD0_C48','CD0_C51', 'CD0_C55', 'CD0_C63','CD0_C69','CD21_C42','CD21_C43','CD21_C79',
#               'HFD0_H10','HFD0_H14','HFD0_H18','HFD0_H19','HFD0_H57','HFD21_H26', 'HFD21_H35','HFD21_H39','HFD21_H49')
#write.table(fc, "FC all_data_CDandHFD.txt",sep="\t",dec = ".", row.names = F,col.names = T)
#paper <- c('Pomc','Agrp','Aldh1l1','Nkx2-1','Isl1','Mash1','Tbx3', 'Th')
#paper <- c('Mag','Gfap','Cdk5','Pmp22', 'Pomc','Agrp','Aldh1l1','Nkx2-1','Isl1','Mash1','Tbx3', 'Th')
#fc2_l <- fc[fc$Probes %in% probes_l,]

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
matriu <- t(matriu)
colnames(matriu) <- fc2_l$Gene
matriu <- matriu[c(2,7),]
matriu <- matriu[c(2,3:5),]
dim(matriu) # 4 40

library(pheatmap)
pdf(file = "Robertas HFD signif genes in common with Dropseq pomc data.pdf", 
#pdf(file = "Robertas HFD signif genes in common with Dropseq pomc population pools.pdf", 
        width = 15, height = 10)
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
        main = "Overlapping with DropSeq_HFD. \n Pomc's FCs. Not scaled.", 
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

