setwd("~/Desktop/Elena/Trabajo/IDIBAPS/Roberta")

#install.packages('ggfittext')
library(ggalluvial) # ggplot2 loaded too
library(ggfittext)
library(dplyr)
library(tidyr)
library(readxl)
library(gplots)
library(RColorBrewer)
library(ggrepel)

library(readxl)
RNAseq_POMCRibotag_P0vsP21_HFD <- read_excel("RNAseq_POMCRibotag_P0vsP21_HFD.xlsx", 
                                               sheet = "CD0vsCD21 FDR>0,05 FC>1,5Exp>20", 
                                               skip = 1)
rownames(RNAseq_POMCRibotag_P0vsP21_HFD) <- RNAseq_POMCRibotag_P0vsP21_HFD$gene...1
gex <- RNAseq_POMCRibotag_P0vsP21_HFD[,1:18]
names(gex)<- c('Probes', 
               'CD0_C48','CD0_C51', 'CD0_C55', 'CD0_C63','CD0_C69',
               'CD21_42','CD21_C43','CD21_79',
               'HFD0_H10','HFD0_H14','HFD0_H18','HFD0_H19','HFD0_H57',
               'HFD21_H26', 'HFD21_H35','HFD21_H39','HFD21_H49')
write.table(gex, "gex data_CD.txt",sep="\t",dec = ".", row.names = F,col.names = T)
rownames(gex) <- gex$Probes
gex_stats <- RNAseq_POMCRibotag_P0vsP21_HFD[,36:44]
names(gex_stats) <- c('Probes',
                      'CD0_mean','CD21_mean','HFD0_mean','HFD21_mean',
                      'CD0_sem','CD21_sem','HFD0_sem','HFD21_sem')
rownames(gex_stats) <- gex_stats$Probes
abs_pres <- RNAseq_POMCRibotag_P0vsP21_HFD[,c(1,45:48)]
names(abs_pres) <- c('Probes','CD0','CD21','HFD0','HFD21')
rownames(abs_pres) <- abs_pres$Probes
abs_pres_comp <- RNAseq_POMCRibotag_P0vsP21_HFD[,c(1,49:52)]
names(abs_pres_comp) <- c('Probes','CD0vsCD21',	'HFD0vsHFD21',	'CD0vsHFD0',	'CD21vsHFD21CD0')
rownames(abs_pres_comp) <- abs_pres_comp$Probes
foldchange <- RNAseq_POMCRibotag_P0vsP21_HFD[,c(1,53:56)]
names(foldchange) <- c('Probes','CD0vsCD21',	'HFD0vsHFD21',	'CD0vsHFD0',	'CD21vsHFD21CD0')
rownames(foldchange) <- foldchange$Probes
diet_t_test <- RNAseq_POMCRibotag_P0vsP21_HFD[,c(1,57:60)]
names(diet_t_test) <- c('Probes','CD0vsCD21',	'HFD0vsHFD21',	'CD0vsHFD0',	'CD21vsHFD21CD0')
rownames(diet_t_test) <- diet_t_test$Probes

diet_groups <- cbind(abs_pres, abs_pres_comp)
                      
is_alluvia_form(as.data.frame(diet_groups), axes = 1:3, silent = TRUE)

ggplot(as.data.frame(diet_groups),
       aes(y = abs_pres,
           axis1 = abs_pres_comp)) +
  geom_alluvium(aes(fill = abs_pres),
                width = 0, knot.pos = 0, reverse = FALSE) +
  guides(fill = FALSE) +
  geom_stratum(width = 1/8, reverse = FALSE) +
  geom_text(stat = "stratum", label.strata = TRUE, reverse = FALSE) +
#  scale_x_continuous(breaks = 1:3, labels = c("OS_bin", "MYCN", "New_clusters")) +
  coord_flip() +
  ggtitle("Training with alluvial plots")



ggplot(clinical,
       aes(x = New_clusters, stratum = MYCN, alluvium = AGE,
           fill = MYCN, label = MYCN)) +
  scale_fill_brewer(type = "qual", palette = "Set2") +
  geom_flow(stat = "alluvium", lode.guidance = "rightleft",
            color = "darkgray") +
  geom_stratum() +
  theme(legend.position = "bottom") +
  ggtitle("Training with alluvial plots")


ggplot(clinical,
       aes(x = STAGE, stratum = MYCN, alluvium = OS_bin,
           y = OS_d,
           fill = MYCN, label = MYCN)) +
  scale_x_discrete(expand = c(.1, .1)) +
  geom_flow() +
  geom_stratum(alpha = .5) +
  geom_text(stat = "stratum", size = 3) +
  theme(legend.position = "none") +
  ggtitle("Training with alluvial plots")

ggplot(data = clinical,
       aes(axis1 = STAGE, axis2 = MYCN, axis3 =AGE,axis4= New_clusters,
           y = sort)) +
  scale_x_discrete(limits = c("STAGE", "MYCN", "AGE","New Clusters"), expand = c(.1, .05)) +
  xlab("Clinical Parameters") +
  #geom_alluvium(aes(fill = as.factor(OS_bin),colour=as.factor(OS_bin))) +
  geom_alluvium(aes(fill = New_clusters,colour=New_clusters)) +
  #geom_alluvium(aes(fill = New_clusters,colour=c("Alpha"= "#00AFBB","Beta"= "green","Delta"= "orange", "Omega"="red"))) + ## ERROR
  geom_stratum() + 
  geom_text(stat = "stratum", label.strata = TRUE, size=3) +
  theme_minimal() +
  ggtitle("Kocak, GSE45547",
          "Distribution of the patients")


 
