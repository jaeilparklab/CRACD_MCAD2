library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(sf)
library(dplyr)

options(future.globals.maxSize = 8000 * 1024^2)

## For the B6-transplanted PC tumor ##

setwd("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067741__Region_1__20250430__175356/analysis/clustering/gene_expression_graphclust_annotated")
cluster = read.csv("clusters.csv")
head(cluster)
cluster$group <- as.character(cluster$group)
cluster$group <- recode(cluster$group, `1` = "Basal/Nfib+_1", `2` = "Fibroblast_2", `3` = "Stratified/Nfib+", 
                               `4` = "Keratin_layer_1", `5` = "CD8 T_1/exhausted", `6` = "Basal_1", 
                               `7` = "CD4 T_1", `8` = "CD8 T_5",`9` = "Basal_2", `10` = "Fibroblast_3",
                            `11` = "Dendritic_1", `12` = "Fibroblast_1", `13` = "CD4 T_3", 
                            `14` = "CD4 T_6", `15` = "Stratified_3", `16` = "Keratinocyte_2", 
                            `17` = "B cell", `18` = "Fibroblast_4",`19` = "Stratified_1", `20` = "Proliferating_1",
                            `21` = "Basal/Nfib+_3", `22` = "Macrophage", `23` = "Stratified_2", 
                            `24` = "Keratinocyte_1", `25` = "Keratinocyte_4", `26` = "CD4 T_4", 
                            `27` = "Basal/Nfib+_2", `28` = "CD4 T_5",`29` = "Keratinocyte_6", `30` = "Keratinocyte_3",
                            `31` = "Keratinocyte_2", `32` = "CD8 T_3", `33` = "CD8 T_6",`34` = "CD8 T_7", `35` = "Keratinocyte_5", `36` = "Stratified_4",
                            `37` = "Stratified_5", `Unassigned` = "Unassigned"
                            )
write.csv(cluster, "cluster_annotated.csv")


## For the normal esophagus ##

setwd("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067741__Region5__20250430__175356/analysis/clustering/gene_expression_graphclust_annotated")
cluster = read.csv("clusters.csv")
head(cluster)
cluster <- rename(cluster, cell_id = Barcode, group = Cluster)
head(cluster)

cluster$group <- as.character(cluster$group)
cluster$group <- recode(cluster$group, `1` = "Basal_1", `2` = "Intermediate_supra", `3` = "Early_supra", 
                        `4` = "Acinar-like/Amy2a5+", `5` = "Fibroblast_1", `6` = "Stratified_1", 
                        `7` = "Stratified_2", `8` = "Endothelial_cell",`9` = "Fibroblast_2", `10` = "Smooth_muscle",
                        `11` = "Stratified_3", `12` = "Macrophage", `Unassigned` = "Unassigned"
)
write.csv(cluster, "cluster_annotated.csv")
