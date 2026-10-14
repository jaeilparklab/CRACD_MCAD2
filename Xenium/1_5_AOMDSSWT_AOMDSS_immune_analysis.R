rm(list=ls())

install.packages("arrow")
install.packages('Seurat')


library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(sf)
library(reshape2)


options(future.globals.maxSize = 8000 * 1024^2)


setwd("D:/KP/Colon_Xenium/")

Xenium_WT.obj = readRDS("WT_clustered.rds")
Xenium_KO.obj = readRDS("KO_clustered.rds")




##########################################################################################################################
############################      Tumor cells of A and CA integration ####################################################            
############################  Merge can be done in other workstation by Jinho ########################################### 


options(future.globals.maxSize = 2000000000000)

Tumor <- merge(Xenium_WT.obj, y=Xenium_KO.obj, add.cell.ids=c("WT","KO")) # Unable to allocate vector due to the size, I merged tumor clusters instead.

Tumor <- merge(Xenium_WT_tumor, y=Xenium_KO_tumor, add.cell.ids=c("WT","KO"))
Tumor$orig.ident
Tumor$celltype
table(Tumor$celltype)
head(colnames(Tumor))
tail(colnames(Tumor))

saveRDS(Tumor, "WT_KO_tumor.rds")

#########################################################################################

#Tumor = readRDS("WT_KO.rds")
table(Tumor$celltype)
table(Idents(Tumor))
head(colnames(Tumor))      
tail(colnames(Tumor))
unique(sapply(X = strsplit(colnames(Tumor), split = "_"), FUN = "[", 1))

head(Tumor@meta.data)
split_names = strsplit(colnames(Tumor), split = "_")
genotype <- sapply(split_names, '[',1)
head(genotype)
Tumor <- AddMetaData(Tumor, metadata = genotype, col.name = "genotype")
head(Tumor@meta.data)
tail(Tumor@meta.data)

combined_metadata <- paste(Tumor@meta.data$genotype, Tumor@meta.data$celltype, sep = "_")
Tumor <- AddMetaData(Tumor, metadata = combined_metadata, col.name = "genotype_celltype")
head(Tumor@meta.data)


###################### sctransform is left but no xenium in assay ###########################

Tumor <- RunPCA(Tumor, npcs = 30, features = rownames(Tumor))
Tumor <- RunUMAP(Tumor, dims = 1:30, n.neighbors = 5)
Tumor <- FindNeighbors(Tumor, reduction = "pca", dims = 1:30)
Tumor <- FindClusters(Tumor, resolution = 0.3)


DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "celltype",label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "genotype",label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap", label = F,label.size = 5)

Idents(Tumor) = "genotype_celltype"
table(Idents(Tumor))
table(Tumor$celltype)


#Tumor <- RenameIdents(Tumor, 
#                               `WT_Neutrophil_4` = "Neutrophil", `WT_Fibroblast_1` = "Fibroblast", 
#                               `WT_Neutrophil_1` = "Neutrophil",`WT_Fibroblast_2` = "Fibroblast",
#                               `WT_Neutrophil_2` = "Neutrophil", `WT_Macrophage_2` = "Macrophage", 
#                               `WT_Macrophage` = "Macrophage", `WT_Endothelial` = "Endothelial",
#                               `WT_Neutrophil_1` = "Neutrophil", `WT_Fibroblast_2` = "Fibroblast",
#                               
#                               `WT_Neutrophil_2` = "Neutrophil", `WT_Macrophage_2` = "Macrophage", 
#                               `WT_Fibroblast_3` = "Fibroblast",`WT_Neutrophil_3` = "Neutrophil",
#                               `KO_Others_not_clear` = "Others", `KO_Macrophage` = "Macrophage", 
#                               `KO_Endothelial` = "Endothelial", `KO_Fibroblast_1` = "Fibroblast"
#                               )

table(Idents(Tumor))
Tumor$genotype_celltype = Idents(Tumor)

my_levels <- c("WT_Epi_1","WT_Epi_2","WT_Epi_3","WT_Epi_4","WT_Epi_5","WT_Epi_6","WT_Epi_7","WT_Epi_8","WT_Epi_9","WT_Epi_10","WT_Epi_11","WT_Epi_12","WT_Epi_13","WT_Epi_14","WT_Epi_15",
               "KO_Epi_1","KO_Epi_2","KO_Epi_3","KO_Epi_4","KO_Epi_5","KO_Epi_6","KO_Epi_7","KO_Epi_8","KO_Epi_9","KO_Epi_10","KO_Epi_11","KO_Epi_12","KO_Epi_13"
               )
Tumor$genotype_celltype <- factor(x = Tumor$genotype_celltype, levels = my_levels)

color_list = c("#23BD3D","#5DDEB1","#096B08","#bef3c7","#33ff33","#59b300",
               "#ff6600","#ff704d","#ff6633","#cccc00","#cc9900","#ffff4d",
               "#ffc34d","#e6ac00","#ffbf80",
               "#003d99","#9999ff","#944dff","#3366ff","#3333ff","#66ccff",
               "#3399ff","#0066ff","#0000ff","#33ccff","#0099ff","#0066cc",
               "#0099cc")
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)


DotPlot(Tumor, features = c("Muc2","Spdef"), cols = "RdYlBu",scale.min = 0,scale.max = 60, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor, features = c("Muc2","Spdef"), cols = "RdYlBu",scale.min = 0,scale.max = 60, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")





saveRDS(Tumor, "WT_KO_Tumor_merged.rds")
#Tumor = readRDS("Tumor_merged.rds")



