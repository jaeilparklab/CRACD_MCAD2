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

# Load the A Xenium data
path <- "X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_3__20250430__175355/output-XETG00320__0067739__Region_3__20250430__175355"
Xenium_A.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
Xenium_A.obj <- subset(Xenium_A.obj, subset = nCount_Xenium > 0)
VlnPlot(Xenium_A.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(Xenium_A.obj, fov = "fov", molecules = c("Agr2","Muc2", "Lyz", "Spdef", "Tff3"), nmols = 20000)


# Load the CA Xenium data
path <- "X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_5__20250430__175355/output-XETG00320__0067739__Region_5__20250430__175356"
Xenium_CA.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
Xenium_CA.obj <- subset(Xenium_CA.obj, subset = nCount_Xenium > 0)
VlnPlot(Xenium_CA.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(Xenium_CA.obj, fov = "fov", molecules = c("Agr2","Muc2", "Lyz", "Spdef", "Tff3"), nmols = 20000)


# Load the AOMDSS-WT Xenium data
path <- "X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_2__20250430__175355/output-XETG00320__0067739__Region_2__20250430__175355"
Xenium_WT.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
Xenium_WT.obj <- subset(Xenium_WT.obj, subset = nCount_Xenium > 0)
VlnPlot(Xenium_WT.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(Xenium_WT.obj, fov = "fov", molecules = c("Agr2","Muc2", "Lyz", "Spdef", "Tff3"), nmols = 20000)


# Load the AOMDSS-KO Xenium data
path <- "X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_1__20250430__175355 (2)/output-XETG00320__0067739__Region_1__20250430__175355"
Xenium_KO.obj <- LoadXenium(path, fov = "fov")
# remove cells with 0 counts
Xenium_KO.obj <- subset(Xenium_KO.obj, subset = nCount_Xenium > 0)
VlnPlot(Xenium_KO.obj, features = c("nFeature_Xenium", "nCount_Xenium"), ncol = 2, pt.size = 0)
ImageDimPlot(Xenium_KO.obj, fov = "fov", molecules = c("Agr2","Muc2", "Lyz", "Spdef", "Tff3"), nmols = 20000)







# Expression level
ImageFeaturePlot(Xenium_A.obj, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(25,5), size = 1.5, cols = c("white", "red"))                                                                                          
ImageFeaturePlot(Xenium_CA.obj, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(25,5), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(Xenium_WT.obj, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(25,5), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(Xenium_KO.obj, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(25,5), size = 1.5, cols = c("white", "red"))



# Zoomed image by crop-RKP
cropped.coords <- Crop(Xenium_A.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
Xenium_A.obj[["zoom"]] <- cropped.coords
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(Xenium_A.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(Xenium_A.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Muc2", "Spdef"), nmols = 10000)+
  theme(panel.grid = element_blank())

# Zoomed image by crop-EKP
cropped.coords_CA <- Crop(Xenium_CA.obj[["fov"]], x = c(1200, 2900), y = c(3750, 4550), coords = "plot")
Xenium_CA.obj[["zoom"]] <- cropped.coords_CA
# visualize cropped area with cell segmentations & selected molecules
DefaultBoundary(Xenium_CA.obj[["zoom"]]) <- "segmentation"
ImageDimPlot(Xenium_CA.obj, fov = "zoom", axes = TRUE, border.color = "white", border.size = 0.1, cols = "polychrome",
             coord.fixed = FALSE, molecules = c("Muc2", "Spdef"), nmols = 10000)+
  theme(panel.grid = element_blank())



# Normalization and clustering - A
Xenium_A.obj <- SCTransform(Xenium_A.obj, assay = "Xenium")
Xenium_A.obj <- RunPCA(Xenium_A.obj, npcs = 30, features = rownames(Xenium_A.obj))
Xenium_A.obj <- RunUMAP(Xenium_A.obj, dims = 1:30)
Xenium_A.obj <- FindNeighbors(Xenium_A.obj, reduction = "pca", dims = 1:30)
Xenium_A.obj <- FindClusters(Xenium_A.obj, resolution = 0.3)
DimPlot(Xenium_A.obj)


# Normalization and clustering - CA
Xenium_CA.obj <- SCTransform(Xenium_CA.obj, assay = "Xenium")
Xenium_CA.obj <- RunPCA(Xenium_CA.obj, npcs = 30, features = rownames(Xenium_CA.obj))
Xenium_CA.obj <- RunUMAP(Xenium_CA.obj, dims = 1:30)
Xenium_CA.obj <- FindNeighbors(Xenium_CA.obj, reduction = "pca", dims = 1:30)
Xenium_CA.obj <- FindClusters(Xenium_CA.obj, resolution = 0.3)
DimPlot(Xenium_CA.obj)


# Normalization and clustering - AOMDSS WT
Xenium_WT.obj <- SCTransform(Xenium_WT.obj, assay = "Xenium")
Xenium_WT.obj <- RunPCA(Xenium_WT.obj, npcs = 30, features = rownames(Xenium_CA.obj))
Xenium_WT.obj <- RunUMAP(Xenium_WT.obj, dims = 1:30)
Xenium_WT.obj <- FindNeighbors(Xenium_WT.obj, reduction = "pca", dims = 1:30)
Xenium_WT.obj <- FindClusters(Xenium_WT.obj, resolution = 0.3)
DimPlot(Xenium_WT.obj)


# Normalization and clustering - AOMDSS KO
Xenium_KO.obj <- SCTransform(Xenium_KO.obj, assay = "Xenium")
Xenium_KO.obj <- RunPCA(Xenium_KO.obj, npcs = 30, features = rownames(Xenium_CA.obj))
Xenium_KO.obj <- RunUMAP(Xenium_KO.obj, dims = 1:30)
Xenium_KO.obj <- FindNeighbors(Xenium_KO.obj, reduction = "pca", dims = 1:30)
Xenium_KO.obj <- FindClusters(Xenium_KO.obj, resolution = 0.3)
DimPlot(Xenium_KO.obj)



ImageDimPlot(Xenium_A.obj, fov = "fov",border.size = 0, axes = TRUE, , cols = "polychrome",
             coord.fixed = FALSE, nmols = 10000, group.by = 'seurat_clusters')+ theme (panel.grid = element_blank())
ImageDimPlot(Xenium_CA.obj, fov = "fov",border.size = 0, axes = TRUE, , cols = "polychrome",
             coord.fixed = FALSE, nmols = 10000, group.by = 'seurat_clusters')+theme (panel.grid = element_blank())


FeaturePlot(Xenium_A.obj, features = c("Muc2", "Spdef","Cdh1"),max.cutoff = 2)
FeaturePlot(Xenium_CA.obj, features = c("Muc2", "Spdef","Cdh1"),max.cutoff = 2)
FeaturePlot(Xenium_WT.obj, features = c("Muc2", "Spdef","Cdh1"),max.cutoff = 2)
FeaturePlot(Xenium_KO.obj, features = c("Muc2", "Spdef","Cdh1"),max.cutoff = 2)

saveRDS(Xenium_A.obj, "A_clustered.rds")
saveRDS(Xenium_CA.obj, "CA_clustered.rds")
saveRDS(Xenium_WT.obj, "WT_clustered.rds")
saveRDS(Xenium_KO.obj, "KO_clustered.rds")


# Importing the cluster information used in Xenium Explorer
cluster_A = read.csv("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_3__20250430__175355/output-XETG00320__0067739__Region_3__20250430__175355/analysis/clustering/gene_expression_graphclust/clusters.csv")
head(cluster_A)
Xenium_A.obj@meta.data[cluster_A$Barcode, "celltype"] = cluster_A$Cluster
head(Xenium_A.obj@meta.data)
str(Xenium_A.obj@meta.data)
DimPlot(Xenium_A.obj, group.by = "celltype")
# rename clusters by as it was annotated in Xenium explorer
length(Idents(Xenium_A.obj))
Idents(Xenium_A.obj) = "celltype"
table(Idents(Xenium_A.obj))
Xenium_A.obj <- RenameIdents(Xenium_A.obj, `1` = "Epi_1", `2` = "Epi_2", 
                             `3` = "Epi_3",`4` = "Fibroblast_1", `5` = "Epi_4", 
                             `6` = "Others_not_clear", 
                             `7` = "Macrophage", `8` = "Epi_5", 
                             `9` = "Epi_6", `10` = "Endothelial" )
Xenium_A.obj$celltype <- Idents(Xenium_A.obj)

saveRDS(Xenium_A.obj, "A_clustered.rds")
#Xenium_A.obj = readRDS("A_clustered.rds")


# Importing the cluster information used in Xenium Explorer
cluster_CA = read.csv("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_5__20250430__175355/output-XETG00320__0067739__Region_5__20250430__175356/analysis/clustering/gene_expression_graphclust/clusters.csv")
head(cluster_CA)
Xenium_CA.obj@meta.data[cluster_CA$Barcode, "celltype"] = cluster_CA$Cluster
head(Xenium_CA.obj@meta.data)
str(Xenium_CA.obj@meta.data)
DimPlot(Xenium_CA.obj, group.by = "celltype")
# rename clusters by as it was annotated in Xenium explorer
length(Idents(Xenium_CA.obj))
Idents(Xenium_CA.obj) = "celltype"
table(Idents(Xenium_CA.obj))
Xenium_CA.obj <- RenameIdents(Xenium_CA.obj, `1` = "Fibroblast_1", `2` = "Epi_1", `3` = "Epi_2",`4` = "Epi_3", `5` = "Epi_4", 
                             `6` = "Epi_5",`7` = "Macrophage", `8` = "Epi_6", `9` = "Epi_7", `10` = "Epi_8",
                             `11` = "Epi_9", `12` = "Neutrophil_1", `13` = "Endothelial",`14` = "Epi_10", `15` = "Epi_11", 
                             `16` = "T_cell",`17` = "Neutrophil_2", `18` = "Fibroblast_2", `19` = "B_cell", `20` = "Fibroblast_3",
                             `21` = "Epi_12", `22` = "Neutrophil_4", `23` = "Macrophage_2",`24` = "Neutrophil_3"
                             )
Xenium_CA.obj$celltype <- Idents(Xenium_CA.obj)

saveRDS(Xenium_CA.obj, "CA_clustered.rds")
#Xenium_CA.obj = readRDS("CA_clustered.rds")


# Importing the cluster information used in Xenium Explorer
cluster_WT = read.csv("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_2__20250430__175355/output-XETG00320__0067739__Region_2__20250430__175355/analysis/clustering/gene_expression_graphclust/clusters.csv")
head(cluster_WT)
Xenium_WT.obj@meta.data[cluster_WT$Barcode, "celltype"] = cluster_WT$Cluster
head(Xenium_WT.obj@meta.data)
str(Xenium_WT.obj@meta.data)
DimPlot(Xenium_WT.obj, group.by = "celltype")
# rename clusters by as it was annotated in Xenium explorer
length(Idents(Xenium_WT.obj))
Idents(Xenium_WT.obj) = "celltype"
table(Idents(Xenium_WT.obj))
Xenium_WT.obj <- RenameIdents(Xenium_WT.obj, `1` = "Fibroblast_1", `2` = "Macrophage", `3` = "Epi_1",`4` = "Epi_2", `5` = "Epi_3", 
                              `6` = "Epi_4",`7` = "B_cell_2", `8` = "Epi_5", `9` = "T_cell", `10` = "Fibroblast_2",
                              `11` = "Epi_6", `12` = "B_cell_1", `13` = "Epi_7",`14` = "Neutrophil_1", `15` = "Epi_8", 
                              `16` = "Fibroblast_3",`17` = "Endothelial", `18` = "Epi_9", `19` = "Epi_10", `20` = "Epi_11",
                              `21` = "Epi_12", `22` = "Epi_13", `23` = "Epi_14",`24` = "Epi_15",`25` = "Neutrophil_2"
)
Xenium_WT.obj$celltype <- Idents(Xenium_WT.obj)

saveRDS(Xenium_WT.obj, "WT_clustered.rds")
#Xenium_CA.obj = readRDS("WT_clustered.rds")


# Importing the cluster information used in Xenium Explorer
cluster_KO = read.csv("X:/Dept. Exp Rad Oncology/Lab Files/Lab Park/Xenium_05-05-2025/output-XETG00320__0067739__Region_1__20250430__175355 (2)/output-XETG00320__0067739__Region_1__20250430__175355/analysis/clustering/gene_expression_graphclust/clusters.csv")
head(cluster_KO)
Xenium_KO.obj@meta.data[cluster_KO$Barcode, "celltype"] = cluster_KO$Cluster
head(Xenium_KO.obj@meta.data)
str(Xenium_KO.obj@meta.data)
DimPlot(Xenium_KO.obj, group.by = "celltype")
# rename clusters by as it was annotated in Xenium explorer
length(Idents(Xenium_KO.obj))
Idents(Xenium_KO.obj) = "celltype"
table(Idents(Xenium_KO.obj))
Xenium_KO.obj <- RenameIdents(Xenium_KO.obj, `1` = "Epi_1", `2` = "Epi_2", `3` = "Epi_3",`4` = "Epi_4", `5` = "Epi_5", 
                              `6` = "B_cell_4",`7` = "Epi_6", `8` = "Fibroblast_1", `9` = "Epi_7", `10` = "T_cell",
                              `11` = "Endothelial_1", `12` = "B_cell_2", `13` = "Epi_8",`14` = "Fibroblast_3", `15` = "B_cell_3", 
                              `16` = "Epi_10",`17` = "Epi_11", `18` = "Epi_12", `19` = "Mucin_lake_1", `20` = "Fibroblast_2",
                              `21` = "Epi_13", `22` = "Mucin_lake_2", `23` = "Macrophage_1",`24` = "Neutrophil",`25` = "Epi_14",
                              `26` = "Fibroblast_4", `27` = "Fibroblast_5", `28` = "Endothelial_2")

Xenium_KO.obj <- RenameIdents(Xenium_KO.obj, `Epi_10` = "Epi_9", `Epi_11` = "Epi_10", `Epi_12` = "Epi_11",`Epi_13` = "Epi_12", `Epi_13` = "Epi_12", 
                              `Epi_14` = "Epi_13")
Xenium_KO.obj$celltype <- Idents(Xenium_KO.obj)
DimPlot(Xenium_KO.obj, group.by = "celltype")

saveRDS(Xenium_KO.obj, "KO_clustered.rds")
#Xenium_KO.obj = readRDS("KO_clustered.rds")



#### Tumor cell subset #####
A_tumor <- c("Epi_1","Epi_2","Epi_3", "Epi_4","Epi_5","Epi_6")
CA_tumor <- c("Epi_1","Epi_2","Epi_3", "Epi_4","Epi_5","Epi_6", "Epi_7","Epi_8","Epi_9", "Epi_10","Epi_11","Epi_12")
WT_tumor <- c("Epi_1","Epi_2","Epi_3", "Epi_4","Epi_5","Epi_6", "Epi_7","Epi_8","Epi_9", "Epi_10","Epi_11","Epi_12","Epi_13", "Epi_14","Epi_15")
KO_tumor <- c("Epi_1","Epi_2","Epi_3", "Epi_4","Epi_5","Epi_6", "Epi_7","Epi_8", "Epi_9","Epi_10","Epi_11","Epi_12", "Epi_13")

head(Xenium_A.obj)
Idents(Xenium_A.obj)

Xenium_A_tumor <- subset(Xenium_A.obj, idents = A_tumor)
Xenium_CA_tumor <- subset(Xenium_CA.obj, idents = CA_tumor)
Xenium_WT_tumor <- subset(Xenium_WT.obj, idents = WT_tumor)
Xenium_KO_tumor <- subset(Xenium_KO.obj, idents = KO_tumor)

DimPlot(Xenium_A.obj)
DimPlot(Xenium_A_tumor)

DimPlot(Xenium_CA.obj)
DimPlot(Xenium_CA_tumor)

DimPlot(Xenium_WT.obj)
DimPlot(Xenium_WT_tumor)

DimPlot(Xenium_KO.obj)
DimPlot(Xenium_KO_tumor)

DotPlot(Xenium_A_tumor, features = c("Muc2","Spdef"), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))
DotPlot(Xenium_CA_tumor, features = c("Muc2","Spdef"), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))
DotPlot(Xenium_WT_tumor, features = c("Muc2","Spdef"), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))
DotPlot(Xenium_KO_tumor, features = c("Muc2","Spdef"), dot.scale = 8, cols = "RdYlBu") +RotatedAxis()+ labs(y="group", x = "Genes") + theme(axis.text=element_text(size=13)) +theme(panel.grid.minor = element_blank(), panel.border = element_rect(colour = "black", fill=NA))



Xenium_A_tumor <- RenameIdents(Xenium_A_tumor, `Epi_1` = "A_Tumor_1", `Epi_2` = "A_Tumor_2",`Epi_3` = "A_Tumor_3",`Epi_4` = "A_Tumor_4", `Epi_5` = "A_Tumor_5",`Epi_6` = "A_Tumor_6")
Xenium_A_tumor$celltype_A <- Idents(Xenium_A_tumor)

Xenium_CA_tumor <- RenameIdents(Xenium_CA_tumor,`Epi_1` = "CA_Tumor_1", `Epi_2` = "CA_Tumor_2",`Epi_3` = "CA_Tumor_3",`Epi_4` = "CA_Tumor_4", `Epi_5` = "CA_Tumor_5",`Epi_6` = "CA_Tumor_6",
                                `Epi_7` = "CA_Tumor_7", `Epi_8` = "CA_Tumor_8",`Epi_9` = "CA_Tumor_9",`Epi_10` = "CA_Tumor_10", `Epi_11` = "CA_Tumor_11",`Epi_12` = "CA_Tumor_12")
Xenium_CA_tumor$celltype_CA <- Idents(Xenium_CA_tumor)




saveRDS(Xenium_CA_tumor, "CA_tumor.rds")
saveRDS(Xenium_A_tumor, "A_tumor.rds")

saveRDS(Xenium_WT_tumor, "WT_tumor.rds")
saveRDS(Xenium_KO_tumor, "KO_tumor.rds")


#Xenium_CA_tumor = readRDS("CA_tumor.rds")
#Xenium_A_tumor = readRDS("A_tumor.rds")






##########################################################################################################################
############################      Tumor cells of A and CA integration ####################################################            
############################  Merge can be done in other workstation by Jinho ########################################### 


options(future.globals.maxSize = 2000000000000)

Tumor <- merge(Xenium_CA_tumor, y=Xenium_A_tumor, add.cell.ids=c("CA","A"))
TUmor$orig.ident
Tumor$celltype
table(Tumor$celltype)
head(colnames(Tumor))
tail(colnames(Tumor))

A_CA <- merge(Xenium_CA.obj, y=Xenium_A.obj, add.cell.ids=c("CA","A"))
#########################################################################################

#A_CA = readRDS("A_CA.rds")
table(A_CA$celltype)
table(Idents(A_CA))
head(colnames(A_CA))      
tail(colnames(A_CA))
unique(sapply(X = strsplit(colnames(A_CA), split = "_"), FUN = "[", 1))

head(A_CA@meta.data)
split_names = strsplit(colnames(A_CA), split = "_")
genotype <- sapply(split_names, '[',1)
head(genotype)
A_CA <- AddMetaData(A_CA, metadata = genotype, col.name = "genotype")
head(A_CA@meta.data)
tail(A_CA@meta.data)

combined_metadata <- paste(A_CA@meta.data$genotype, A_CA@meta.data$celltype, sep = "_")
A_CA <- AddMetaData(A_CA, metadata = combined_metadata, col.name = "genotype_celltype")
head(A_CA@meta.data)


###################### sctransform is left but no xenium in assay ###########################

A_CA <- RunPCA(A_CA, npcs = 30, features = rownames(A_CA))
A_CA <- RunUMAP(A_CA, dims = 1:30)
A_CA <- FindNeighbors(A_CA, reduction = "pca", dims = 1:30)
A_CA <- FindClusters(A_CA, resolution = 0.3)


DimPlot(A_CA, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(A_CA, reduction = "umap",group.by = "celltype",label = F,label.size = 5)
DimPlot(A_CA, reduction = "umap",group.by = "genotype",label = F,label.size = 5)
DimPlot(A_CA, reduction = "umap", label = F,label.size = 5)

Idents(A_CA) = "genotype_celltype"
table(Idents(A_CA))
table(A_CA$celltype)


A_CA <- RenameIdents(A_CA, 
                               `CA_Neutrophil_4` = "Neutrophil", `CA_Fibroblast_1` = "Fibroblast", 
                               `CA_Neutrophil_1` = "Neutrophil",`CA_Fibroblast_2` = "Fibroblast",
                               `CA_Neutrophil_2` = "Neutrophil", `CA_Macrophage_2` = "Macrophage", 
                               `CA_Macrophage` = "Macrophage", `CA_Endothelial` = "Endothelial",
                               `CA_Neutrophil_1` = "Neutrophil", `CA_Fibroblast_2` = "Fibroblast",
                               
                               `CA_Neutrophil_2` = "Neutrophil", `CA_Macrophage_2` = "Macrophage", 
                               `CA_Fibroblast_3` = "Fibroblast",`CA_Neutrophil_3` = "Neutrophil",
                               `A_Others_not_clear` = "Others", `A_Macrophage` = "Macrophage", 
                               `A_Endothelial` = "Endothelial", `A_Fibroblast_1` = "Fibroblast"
                              
                               )
A_CA <- RenameIdents(A_CA, 
                     `CA_B_cell` = "B cell", `CA_T_cell` = "T cell"
                     )
A_CA <- RenameIdents(A_CA, 
                     `B cell` = "B_cell", `T cell` = "T_cell"
)


table(Idents(A_CA))
A_CA$genotype_celltype = Idents(A_CA)

my_levels <- c("A_Epi_1","A_Epi_2","A_Epi_3","A_Epi_4","A_Epi_5","A_Epi_6",
               "CA_Epi_1","CA_Epi_2","CA_Epi_3","CA_Epi_4","CA_Epi_5","CA_Epi_6","CA_Epi_7","CA_Epi_8","CA_Epi_9","CA_Epi_10","CA_Epi_11","CA_Epi_12",
               "T cell","B cell","Neutrophil","Macrophage","Fibroblast","Endothelial","Others")
A_CA$genotype_celltype <- factor(x = A_CA$genotype_celltype, levels = my_levels)

color_list = c("#23BD3D","#5DDEB1","#096B08","#bef3c7","#33ff33","#59b300",
               "#2FBDEC","#346CF0","#4136EE","#4d94ff","#4d79ff","#6666ff","#0080ff","#003366","#99bbff","#003d99","#9999ff","#944dff",
               "#FF0408","#F2E827","#9F6668","#FF7400","#ff66cc","#A9A7A9","#838685")
DimPlot(A_CA, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)

saveRDS(A_CA, "A_CA_merged.rds")
#A_CA = readRDS("A_CA_merged.rds")

################### Subsetting Tumor and Immune cells because they are not separated well ####################
Tumor <- subset(A_CA, idents = c("A_Epi_1","A_Epi_2","A_Epi_3","A_Epi_4","A_Epi_5","A_Epi_6",
                                 "CA_Epi_1","CA_Epi_2","CA_Epi_3","CA_Epi_4","CA_Epi_5","CA_Epi_6","CA_Epi_7","CA_Epi_8","CA_Epi_9","CA_Epi_10","CA_Epi_11","CA_Epi_12"))
#Immn <- subset(A_CA, idents = c("T_cell","B_cell","Neutrophil","Macrophage","Fibroblast"))
#table(Immn$genotype_celltype)

Tumor <- RunPCA(Tumor, npcs = 30, features = rownames(Tumor))
Tumor <- RunUMAP(Tumor, dims = 1:30)
Tumor <- FindNeighbors(Tumor, reduction = "pca", dims = 1:30)
Tumor <- FindClusters(Tumor, resolution = 0.3)
DimPlot(Tumor, reduction = "umap", label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype", label = T,label.size = 5)

#table(Tumor$genotype_celltype)
#my_levels <- c("EKP_Tumor_1","EKP_Tumor_2","EKP_Tumor_3",
#               "RKP_Tumor_1","RKP_Tumor_2","RKP_Tumor_3",
#               "T cell","Myeloid","Fibroblast","Endothelial","Others_1","Others_2")
#Tumor$genotype_celltype <- factor(x = Tumor$genotype_celltype, levels = my_levels)

#color_list = c("#23BD3D","#5DDEB1","#096B08",
#               "#2FBDEC","#346CF0","#4136EE",
#               "#FF0408","#F2E827","#9F6668","#FF7400","#A9A7A9","#838685")
#DimPlot(Tumor, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)



saveRDS(Tumor, "A_CA_tumor.rds")
#Tumor = readRDS("A_CA_tumor.rds")

Immn <- RunPCA(Immn, npcs = 30, features = rownames(Tumor))
Immn <- RunUMAP(Immn, dims = 1:30)
Immn <- FindNeighbors(Immn, reduction = "pca", dims = 1:30)
Immn <- FindClusters(Immn, resolution = 0.3)
DimPlot(Immn, reduction = "umap", label = F,label.size = 5)
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype", label = F,label.size = 5)
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype", label = T,label.size = 5)

table(Idents(Immn))
table(Immn$genotype_celltype)
Idents(Immn)= Immn$genotype_celltype


color_list = c("#FFBD22","#FF39BC")
DimPlot(Immn, reduction = "umap",group.by = "genotype_celltype",label = F,label.size = 5, cols = color_list)





saveRDS(Immn, "EKP_RKP_immune.rds")
#Immn = readRDS("EKP_RKP_immune.rds")



######################################################################################
########################## AddModuleScore ###########################################

DefaultAssay(Tumor)

#wnt = list(c("Axin2","Tbx3"))
#tgfb = list(c("Serpine1","Smad7"))
#hh = list(c("Gli","Ptch1"))
hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))


APP_hippo <- AddModuleScore(object = Tumor, features = hippo, name = "hippo", assay = "SCT", nbin=4)
head(APP_hippo)
DotPlot(APP_hippo, features = "hippo1", cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
VlnPlot(APP_hippo, features = "hippo1")

APP_secret <- AddModuleScore(object = Tumor, features = secret, name = "secret", assay = "SCT", nbin=4)
head(APP_secret)
DotPlot(APP_secret, features = "secret1", cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")


APP_notch <- AddModuleScore(object = Tumor, features = notch, name = "notch", assay = "SCT", nbin=4)
head(APP_notch)
DotPlot(APP_notch, features = "notch1", cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")




A_tumor$hippo1 <- APP_hippo$hippo1
A_tumor$secret1 <- APP_secret$secret1
A_tumor$notch1 <- APP_notch$notch1
CA_tumor$hippo1 <- APP_hippo$hippo1
CA_tumor$secret1 <- APP_secret$secret1
CA_tumor$notch1 <- APP_notch$notch1





ImageFeaturePlot(A_tumor, fov="fov.2", features = c("hippo1", "secret1", "notch1"), max.cutoff = c(0.15,0.8,0.05), min.cutoff = -0.5, size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(CA_tumor, fov="fov", features = c("hippo1", "secret1", "notch1"), max.cutoff = c(0.15,0.8,0.05), min.cutoff = -0.5, size = 1.5, cols = c("white", "red"))

p <- ImageFeaturePlot(
  object = A_tumor,
  features = "hippo1",
  min.cutoff = 0,
  max.cutoff = 0.5,size=1,
  combine = FALSE
)[[1]]
p + theme(
  panel.background = element_rect(fill = "black", color = NA),
  plot.background = element_rect(fill = "black", color = NA),  panel.grid = element_blank(),
  axis.ticks = element_blank(),
  axis.text = element_blank(),
  axis.title = element_blank(),
)

####################################################################################################################################################################################


head(Tumor@meta.data)
Idents(Tumor) <- "genotype_celltype"
table(Idents(Tumor))


A_tumor <- subset(x= Tumor, idents = c("A_Epi_1","A_Epi_2","A_Epi_3","A_Epi_4","A_Epi_5","A_Epi_6"))
CA_tumor <- subset(x= Tumor, idents = c("CA_Epi_1","CA_Epi_2","CA_Epi_3","CA_Epi_4","CA_Epi_5","CA_Epi_6","CA_Epi_7","CA_Epi_8","CA_Epi_9","CA_Epi_10","CA_Epi_11","CA_Epi_12"))

A_tumor@images
CA_tumor@images
Tumor@images

ImageFeaturePlot(A_tumor, fov="fov.2", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(CA_tumor, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))



ImageDimPlot(A_tumor, fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(CA_tumor, fov="fov", size=1.0, cols = "polychrome", axes=TRUE)

Idents(A_tumor) <- "genotype_celltype"
Idents(CA_tumor) <- "genotype_celltype"

A <- ImageDimPlot(
  object = A_tumor, size=1.5,
  group.by = "genotype_celltype",    # or whatever metadata you're plotting
  combine = FALSE
)[[1]]  # if combine = FALSE returns a list

# Step 2: Remove grid and axis elements
A + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank()
)

CA <- ImageDimPlot(
  object = CA_tumor, size=1.5,
  group.by = "genotype_celltype",    # or whatever metadata you're plotting
  combine = FALSE
)[[1]]  # if combine = FALSE returns a list

# Step 2: Remove grid and axis elements
CA + theme(
  panel.background = element_rect(fill = "white", color = NA),
  plot.background = element_rect(fill = "white", color = NA),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  axis.title = element_blank(),
  axis.line = element_blank()
)




############################### pValue (t-test) ####################################
#################### Tumor cells ############################
A_tumor.mat = A_tumor@assays$SCT@counts
A.cells = colnames(A_tumor@assays$SCT@counts)
A.genes = rownames(A_tumor@assays$SCT@counts)
rownames(A_tumor.mat) = A.genes
colnames(A_tumor.mat) = paste0("A_", A.cells)

CA_tumor.mat = CA_tumor@assays$SCT@counts
CA.cells = colnames(CA_tumor@assays$SCT@counts)
CA.genes = rownames(CA_tumor@assays$SCT@counts)
rownames(CA_tumor.mat) = CA.genes
colnames(CA_tumor.mat) = paste0("CA_", CA.cells)

output.df = data.frame()
for(i in 1:length(CA.genes)) {
  cat(i, "\t", length(CA.genes), "\n")
  tmp.gene = CA.genes[i]
  A.array = as.numeric(A_tumor.mat[tmp.gene,])
  CA.array = as.numeric(CA_tumor.mat[tmp.gene,])
  test.result = t.test(A.array, CA.array)
  pval = test.result$p.value
  output.df[i,1] = tmp.gene
  output.df[i,2] = pval
  
  if(mean(A.array) > mean(CA.array)) {
    output.df[i,3] = "A"
  }
  if(mean(A.array) < mean(CA.array)) {
    output.df[i,3] = "CA"
  }
}

colnames(output.df) = c("Gene", "pval", "Group")
output.df[which(output.df$Gene == "Muc2"),]
output.df[which(output.df$Gene == "Spdef"),]



########################## p Values of scores ##########################################################################

# Step 1: Extract module scores from meta.data
A_scores <- A_tumor@meta.data[, c("hippo1", "secret1", "notch1")]
CA_scores <- CA_tumor@meta.data[, c("hippo1", "secret1", "notch1")]

# Step 2: Add row labels for tracking
rownames(A_scores) <- paste0("A_", rownames(A_scores))
rownames(CA_scores) <- paste0("CA_", rownames(CA_scores))

# Step 3: Merge the data (if needed later for plotting)
all_scores <- rbind(A_scores, CA_scores)
group <- c(rep("A", nrow(A_scores)), rep("CA", nrow(CA_scores)))

# Step 4: Run t-tests and collect means
features <- colnames(A_scores)
output.df <- data.frame()

for (i in seq_along(features)) {
  feature <- features[i]
  A_vals <- A_scores[[feature]]
  CA_vals <- CA_scores[[feature]]
  
  test.result <- t.test(A_vals, CA_vals)
  pval <- test.result$p.value
  A_mean <- mean(A_vals)
  CA_mean <- mean(CA_vals)
  
  output.df[i, "Feature"] <- feature
  output.df[i, "pval"] <- pval
  output.df[i, "A_mean"] <- A_mean
  output.df[i, "CA_mean"] <- CA_mean
  output.df[i, "Higher_in"] <- ifelse(A_mean > CA_mean, "A", "CA")
}

# View output
print(output.df)


########################## Vioin plots of scores by genotype ##########################################################################

Tumor$hippo1 <- APP_hippo$hippo1
Tumor$secret1 <- APP_secret$secret1
Tumor$notch1 <- APP_notch$notch1

VlnPlot(Tumor, features = c("hippo1",'secret1',"notch1"), group.by = "genotype", split.by = "genotype", ncol = 2, pt.size = 0,
        cols = c("coral", "deepskyblue2"))
