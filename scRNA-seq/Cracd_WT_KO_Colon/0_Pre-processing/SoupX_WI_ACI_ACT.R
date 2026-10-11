library(SoupX)
library(Seurat)

tod <- Read10X(data.dir = "D:/BaiduSyncdisk/HYJ/Seq_data/2021_6_25_scRNA-seq/WT_Colon/raw_feature_bc_matrix/")
toc <- Read10X(data.dir = "D:/BaiduSyncdisk/HYJ/Seq_data/2021_6_25_scRNA-seq/WT_Colon/filtered_feature_bc_matrix/")
sc = SoupChannel(tod, toc, calcSoupProfile=TRUE)

set.seed(1)
srat = CreateSeuratObject(sc$toc)
srat = ScaleData(srat)
srat = FindVariableFeatures(srat)
srat = RunPCA(srat,pcs.compute=30)
srat = RunUMAP(srat, dims=seq(30))
srat = FindNeighbors(srat, reductions.type="umap", dims = seq(30), resolution = 1)
srat = FindClusters(srat,reductions.type="umap", dims = seq(30), resolution = 1)

metadata = as.data.frame(srat[['umap']]@cell.embeddings)
colnames(metadata) = c('RD1','RD2')
metadata$Cluster = factor(srat@meta.data[rownames(metadata),'RNA_snn_res.1'])
sc = setClusters(sc, setNames(metadata$Cluster, rownames(metadata)))
sc = setDR(sc, metadata[colnames(sc$toc), c("RD1", "RD2")])

sc = autoEstCont(sc)
out = adjustCounts(sc, roundToInt = TRUE)    # f roundToInt = TRUE, will produce an output count matrix that consists only integers.
# if roundToInt = FALSE, will result in non-integer counts as representing probabilities that a count occurred (similar to the output from pseudo-aligners), making scVelo in turn interprets as the data having been normalized.
DropletUtils:::write10xCounts("D:/BaiduSyncdisk/HYJ/Seq_data/2021_6_25_scRNA-seq/WT_Colon/filtered_feature_bc_matrix_SoupX_out/", x=out)
