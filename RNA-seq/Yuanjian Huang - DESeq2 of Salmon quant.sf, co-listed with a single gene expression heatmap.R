Sys.setenv(LANG = "en_US")
library("tximport")
library("readr")
library("tximportData")
library("DESeq2")
library("apeglm")
library("pheatmap")
library("vsn")
library("RColorBrewer")

# Read sample information
coldata <- read.csv("C:/Users/hyj/OneDrive/MDACC/P1. CRAD in mucinous CRC/GSEA/NCS RhoA treat PDO and PDC/2024_3_1_RNA-seq_PDC_C1144_NCS_vs_DMSO/Seq_info_PDC_C1144.csv", row.names=1)    # check Love - Analyzing RNA-seq data with DESeq2
coldata <- coldata[, "Group", drop = FALSE]
coldata$Group <- factor(coldata$Group)    # convert numeric variables to factors

# Read table linking transcripts to genes
tx2gene <- read.delim(file="J:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/tx2gene.tsv", header=FALSE)[, 1:3]

# Read gene matrix in the same sample order of Seq_info.csv
txi <- tximport(files=c("C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/Vehicle_REP1/quant.sf",
                        "C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/Vehicle_REP2/quant.sf",
                        "C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/Vehicle_REP3/quant.sf",
                        "C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/NCS_REP1/quant.sf",
                        "C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/NCS_REP2/quant.sf",
                        "C:/HYJ/Seq_data/2024_3_1_RNA-seq/Analysis/PDC_C1144/salmon/NCS_REP3/quant.sf"),
                type="salmon", tx2gene=tx2gene)

# Construct DESeqDataSet
dds <- DESeqDataSetFromTximport(txi=txi, colData=coldata, design=~Group)
assay(dds)

# Pre-filtering
smallestGroupSize <- 3
keep <- rowSums(counts(dds) >= 10) >= smallestGroupSize
dds <- dds[keep,]

# Extract filtered raw read counts
raw_counts <- assay(dds)
write.csv(x = raw_counts, file = "C:/Users/hyj/Documents/dds_raw_counts_filtered.csv", row.names = TRUE)

# Set comparison reference
dds$Group <- relevel(dds$Group, ref = "DMSO")

# Differential expression analysis
dds <- DESeq(dds)
norm_counts <- counts(dds, normalized = TRUE)    # Extract filtered and normalized read counts, for gene expression comparison
write.csv(x = norm_counts, file = "C:/Users/hyj/Documents/dds_normalized_counts_filtered.csv", row.names = TRUE)
res <- results(dds)
res
write.csv(res, file="C:/Users/hyj/Documents/res_PDC_C1144.csv", sep = ",", row.names=TRUE, col.names = TRUE)

# LFC estimates
resultsNames(dds)    # copy the comparison name to coef= in the next step
resLFC <- lfcShrink(dds, coef="Group_NCS_vs_DMSO", type="apeglm")
resLFC

# Data transformations
vsd <- vst(dds, blind=FALSE)    # vst just transforms all the counts in the matrix. The transformation will not be affected by how you set up the reference or the constrast

# Exporting vst matrix
vst_matrix=assay(vsd)
write.csv(vst_matrix, file="C:/Users/hyj/Documents/vst_matrix_PDC_C1144.csv", sep = ",", row.names=TRUE, col.names = TRUE)

# Heatmap of the count matrix
select <- order(rowMeans(counts(dds,normalized=TRUE)), decreasing=TRUE)[1:20]    # select the top 20 genes according to the mean normalized counts
df <- as.data.frame(colData(dds)[,"Group", drop = FALSE])
Matrix_heatmap=pheatmap(vst_matrix[select,],
                        cluster_rows=FALSE,
                        show_rownames=TRUE,
                        cluster_cols=FALSE,
                        annotation_col=df)

# Sample-to-sample distances
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(colnames(vsd))    # copy patient IDs
colnames(sampleDistMatrix) <- paste(colnames(vsd))
PuOr <- colorRampPalette( rev(brewer.pal(10, "PuOr")) )(1024)
Distance_heatmap=pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=PuOr,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 4,
         fontsize_col = 4,
         annotation_col=coldata)

# PCA plot
plotPCA(vsd, intgroup=c("Group"))

# SPDEF expression heatmap in the same order of Distance_heatmap
col_order <- Distance_heatmap$tree_col$order    # Retrieve the order of columns
col_labels_ordered <- Distance_heatmap$tree_col$labels[col_order]    # Retrieve corresponding column labels using the order
vst_matrix_reordered <- vst_matrix[, col_labels_ordered]    # Reorder the vst-transformed matrix according to col_labels_ordered
write.csv(vst_matrix_reordered, file="C:/Users/hyj/Documents/vst_matrix_reordered_PDC_C1144.csv", sep = ",", row.names=TRUE, col.names = TRUE)

gene_row <- vst_matrix_reordered["SPDEF", , drop = FALSE]    # Extract the row for SPDEF gene
gene_row <- as.data.frame(gene_row)

RdBu <- colorRampPalette( rev(brewer.pal(10, "RdBu")) )(1024)
Gene_heatmap=pheatmap(gene_row, 
         cluster_rows = FALSE, 
         cluster_cols = FALSE,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 4,
         fontsize_col = 4,
         col = RdBu)
