library(SoupX)
library(Seurat)
library(DropletUtils)

#### GSE236698, 12 samples
# Set seed for reproducibility of random processes
set.seed(1)

# Define sample IDs to be processed
samples <- c('P1_N', 'P2_N', 'P3_N', 'P4_N', 'P5_N', 'P6_N',
             'P1_T', 'P2_T', 'P3_T', 'P4_T', 'P5_T', 'P6_T')

# Base directory for storing scRNA-seq data
base_dir <- "E:/GSE236698/Matrix/scRNA-seq/"

# Main processing loop for each sample
for (sample in samples) {
  cat(paste0("\n==== Processing sample: ", sample, " ====\n"))
  
  # Construct directory paths for input and output data
  raw_matrix_dir <- file.path(base_dir, sample, "outs", "raw_feature_bc_matrix")
  filtered_matrix_dir <- file.path(base_dir, sample, "outs", "filtered_feature_bc_matrix")
  output_dir <- file.path(base_dir, sample, "outs", "filtered_feature_bc_matrix_SoupX")
  
  # Validate input directories exist
  if (!dir.exists(raw_matrix_dir)) {
    cat(paste0("Warning: Raw matrix directory not found: ", raw_matrix_dir, "\n"))
    next
  }
  
  if (!dir.exists(filtered_matrix_dir)) {
    cat(paste0("Warning: Filtered matrix directory not found: ", filtered_matrix_dir, "\n"))
    next
  }
  
  # Process sample with error handling
  tryCatch({
    # Read 10X Genomics count matrices
    tod <- Read10X(data.dir = raw_matrix_dir)    # Raw counts (all barcodes)
    toc <- Read10X(data.dir = filtered_matrix_dir)  # Filtered counts (cell barcodes)
    
    # Check and align gene sets between raw and filtered matrices
    if (!identical(rownames(tod), rownames(toc))) {
      cat("Warning: Gene lists do not match. Attempting to align...\n")
      
      # Find common genes between both matrices
      common_genes <- intersect(rownames(tod), rownames(toc))
      
      if (length(common_genes) == 0) {
        stop("No common genes found between raw and filtered matrices!")
      }
      
      # Subset matrices to common genes and ensure identical order
      tod <- tod[common_genes, ]
      toc <- toc[common_genes, ]
      
      cat(paste0("Aligned ", length(common_genes), " common genes\n"))
    }
    
    # Create SoupChannel object for ambient RNA correction
    sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
    
    # Create Seurat object and perform standard workflow
    srat <- CreateSeuratObject(counts = sc$toc)
    
    # Normalize data (required before scaling)
    srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
    
    # Identify highly variable features and scale data
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
    srat <- ScaleData(srat)
    
    # Dimensionality reduction and clustering
    srat <- RunPCA(srat, npcs = 30, verbose = FALSE)
    srat <- RunUMAP(srat, reduction = "pca", dims = 1:30)
    srat <- FindNeighbors(srat, reduction = "pca", dims = 1:30)
    srat <- FindClusters(srat, resolution = 1)
    
    # Extract UMAP coordinates and cluster assignments
    metadata <- as.data.frame(srat[["umap"]]@cell.embeddings)
    colnames(metadata) <- c("UMAP_1", "UMAP_2")  # More descriptive names
    metadata$Cluster <- factor(srat@meta.data[rownames(metadata), "seurat_clusters"])
    
    # Integrate Seurat clustering results into SoupX
    sc <- setClusters(sc, setNames(metadata$Cluster, rownames(metadata)))
    sc <- setDR(sc, metadata[colnames(sc$toc), c("UMAP_1", "UMAP_2")])
    
    # Estimate ambient RNA contamination and adjust counts
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc, roundToInt = TRUE)
    
    # Write corrected counts in 10X HDF5 format
    DropletUtils:::write10xCounts(output_dir, x = out)
    
    cat(paste0("Sample ", sample, " processed successfully!\n"))
    
  }, error = function(e) {
    cat(paste0("Failed to process sample ", sample, ": ", e$message, "\n"))
  })
}

cat("\nBatch processing completed for all samples!\n")



#### SMC cohort, 33 samples
# Set seed for reproducibility of random processes
set.seed(1)

# Define sample IDs to be processed
samples <- c('SMC01N', 'SMC02N', 'SMC03N', 'SMC04N', 'SMC05N', 'SMC06N', 'SMC07N', 'SMC08N', 'SMC09N', 'SMC10N',
             'SMC01T', 'SMC02T', 'SMC03T', 'SMC04T', 'SMC05T', 'SMC06T', 'SMC07T', 'SMC08T', 'SMC09T', 'SMC10T', 'SMC11T', 'SMC14T', 'SMC15T', 'SMC16T', 'SMC18T', 'SMC19T', 'SMC21T', 'SMC22T', 'SMC23T', 'SMC24T', 'SMC25T',
             'SMC17T', 'SMC20T')

# Base directory for storing scRNA-seq data
base_dir <- "H:/EGAD00001005191_EGAD00001005198/Matrix_v2024/scRNA-seq/EGAD00001005198/"

# Main processing loop for each sample
for (sample in samples) {
  cat(paste0("\n==== Processing sample: ", sample, " ====\n"))
  
  # Construct directory paths for input and output data
  raw_matrix_dir <- file.path(base_dir, sample, "outs", "raw_feature_bc_matrix")
  filtered_matrix_dir <- file.path(base_dir, sample, "outs", "filtered_feature_bc_matrix")
  output_dir <- file.path(base_dir, sample, "outs", "filtered_feature_bc_matrix_SoupX")
  
  # Validate input directories exist
  if (!dir.exists(raw_matrix_dir)) {
    cat(paste0("Warning: Raw matrix directory not found: ", raw_matrix_dir, "\n"))
    next
  }
  
  if (!dir.exists(filtered_matrix_dir)) {
    cat(paste0("Warning: Filtered matrix directory not found: ", filtered_matrix_dir, "\n"))
    next
  }
  
  # Process sample with error handling
  tryCatch({
    # Read 10X Genomics count matrices
    tod <- Read10X(data.dir = raw_matrix_dir)    # Raw counts (all barcodes)
    toc <- Read10X(data.dir = filtered_matrix_dir)  # Filtered counts (cell barcodes)
    
    # Check and align gene sets between raw and filtered matrices
    if (!identical(rownames(tod), rownames(toc))) {
      cat("Warning: Gene lists do not match. Attempting to align...\n")
      
      # Find common genes between both matrices
      common_genes <- intersect(rownames(tod), rownames(toc))
      
      if (length(common_genes) == 0) {
        stop("No common genes found between raw and filtered matrices!")
      }
      
      # Subset matrices to common genes and ensure identical order
      tod <- tod[common_genes, ]
      toc <- toc[common_genes, ]
      
      cat(paste0("Aligned ", length(common_genes), " common genes\n"))
    }
    
    # Create SoupChannel object for ambient RNA correction
    sc <- SoupChannel(tod, toc, calcSoupProfile = TRUE)
    
    # Create Seurat object and perform standard workflow
    srat <- CreateSeuratObject(counts = sc$toc)
    
    # Normalize data (required before scaling)
    srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
    
    # Identify highly variable features and scale data
    srat <- FindVariableFeatures(srat, selection.method = "vst", nfeatures = 2000)
    srat <- ScaleData(srat)
    
    # Dimensionality reduction and clustering
    srat <- RunPCA(srat, npcs = 30, verbose = FALSE)
    srat <- RunUMAP(srat, reduction = "pca", dims = 1:30)
    srat <- FindNeighbors(srat, reduction = "pca", dims = 1:30)
    srat <- FindClusters(srat, resolution = 1)
    
    # Extract UMAP coordinates and cluster assignments
    metadata <- as.data.frame(srat[["umap"]]@cell.embeddings)
    colnames(metadata) <- c("UMAP_1", "UMAP_2")  # More descriptive names
    metadata$Cluster <- factor(srat@meta.data[rownames(metadata), "seurat_clusters"])
    
    # Integrate Seurat clustering results into SoupX
    sc <- setClusters(sc, setNames(metadata$Cluster, rownames(metadata)))
    sc <- setDR(sc, metadata[colnames(sc$toc), c("UMAP_1", "UMAP_2")])
    
    # Estimate ambient RNA contamination and adjust counts
    sc <- autoEstCont(sc)
    out <- adjustCounts(sc, roundToInt = TRUE)
    
    # Write corrected counts in 10X HDF5 format
    DropletUtils:::write10xCounts(output_dir, x = out)
    
    cat(paste0("Sample ", sample, " processed successfully!\n"))
    
  }, error = function(e) {
    cat(paste0("Failed to process sample ", sample, ": ", e$message, "\n"))
  })
}

cat("\nBatch processing completed for all samples!\n")
