library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

library(reticulate)
use_python("C:/Users/hyj/Documents/.virtualenvs/r-reticulate/Scripts/python.exe")    # numpy and anndata were installed in virtual env. Use py_config() to check pythonhome path, should also be the virtual env. 
adata <- import("anndata", convert = FALSE)
adata_object <- adata$read_h5ad('C:/Users/hyj/Downloads/human_NMCAD_cellchat.h5ad')
# access normalized data matrix
data.input <- t(py_to_r(adata_object$X))
rownames(data.input) <- rownames(py_to_r(adata_object$var))
colnames(data.input) <- rownames(py_to_r(adata_object$obs))
# access meta data
meta.data <- py_to_r(adata_object$obs)
meta <- meta.data

library(Matrix)
data.input <- as(data.input, "CsparseMatrix")

# Create a CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "cluster")
levels(cellchat@idents)    # show the cell labels
# [1] "B"           "DC"          "Endothelial" "Fibroblasts" "Macrophages" "Mast"        "Monocytes"   "NK"          "Neurons"     "Neutrophils"
# [11] "Parenchymal" "Plasma"      "T"

cellchat@DB <- CellChatDB.human    # use all datasets
cellchat <- subsetData(cellchat)

cellchat <- identifyOverExpressedGenes(object=cellchat, thresh.pc = 0, thresh.fc = 0.322, thresh.p = 0.05)    # https://github.com/sqjin/CellChat/issues/590
cellchat <- identifyOverExpressedInteractions(object=cellchat)

cellchat <- computeCommunProb(object=cellchat, type = c("triMean"), trim = 0.25, LR.use = NULL, raw.use = TRUE, population.size = FALSE)
cellchat <- filterCommunication(object=cellchat, min.cells = 25)

cellchat <- computeCommunProbPathway(object=cellchat, thresh = 0.05)
cellchat <- aggregateNet(object=cellchat, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, thresh = 0.05)

cluster_colors=c('#1f77b4',
                 '#aec7e8',
                 '#ffbb78',
                 '#98df8a',
                 '#d62728',
                 '#9467bd',
                 '#8c564b',
                 '#c49c94',
                 '#f7b6d2',
                 '#c7c7c7',
                 '#bcbd22',
                 '#17becf',
                 '#9edae5')

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, color.use =cluster_colors, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, color.use =cluster_colors, vertex.weight = groupSize, weight.scale = TRUE, label.edge= FALSE, title.name = "Interaction weights/strength")

mat <- cellchat@net$weight
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

netVisual_aggregate(object=cellchat, signaling= c("NOTCH"), color.use =cluster_colors, thresh = 0.05, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("hierarchy"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
netVisual_aggregate(object=cellchat, signaling= c("NOTCH"), color.use =cluster_colors, thresh = 0.05, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("circle"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
netVisual_aggregate(object=cellchat, signaling= c("NOTCH"), color.use =cluster_colors, thresh = 0.05, vertex.receiver = seq(1,13), sources.use = c(3,4,10,11), targets.use = c(11), idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("circle"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
netVisual_aggregate(object=cellchat, signaling= c("NOTCH"), color.use =cluster_colors, thresh = 0.05, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("chord"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
netVisual_aggregate(object=cellchat, signaling= c("NOTCH"), color.use =cluster_colors, thresh = 0.05, vertex.receiver = seq(1,13), sources.use = c(3,4,10,11), targets.use = c(11), idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("chord"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)

netAnalysis_contribution(object=cellchat, signaling = c("NOTCH"), thresh = 0.05, font.size = 8, font.size.title = 8)
pairLR.NOTCH <- extractEnrichedLR(object=cellchat, signaling = c("NOTCH"), geneLR.return = FALSE, enriched.only = TRUE, thresh = 0.05)
pairLR.NOTCH
LR.show <- pairLR.NOTCH[5,] # choose JAG1_NOTCH1
LR.show
netVisual_individual(object=cellchat, signaling= c("NOTCH"), pairLR.use = LR.show, color.use =cluster_colors, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("hierarchy"), thresh = 0.05)
netVisual_individual(object=cellchat, signaling= c("NOTCH"), pairLR.use = LR.show, color.use =cluster_colors, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("circle"), thresh = 0.05)
netVisual_individual(object=cellchat, signaling= c("NOTCH"), pairLR.use = LR.show, color.use =cluster_colors, vertex.receiver = seq(1,13), sources.use = c(3,4,10,11), targets.use = c(11), top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("circle"), thresh = 0.05)
netVisual_individual(object=cellchat, signaling= c("NOTCH"), pairLR.use = LR.show, color.use =cluster_colors, vertex.receiver = seq(1,13), sources.use = NULL, targets.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("chord"), thresh = 0.05)
netVisual_individual(object=cellchat, signaling= c("NOTCH"), pairLR.use = LR.show, color.use =cluster_colors, vertex.receiver = seq(1,13), sources.use = c(3,4,10,11), targets.use = c(11), top = 1, remove.isolate = FALSE, vertex.weight = NULL, weight.scale = TRUE, layout = c("chord"), thresh = 0.05)

netVisual_bubble(object=cellchat, sources.use = c(3,4,10,11), targets.use = c(11), signaling = NULL, pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = NULL, remove.isolate = FALSE)
netVisual_bubble(object=cellchat, sources.use = c(3,4,10,11), targets.use = c(11), signaling = c("NOTCH"), pairLR.use = NULL, color.heatmap = c("viridis"), thresh = 0.05, comparison = NULL, remove.isolate = FALSE)
pairLR.NOTCH <- extractEnrichedLR(object=cellchat, signaling = c("NOTCH"), geneLR.return = FALSE, enriched.only = TRUE, thresh = 0.05)
pairLR.NOTCH
netVisual_bubble(object=cellchat, sources.use = c(3,4,10,11), targets.use = c(11), signaling = NULL, pairLR.use = pairLR.NOTCH, color.heatmap = c("viridis"), thresh = 0.05, comparison = NULL, remove.isolate = FALSE)
netVisual_chord_gene(object=cellchat, slot.name = "net", color.use =cluster_colors, signaling = NULL, pairLR.use = NULL, sources.use = c(3,4,10,11), targets.use = c(11), lab.cex = 0.5, thresh = 0.05)
netVisual_chord_gene(object=cellchat, slot.name = "net", color.use =cluster_colors, signaling = c("NOTCH"), pairLR.use = NULL, sources.use = c(3,4,10,11), targets.use = c(11), lab.cex = 0.5, thresh = 0.05)
netVisual_chord_gene(object=cellchat, slot.name = "netP", color.use =cluster_colors, signaling = NULL, pairLR.use = NULL, sources.use = c(3,4,10,11), targets.use = c(11), lab.cex = 0.5, thresh = 0.05)
netVisual_chord_gene(object=cellchat, slot.name = "netP", color.use =cluster_colors, signaling = c("NOTCH"), pairLR.use = NULL, sources.use = c(3,4,10,11), targets.use = c(11), lab.cex = 0.5, thresh = 0.05)

plotGeneExpression(object=cellchat, signaling = c("NOTCH"), enriched.only = TRUE, type = c("violin"), color.use =cluster_colors, group.by = NULL)

cellchat <- netAnalysis_computeCentrality(object =cellchat, slot.name = "netP", thresh = 0.05)
netAnalysis_signalingRole_network(object =cellchat, signaling = c("NOTCH"), slot.name = "netP", measure = c("outdeg", "indeg", "flowbet", "info"), measure.name = c("Sender", "Receiver", "Mediator", "Influencer"), color.use =cluster_colors, color.heatmap = "BuGn", width = 8, height = 2.5, font.size = 8)
netAnalysis_signalingRole_scatter(object=cellchat, signaling = NULL, color.use =cluster_colors, slot.name = "netP")
netAnalysis_signalingRole_scatter(object=cellchat, signaling = c("NOTCH"), color.use =cluster_colors, slot.name = "netP")

ht1=netAnalysis_signalingRole_heatmap(object=cellchat, signaling = NULL, pattern = c("outgoing"), slot.name = "netP", color.use =cluster_colors, color.heatmap = "BuGn", font.size = 8, font.size.title = 8)
ht2=netAnalysis_signalingRole_heatmap(object=cellchat, signaling = NULL, pattern = c("incoming"), slot.name = "netP", color.use =cluster_colors, color.heatmap = "BuGn", font.size = 8, font.size.title = 8)
ht3=netAnalysis_signalingRole_heatmap(object=cellchat, signaling = NULL, pattern = c("all"), slot.name = "netP", color.use =cluster_colors, color.heatmap = "BuGn", font.size = 8, font.size.title = 8)
ht1+ht2+ht3

library(NMF)
library(ggalluvial)
selectK(object=cellchat, slot.name = "netP", pattern = "outgoing")
nPatterns = 5
cellchat <- identifyCommunicationPatterns(object=cellchat, slot.name = "netP", pattern = "outgoing", k = nPatterns, color.use =cluster_colors, color.heatmap = "Spectral", font.size = 8)
netAnalysis_river(object=cellchat, slot.name = "netP", pattern = "outgoing", cutoff = 0.5, sources.use = NULL, targets.use = NULL, signaling = NULL, color.use =cluster_colors, font.size = 2.5)
netAnalysis_dot(object=cellchat, slot.name = "netP", pattern = "outgoing", cutoff = NULL, color.use =cluster_colors, font.size = 8)

selectK(object=cellchat, slot.name = "netP", pattern = "incoming")
nPatterns = 4
cellchat <- identifyCommunicationPatterns(object=cellchat, slot.name = "netP", pattern = "incoming", k = nPatterns, color.use =cluster_colors, color.heatmap = "Spectral", font.size = 8)
netAnalysis_river(object=cellchat, slot.name = "netP", pattern = "incoming", cutoff = 0.5, sources.use = NULL, targets.use = NULL, signaling = NULL, color.use =cluster_colors, font.size = 2.5)
netAnalysis_dot(object=cellchat, slot.name = "netP", pattern = "incoming", cutoff = NULL, color.use =cluster_colors, font.size = 8)

cellchat <- computeNetSimilarity(object=cellchat, slot.name = "netP", type = "functional", thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "functional", umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "functional", do.parallel = TRUE, nCores = 32)
netVisual_embedding(object=cellchat, slot.name = "netP", type = "functional", color.use =cluster_colors, label.size = 3, font.size = 8, font.size.title = 8)

cellchat <- computeNetSimilarity(object=cellchat, slot.name = "netP", type = "structural", thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "structural", umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "structural", do.parallel = TRUE, nCores = 32)
netVisual_embedding(object=cellchat, slot.name = "netP", type = "structural", color.use =cluster_colors, label.size = 3, font.size = 8, font.size.title = 8)

saveRDS(cellchat, file = "C:/Users/hyj/Documents/human_NMCAD_cellchat_0322.rds")
cellchat <- readRDS(file = "C:/Users/hyj/Documents/human_NMCAD_cellchat_0322.rds")
