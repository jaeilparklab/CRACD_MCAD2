library(CellChat)
library(patchwork)
library(ComplexHeatmap)

cellchat_WI <- readRDS(file = "C:/Users/hyj/Documents/WI_cellchat_0322.rds")
cellchat_ACI <- readRDS(file = "C:/Users/hyj/Documents/ACI_cellchat_0322.rds")
object.list <- list(WI = cellchat_WI, ACI = cellchat_ACI)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

leiden_colors=c('#ffa500ff', '#ee82eeff', '#9acd32ff', '#ff69b4ff', '#ff1493ff', '#6495edff', '#4169e1ff', '#0000ffff', '#00008bff', '#191970ff', '#808080ff', '#000000ff', '#228b22ff', '#3cb371ff', '#800080ff', '#008000ff')

# [1] "Stem"            "Progenitors_1"   "Progenitors_2"   "Progenitors_3"   "Progenitors_4"   "Colonocytes_1"   "Colonocytes_2"  
# [8] "Colonocytes_3"   "Colonocytes_4"   "Colonocytes_5"   "Colonocytes_6"   "Colonocytes_7"   "Goblet_1"        "Goblet_2"       
# [15] "Enteroendocrine" "Tuft"

compareInteractions(object=cellchat, measure = c("count"), color.use =c('#0000ffff', '#008000ff'), group = c(1,2), show.legend = FALSE)
compareInteractions(object=cellchat, measure = c("weight"), color.use =c('#0000ffff', '#008000ff'), group = c(1,2), show.legend = FALSE)
netVisual_diffInteraction(object=cellchat, comparison = c(1, 2), measure = c("count"), color.use =leiden_colors, color.edge = c("Red", "Blue"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, weight.scale = TRUE)
netVisual_diffInteraction(object=cellchat, comparison = c(1, 2), measure = c("weight"), color.use =leiden_colors, color.edge = c("Red", "Blue"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, weight.scale = TRUE)
netVisual_heatmap(object=cellchat, comparison = c(1, 2), measure = c("count"), signaling = NULL, slot.name = c("netP"), color.use =leiden_colors, color.heatmap = c("Blue", "Red"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
netVisual_heatmap(object=cellchat, comparison = c(1, 2), measure = c("weight"), signaling = NULL, slot.name = c("netP"), color.use =leiden_colors, color.heatmap = c("Blue", "Red"), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link))
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], color.use =leiden_colors, title = names(object.list)[i], weight.MinMax = weight.MinMax)
}
patchwork::wrap_plots(plots = gg)

netAnalysis_signalingChanges_scatter(object=cellchat, idents.use = c('Progenitors_2'), color.use =c('Gray','#0000ffff', '#008000ff'), comparison = c(1, 2), signaling = NULL, slot.name = "netP", font.size = 8)
netAnalysis_signalingChanges_scatter(object=cellchat, idents.use = c('Goblet_1'), color.use =c('Gray','#0000ffff', '#008000ff'), comparison = c(1, 2), signaling = NULL, slot.name = "netP", font.size = 8)

cellchat <- computeNetSimilarityPairwise(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "functional", comparison = c(1, 2), nCores = 32)
netVisual_embeddingPairwise(object=cellchat, slot.name = "netP", type = "functional", color.use =leiden_colors, comparison = c(1, 2), dot.size = c(4, 8), label.size = 4)
netVisual_embeddingPairwiseZoomIn(object=cellchat, slot.name = "netP", type = "functional", color.use =leiden_colors, comparison = c(1, 2), nCol = 2, dot.size = c(4, 8), label.size = 4)

cellchat <- computeNetSimilarityPairwise(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), thresh = 0.05)
cellchat <- netEmbedding(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), umap.method = c("umap-learn"))
cellchat <- netClustering(object=cellchat, slot.name = "netP", type = "structural", comparison = c(1, 2), nCores = 32)
netVisual_embeddingPairwise(object=cellchat, slot.name = "netP", type = "structural", color.use =leiden_colors, comparison = c(1, 2), dot.size = c(4, 8), label.size = 4)
netVisual_embeddingPairwiseZoomIn(object=cellchat, slot.name = "netP", type = "structural", color.use =leiden_colors, comparison = c(1, 2), nCol = 2, dot.size = c(4, 8), label.size = 4)

rankSimilarity(object=cellchat, slot.name = "netP", type = "functional", comparison2 = c(1, 2), font.size = 8)

rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = TRUE, sources.use = NULL, targets.use = NULL, signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)
rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = FALSE, sources.use = NULL, targets.use = NULL, signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)

rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = TRUE, sources.use = NULL, targets.use = c(3), signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)
rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = TRUE, sources.use = c(3,13), targets.use = c(3), signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)

rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = TRUE, sources.use = c(3), targets.use = c(3), signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)
rankNet(object=cellchat, slot.name = "netP", measure = "weight", mode = "comparison", comparison = c(1, 2), color.use =c('#0000ffff', '#008000ff'), stacked = TRUE, sources.use = c(13), targets.use = c(3), signaling = NULL, do.stat = TRUE, cutoff.pvalue = 0.05, thresh = 0.05, font.size = 8, do.flip=FALSE, x.angle=90)

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object=object.list[[i]], signaling = pathway.union, pattern = "outgoing", slot.name = "netP", color.use =leiden_colors, color.heatmap = "BuGn", title = names(object.list)[i], font.size = 8, font.size.title = 8)
ht2 = netAnalysis_signalingRole_heatmap(object=object.list[[i+1]], signaling = pathway.union, pattern = "outgoing", slot.name = "netP", color.use =leiden_colors, color.heatmap = "BuGn", title = names(object.list)[i+1], font.size = 8, font.size.title = 8)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object=object.list[[i]], signaling = pathway.union, pattern = "incoming", slot.name = "netP", color.use =leiden_colors, color.heatmap = "BuGn", title = names(object.list)[i], font.size = 8, font.size.title = 8)
ht2 = netAnalysis_signalingRole_heatmap(object=object.list[[i+1]], signaling = pathway.union, pattern = "incoming", slot.name = "netP", color.use =leiden_colors, color.heatmap = "BuGn", title = names(object.list)[i+1], font.size = 8, font.size.title = 8)
draw(ht1 + ht2, ht_gap = unit(1, "cm"))

netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = FALSE, color.text.use = FALSE, angle.x=45)
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = c("NOTCH"), pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = FALSE, color.text.use = FALSE, angle.x=45)
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 2, color.text.use = FALSE, title.name = "Increased signaling in ACI", angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = c("NOTCH"), pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 2, color.text.use = FALSE, title.name = "Increased signaling in ACI", font.size = 20, angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = NULL, pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 1, color.text.use = FALSE, title.name = "Decreased signaling in ACI", angle.x = 45)
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = c("NOTCH"), pairLR.use = NULL, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, max.dataset = 1, color.text.use = FALSE, title.name = "Decreased signaling in ACI", font.size = 20, angle.x = 45)

cellchat <- identifyOverExpressedGenes(object=cellchat, group.by = NULL, group.dataset = "datasets", pos.dataset = "ACI", features.name = "ACI", only.pos = FALSE, thresh.pc = 0, thresh.fc = 0.322, thresh.p = 0.05)
net <- netMappingDEG(object=cellchat, features.name = "ACI", thresh= 0.05)
net.up <- subsetCommunication(object=cellchat, net = net, slot.name = "net", sources.use = NULL, targets.use = NULL, signaling = NULL, thresh = 0.05, datasets = "ACI", ligand.pvalues = 0.05, ligand.logFC = 0.322, receptor.pvalues = 0.05, receptor.logFC = NULL)
net.down <- subsetCommunication(object=cellchat, net = net, slot.name = "net", sources.use = NULL, targets.use = NULL, signaling = NULL, thresh = 0.05, datasets = "WI", ligand.pvalues = 0.05, ligand.logFC = -0.415, receptor.pvalues = 0.05, receptor.logFC = NULL)
gene.up <- extractGeneSubsetFromPair(pairLR=net.up, object = cellchat)
gene.down <- extractGeneSubsetFromPair(pairLR=net.down, object = cellchat)
pairLR.use.up = net.up[, "interaction_name", drop = FALSE]
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = NULL, pairLR.use = pairLR.use.up, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, title.name = paste0("Up-regulated signaling in ACI", names(object.list)[2]), angle.x = 45)
pairLR.use.down = net.down[, "interaction_name", drop = FALSE]
netVisual_bubble(object=cellchat, sources.use = c(1:16), targets.use = c(13), signaling = NULL, pairLR.use = pairLR.use.down, color.heatmap = c("Spectral"), thresh = 0.05, comparison = c(1, 2), remove.isolate = TRUE, title.name = paste0("Down-regulated signaling in ACI", names(object.list)[2]), angle.x = 45)
netVisual_chord_gene(object=object.list[[2]], slot.name = "net", signaling = NULL, pairLR.use = NULL, net = net.up, sources.use = c(1:16), targets.use = c(13), lab.cex = 0.5, title.name = paste0("Up-regulated signaling in ACI", names(object.list)[2]), thresh = 0.05)
netVisual_chord_gene(object=object.list[[1]], slot.name = "net", signaling = NULL, pairLR.use = NULL, net = net.down, sources.use = c(1:16), targets.use = c(13), lab.cex = 0.5, title.name = paste0("Down-regulated signaling in ACI", names(object.list)[1]), thresh = 0.05)

weight.max <- getMaxWeight(object.list=object.list, slot.name = c("netP"), attribute = c("NOTCH")) # control the edge weights across different datasets

par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_aggregate(object=object.list[[i]], signaling= c("NOTCH"), signaling.name = paste(c("NOTCH"), names(object.list)[i]), color.use =leiden_colors, thresh = 0.05, vertex.receiver = NULL, sources.use = NULL, targets.use = NULL, idents.use = NULL, top = 1, remove.isolate = FALSE, vertex.weight = 1, weight.scale = TRUE, edge.weight.max = weight.max[1], edge.width.max = 10, layout = c("circle"), pt.title = 12, title.space = 1, vertex.label.cex = 0.8)
}

par(mfrow = c(1,2), xpd=TRUE)
ht <- list()
for (i in 1:length(object.list)) {
  ht[[i]] <- netVisual_heatmap(object=object.list[[i]], comparison = c(1, 2), measure = c("count"), signaling = c("NOTCH"), slot.name = c("netP"), color.use =leiden_colors, color.heatmap = c("Reds"), title.name = paste(c("NOTCH"), "signaling ",names(object.list)[i]), sources.use = NULL, targets.use = NULL, remove.isolate = FALSE)
}
draw(ht[[1]] + ht[[2]], ht_gap = unit(1, "cm"))

cellchat@meta$datasets = factor(cellchat@meta$datasets, levels = c("WI", "ACI")) # set factor level
plotGeneExpression(object=cellchat, signaling = c("NOTCH"), enriched.only = TRUE, type = c("violin"), split.by = "datasets",colors.ggplot = TRUE)
