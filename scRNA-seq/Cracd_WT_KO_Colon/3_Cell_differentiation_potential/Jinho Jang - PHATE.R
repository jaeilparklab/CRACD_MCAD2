rm(list = ls())

library(data.table)
library(phateR)
library(ggplot2)
library(readr)
library(viridis)
library(Rmagic)

# Prepare input data
count.df = as.data.frame(fread("/Path_to_count_table"))
gene.df = read.table("/Path_to_feature_text", header = T, as.is = T, sep = ",")
meta.df = read.table("/Path_to_meta_data", header = T, as.is = T, sep = ",")
rownames(meta.df) = meta.df[,1]
rownames(rna.df) = gene.df[,2]

output.dir = "/Path_to_output_directory"

input.df = as.data.frame(t(rna.df))

# Keep genes expressed in at least 10 cells
keep_cols <- colSums(input.df > 0) > 10
input.df <- input.df[,keep_cols]

# Keep cells with at least 1000 UMIs
keep_rows <- rowSums(input.df) > 1000
input.df <- input.df[keep_rows,]

input.df <- library.size.normalize(input.df)
input.df <- sqrt(input.df)

input_PCA <- as.data.frame(prcomp(input.df)$x)

# Run PHATE
input_PHATE <- phate(input_PCA)
output.df = as.data.frame(input_PHATE[[1]])
write.table(output.df, paste0(output.dir, "/PHATE.output.txt"), sep = "\t", row.names = T, col.names = T, quote = F)

# Remove outlier
b.index = which(output.df$PHATE2 == max(output.df$PHATE2))
f.ouptut.df = output.df[-b.index,]

# Add cell information for plotting
output.df$CellType = as.character(sapply(as.character(output.df$leiden), function(x){strsplit(x,"_")[[1]][1]}))
type.array = unique(output.df$CellType)

# Specify color order
cell.index = c('Stem', 'Progenitors_1', 'Progenitors_2', 'Progenitors_3', 'Progenitors_4', 'Clonocytes_1', 'Clonocytes_2', 'Clonocytes_3', 'Clonocytes_4', 'Clonocytes_5', 'Clonocytes_6', 'Clonocytes_7', 'Goblet_1', 'Goblet_2', 'Enteroendocrine', 'Tuft')
output.df$CellType = factor(output.df$leiden, levels = cell.index)
color.index = c('#ffa500ff', '#ee82eeff', '#9acd32ff', '#ff69b4ff', '#ff1493ff','#6495edff', '#4169e1ff', '#0000ffff', '#00008bff', '#191970ff','#808080ff', '#000000ff', '#228b22ff', '#3cb371ff', '#800080ff','#ff0000ff')

# Cell type plot (Merge)
pdf(paste0(output.dir, "/PHATE.CellType.Filter.pdf"), width=6, height=5)
ggplot(f.output.df) +
  geom_point(aes(PHATE1, PHATE2, colour = leiden), size = 0.1) +
  scale_colour_manual(values = color.index) +
  labs(colour="Type") +
  xlim(c(min(f.output.df$PHATE1), max(f.output.df$PHATE1))) +
  ylim(c(min(f.output.df$PHATE2), max(f.output.df$PHATE2))) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
dev.off()

# Cell type plot (ACI)
aci.df = subset(f.output.df, Type == "ACI")
pdf(paste0(output.dir, "/PHATE.CellType.ACI.Filter.pdf"), width=6, height=5)
ggplot(aci.df) +
  geom_point(aes(PHATE1, PHATE2, colour = leiden), size = 0.1) +
  scale_colour_manual(values = color.index) +
  labs(colour="Type") +
  xlim(c(min(f.output.df$PHATE1), max(f.output.df$PHATE1))) +
  ylim(c(min(f.output.df$PHATE2), max(f.output.df$PHATE2))) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
dev.off()

# Cell type plot (WI)
wi.df = subset(f.output.df, Type == "WI")
pdf(paste0(output.dir, "/PHATE.CellType.WI.Filter.v2.pdf"), width=6, height=5)
ggplot(wi.df) +
  geom_point(aes(PHATE1, PHATE2, colour = leiden), size = 0.1) +
  scale_colour_manual(values = color.index) +
  labs(colour="Type") +
  xlim(c(min(f.output.df$PHATE1), max(f.output.df$PHATE1))) +
  ylim(c(min(f.output.df$PHATE2), max(f.output.df$PHATE2))) +
  guides(colour = guide_legend(override.aes = list(size = 3)))
dev.off()