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





saveRDS(Tumor, "WT_KO_Tumor_merged.rds")
Tumor = readRDS("WT_KO_Tumor_merged.rds")
head(Tumor)


######################################################################################
########################## AddModuleScore ###########################################

DefaultAssay(Tumor)
head(rownames(Tumor))

options(warn=1)

#wnt = list(c("Axin2","Tbx3"))
#tgfb = list(c("Serpine1","Smad7"))
#hh = list(c("Gli","Ptch1"))
hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))
hif1a = list(c("Adm",	"Aldoa",	"Aldoc",	"Angpt2",	"Bnip3",	"Epo","Hk1",	"Hk2",	"Igf2",	"Ldha",	"Nos3",	"Pdk1","Pgk1",	"Tfrc",	"Tgfa",	"Vegfa"))



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

APP_hif1a <- AddModuleScore(object = Tumor, features = hif1a, name = "hif1a", assay = "SCT", nbin=4)
head(APP_hif1a)
DotPlot(APP_hif1a, features = "hif1a1", cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(APP_hif1a, features = "hif1a1", cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


Tumor$hippo1 <- APP_hippo$hippo1
Tumor$secret1 <- APP_secret$secret1
Tumor$notch1 <- APP_notch$notch1
Tumor$hif1a1 <- APP_hif1a$hif1a1


WT_tumor <- c("WT_Epi_1","WT_Epi_2","WT_Epi_3", "WT_Epi_4","WT_Epi_5","WT_Epi_6", "WT_Epi_7","WT_Epi_8","WT_Epi_9", "WT_Epi_10","WT_Epi_11","WT_Epi_12","WT_Epi_13", "WT_Epi_14","WT_Epi_15")
KO_tumor <- c("KO_Epi_1","KO_Epi_2","KO_Epi_3", "KO_Epi_4","KO_Epi_5","KO_Epi_6", "KO_Epi_7","KO_Epi_8", "KO_Epi_9","KO_Epi_10","KO_Epi_11","KO_Epi_12", "KO_Epi_13")

WT_tumor <- subset(Tumor, idents = WT_tumor)
KO_tumor <- subset(Tumor, idents = KO_tumor)


WT_tumor@images
KO_tumor@images

ImageFeaturePlot(KO_tumor, fov="fov.2", features = c("hippo1", "secret1", "notch1"), max.cutoff = c(0.15,0.8,0.05), min.cutoff = -0.5, size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(WT_tumor, fov="fov", features = c("hippo1", "secret1", "notch1"), max.cutoff = c(0.15,0.8,0.05), min.cutoff = -0.5, size = 1.5, cols = c("white", "red"))

p <- ImageFeaturePlot(
  object = KO_tumor,
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

DotPlot(KO_tumor, features = c("Muc2","Spdef","Atoh1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, group.by = "genotype", split.by = "celltype")
DotPlot(WT_tumor, features = c("Muc2","Spdef","Atoh1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, group.by = "genotype", split.by = "celltype")


ImageFeaturePlot(KO_tumor, fov="fov.2", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))
ImageFeaturePlot(WT_tumor, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))



ImageDimPlot(KO_tumor, fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(WT_tumor, fov="fov", size=1.0, cols = "polychrome", axes=TRUE)

Idents(KO_tumor) <- "genotype_celltype"
Idents(WT_tumor) <- "genotype_celltype"

A <- ImageDimPlot(
  object = KO_tumor, size=1.5,
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
  object = WT_tumor, size=1.5,
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
KO_tumor.mat = KO_tumor@assays$SCT@counts
A.cells = colnames(KO_tumor@assays$SCT@counts)
A.genes = rownames(KO_tumor@assays$SCT@counts)
rownames(KO_tumor.mat) = A.genes
colnames(KO_tumor.mat) = paste0("KO_", A.cells)

WT_tumor.mat = WT_tumor@assays$SCT@counts
CA.cells = colnames(WT_tumor@assays$SCT@counts)
CA.genes = rownames(WT_tumor@assays$SCT@counts)
rownames(WT_tumor.mat) = CA.genes
colnames(WT_tumor.mat) = paste0("WT_", CA.cells)

output.df = data.frame()
for(i in 1:length(CA.genes)) {
  cat(i, "\t", length(CA.genes), "\n")
  tmp.gene = CA.genes[i]
  A.array = as.numeric(KO_tumor.mat[tmp.gene,])
  CA.array = as.numeric(WT_tumor.mat[tmp.gene,])
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
KO_scores <- KO_tumor@meta.data[, c("hippo1", "secret1", "notch1", "hif1a1")]
WT_scores <- WT_tumor@meta.data[, c("hippo1", "secret1", "notch1","hif1a1")]

# Step 2: Add row labels for tracking
rownames(KO_scores) <- paste0("KO_", rownames(KO_scores))
rownames(WT_scores) <- paste0("WT_", rownames(WT_scores))

# Step 3: Merge the data (if needed later for plotting)
all_scores <- rbind(KO_scores, WT_scores)
group <- c(rep("A", nrow(KO_scores)), rep("CA", nrow(WT_scores)))

# Step 4: Run t-tests and collect means
features <- colnames(KO_scores)
output.df <- data.frame()

for (i in seq_along(features)) {
  feature <- features[i]
  KO_vals <- KO_scores[[feature]]
  WT_vals <- WT_scores[[feature]]
  
  test.result <- t.test(KO_vals, WT_vals)
  pval <- test.result$p.value
  KO_mean <- mean(KO_vals)
  WT_mean <- mean(WT_vals)
  
  output.df[i, "Feature"] <- feature
  output.df[i, "pval"] <- pval
  output.df[i, "KO_mean"] <- KO_mean
  output.df[i, "WT_mean"] <- WT_mean
  output.df[i, "Higher_in"] <- ifelse(KO_mean > WT_mean, "A", "CA")
}

# View output
print(output.df)


########################## Vioin plots of scores by genotype ##########################################################################


VlnPlot(Tumor, features = c("hippo1",'secret1',"notch1"), group.by = "genotype", split.by = "genotype", ncol = 2, pt.size = 0,
        cols = c("coral", "deepskyblue2"))




###################################################################################################################
########################## ROI subsets 1 ##########################################################################

#### KO tumor ROI1 ###

coords <- GetTissueCoordinates(KO_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(KO_tumor)


max(coords$x)
max(coords$y)


selected_cells <- rownames(coords[
  coords$x >= 5000.000 & coords$x <= 6000 &
    coords$y >= 5000 & coords$y <= 6000, ])
length(selected_cells)

KO_tumor@graphs <- list()
KO_tumor@neighbors <- list()
KO_tumor_ROI1 <- subset(KO_tumor, cells = selected_cells)

Images(KO_tumor_ROI1)
ImageDimPlot(KO_tumor,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(KO_tumor_ROI1,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(KO_tumor_ROI1, fov="fov.2", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


#### WT tumor ROI1 ###

coords <- GetTissueCoordinates(WT_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(WT_tumor)


max(coords$x)
max(coords$y)


selected_cells <- rownames(coords[
  coords$x >= 2500 & coords$x <= 3500 &
    coords$y >= 2000 & coords$y <= 3000.000, ])
length(selected_cells)

WT_tumor@graphs <- list()
WT_tumor@neighbors <- list()
WT_tumor_ROI1 <- subset(WT_tumor, cells = selected_cells)

Images(WT_tumor_ROI1)
ImageDimPlot(WT_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(WT_tumor_ROI1,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(WT_tumor_ROI1, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


#KO_tumor_ROI1$hippo1 <- APP_hippo$hippo1
#KO_tumor_ROI1$secret1 <- APP_secret$secret1
#KO_tumor_ROI1$notch1 <- APP_notch$notch1
#WT_tumor_ROI1$hippo1 <- APP_hippo$hippo1
#WT_tumor_ROI1$secret1 <- APP_secret$secret1
#WT_tumor_ROI1$notch1 <- APP_notch$notch1

DotPlot(WT_tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(KO_tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")

Tumor_ROI1 <- merge(WT_tumor_ROI1, y=KO_tumor_ROI1, add.cell.ids=c("WT_ROI1","KO_ROI1"))
table(Idents(Tumor_ROI1))

DotPlot(Tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")
DotPlot(Tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


DotPlot(Tumor_ROI1, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


ImageDimPlot(Tumor_ROI1,fov="fov", size=2.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(Tumor_ROI1,fov="fov.2", size=2.0, cols = "polychrome", axes=TRUE)

ImageFeaturePlot(Tumor_ROI1, fov="fov", features = c("Yap1"), max.cutoff = 0.8, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid = element_blank()
  )
ImageFeaturePlot(Tumor_ROI1, fov="fov.2", features = c("Yap1"), max.cutoff = 0.8, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )
ImageFeaturePlot(Tumor_ROI1, fov="fov", features = c("Hes1"), max.cutoff = 0.8, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )
ImageFeaturePlot(Tumor_ROI1, fov="fov.2", features = c("Hes1"), max.cutoff = 0.8, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )

#hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
#secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
#notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))


#APP_hippo <- AddModuleScore(object = Tumor_ROI1, features = hippo, name = "hippo", assay = "SCT", nbin=4)
#APP_secret <- AddModuleScore(object = Tumor_ROI1, features = secret, name = "secret", assay = "SCT", nbin=4)
#APP_notch <- AddModuleScore(object = Tumor_ROI1, features = notch, name = "notch", assay = "SCT", nbin=4)

#Tumor_ROI1$hippo1 <- APP_hippo$hippo1
#Tumor_ROI1$secret1 <- APP_secret$secret1
#Tumor_ROI1$notch1 <- APP_notch$notch1

head(Tumor_ROI1)



ImageFeaturePlot(Tumor_ROI1, fov=c("fov","fov.2"), features = c("notch1"), max.cutoff =0.02, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )


ImageFeaturePlot(Tumor_ROI1, fov=c("fov","fov.2"), features = c("hippo1"), max.cutoff =0.02, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )



###################################################################################################################
########################## ROI subsets 2 ##########################################################################

#### KO tumor ROI2 ###

coords <- GetTissueCoordinates(KO_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(KO_tumor)


max(coords$x)
max(coords$y)

ImageDimPlot(KO_tumor,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)


selected_cells <- rownames(coords[
  coords$x >= 2000 & coords$x <= 3000 &
    coords$y >= 4000 & coords$y <= 5000, ])
length(selected_cells)

KO_tumor@graphs <- list()
KO_tumor@neighbors <- list()
KO_tumor_ROI2 <- subset(KO_tumor, cells = selected_cells)

Images(KO_tumor_ROI2)
ImageDimPlot(KO_tumor,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(KO_tumor_ROI2,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(KO_tumor_ROI2, fov="fov.2", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


#### WT tumor ROI2 ###

coords <- GetTissueCoordinates(WT_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(WT_tumor)


max(coords$x)
max(coords$y)

ImageDimPlot(WT_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)

selected_cells <- rownames(coords[
  coords$x >= 2500 & coords$x <= 3500 &
    coords$y >= 1000 & coords$y <= 2000, ])
length(selected_cells)

WT_tumor@graphs <- list()
WT_tumor@neighbors <- list()
WT_tumor_ROI2 <- subset(WT_tumor, cells = selected_cells)

Images(WT_tumor_ROI2)
ImageDimPlot(WT_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(WT_tumor_ROI2,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(WT_tumor_ROI2, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))



#KO_tumor_ROI2$hippo1 <- APP_hippo$hippo1
#KO_tumor_ROI2$secret1 <- APP_secret$secret1
#KO_tumor_ROI2$notch1 <- APP_notch$notch1
#WT_tumor_ROI2$hippo1 <- APP_hippo$hippo1
#WT_tumor_ROI2$secret1 <- APP_secret$secret1
#WT_tumor_ROI2$notch1 <- APP_notch$notch1

DotPlot(WT_tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(KO_tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")


Tumor_ROI2 <- merge(WT_tumor_ROI2, y=KO_tumor_ROI2, add.cell.ids=c("WT_ROI2","KO_ROI2"))
table(Idents(Tumor_ROI2))

DotPlot(Tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")
DotPlot(Tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")



DotPlot(Tumor_ROI2, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")



hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))


APP_hippo <- AddModuleScore(object = Tumor_ROI2, features = hippo, name = "hippo", assay = "SCT", nbin=4)
APP_secret <- AddModuleScore(object = Tumor_ROI2, features = secret, name = "secret", assay = "SCT", nbin=4)
APP_notch <- AddModuleScore(object = Tumor_ROI2, features = notch, name = "notch", assay = "SCT", nbin=4)

Tumor_ROI2$hippo1 <- APP_hippo$hippo1
Tumor_ROI2$secret1 <- APP_secret$secret1
Tumor_ROI2$notch1 <- APP_notch$notch1

head(Tumor_ROI2)



ImageFeaturePlot(Tumor_ROI2, fov=c("fov","fov.2"), features = c("notch1"), max.cutoff =0.02, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )


ImageFeaturePlot(Tumor_ROI2, fov=c("fov","fov.2"), features = c("hippo1"), max.cutoff =0.02, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )





####### Merge A_ROI1, CA_ROI1, CA_ROI2 for plot #########

Tumor_ROI_all <- merge(WT_tumor_ROI1, y=c(WT_tumor_ROI2, KO_tumor_ROI1, KO_tumor_ROI2), add.cell.ids=c("WT_ROI1","WT_ROI2","KO_ROI1","KO_ROI2"))
head(Tumor_ROI_all)
Tumor_ROI_all$ROI <- str_extract(colnames(Tumor_ROI_all), "^[^_]+_[^_]+")
table(Tumor_ROI_all$ROI)

DotPlot(Tumor_ROI_all, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "ROI")

saveRDS(Tumor_ROI_all, "AOMDSS_WT_KO_ROIs_merged.rds")
