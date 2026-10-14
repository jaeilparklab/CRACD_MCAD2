rm(list=ls())

install.packages("arrow")
install.packages('Seurat')


library(Seurat)
library(future)
plan("multisession", workers = 10)
library(ggplot2)
library(sf)
library(reshape2)
library(stringr)


options(future.globals.maxSize = 8000 * 1024^2)


setwd("D:/KP/Colon_Xenium/")






Tumor = readRDS("A_CA_tumor.rds")



######################################################################################
########################## AddModuleScore ###########################################

DefaultAssay(Tumor)

table(Idents(Tumor))
Idents(Tumor) = "genotype_celltype"
head(Tumor@meta.data)

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


table(Idents(Tumor))

Tumor$hippo1 <- APP_hippo$hippo1
Tumor$secret1 <- APP_secret$secret1
Tumor$notch1 <- APP_notch$notch1
Tumor$hif1a1 <- APP_hif1a$hif1a1

A_tumor <- c("A_Epi_1","A_Epi_2","A_Epi_3", "A_Epi_4","A_Epi_5","A_Epi_6")
CA_tumor <- c("CA_Epi_1","CA_Epi_2","CA_Epi_3", "CA_Epi_4","CA_Epi_5","CA_Epi_6", "CA_Epi_7","CA_Epi_8", "CA_Epi_9","CA_Epi_10","CA_Epi_11","CA_Epi_12")

A_tumor <- subset(Tumor, idents = A_tumor)
CA_tumor <- subset(Tumor, idents = CA_tumor)




###################################################################################################################
########################## ROI subsets 1 ##########################################################################

#### CA tumor ROI1 ###

coords <- GetTissueCoordinates(CA_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(CA_tumor)


max(coords$x)
max(coords$y)

ImageDimPlot(CA_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)


selected_cells <- rownames(coords[
  coords$x >= 5000.000 & coords$x <= 6000 &
    coords$y >= 2000 & coords$y <= 3000, ])
length(selected_cells)

CA_tumor@graphs <- list()
CA_tumor@neighbors <- list()
CA_tumor_ROI1 <- subset(CA_tumor, cells = selected_cells)

Images(CA_tumor_ROI1)
ImageDimPlot(CA_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(CA_tumor_ROI1,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(CA_tumor_ROI1, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


#### A tumor ROI1 ###

coords <- GetTissueCoordinates(A_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(A_tumor)


max(coords$x)
max(coords$y)

ImageDimPlot(A_tumor,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)

selected_cells <- rownames(coords[
  coords$x >= 100 & coords$x <= 1100 &
    coords$y >= 600 & coords$y <= 1600, ])
length(selected_cells)

A_tumor@graphs <- list()
A_tumor@neighbors <- list()
A_tumor_ROI1 <- subset(A_tumor, cells = selected_cells)

Images(A_tumor_ROI1)
ImageDimPlot(A_tumor,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(A_tumor_ROI1,fov="fov.2", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(A_tumor_ROI1, fov="fov.2", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))

#hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
#secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
#notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))
#hif1a = list(c("Adm",	"Aldoa",	"Aldoc",	"Angpt2",	"Bnip3",	"Epo","Hk1",	"Hk2",	"Igf2",	"Ldha",	"Nos3",	"Pdk1","Pgk1",	"Tfrc",	"Tgfa",	"Vegfa"))

#APP_hippo <- AddModuleScore(object = Tumor, features = hippo, name = "hippo", assay = "SCT", nbin=4)
#head(APP_hippo)
#APP_secret <- AddModuleScore(object = Tumor, features = secret, name = "secret", assay = "SCT", nbin=4)
#head(APP_secret)
#APP_notch <- AddModuleScore(object = Tumor, features = notch, name = "notch", assay = "SCT", nbin=4)
#head(APP_notch)
#APP_hif1a <- AddModuleScore(object = Tumor, features = hif1a, name = "hif1a", assay = "SCT", nbin=4)
#head(APP_hif1a)

#CA_tumor_ROI1$hippo1 <- APP_hippo$hippo1
#CA_tumor_ROI1$secret1 <- APP_secret$secret1
#CA_tumor_ROI1$notch1 <- APP_notch$notch1
#CA_tumor_ROI1$hif1a1 <- APP_hif1a$hif1a1

#A_tumor_ROI1$hippo1 <- APP_hippo$hippo1
#A_tumor_ROI1$secret1 <- APP_secret$secret1
#A_tumor_ROI1$notch1 <- APP_notch$notch1
#A_tumor_ROI1$hif1a1 <- APP_hif1a$hif1a1


DotPlot(A_tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(CA_tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(A_tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(CA_tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")


DotPlot(A_tumor_ROI1, features = c("Muc2","Spdef"), cols = "RdYlBu",scale.min = 0,scale.max = 50, group.by = "genotype", split.by = "celltype")
DotPlot(CA_tumor_ROI1, features = c("Muc2","Spdef"), cols = "RdYlBu",scale.min = 0,scale.max = 50, group.by = "genotype", split.by = "celltype")


Tumor_ROI1 <- merge(A_tumor_ROI1, y=CA_tumor_ROI1, add.cell.ids=c("A_ROI1","CA_ROI1"))
table(Idents(Tumor_ROI1))

DotPlot(Tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")
DotPlot(Tumor_ROI1, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


DotPlot(Tumor_ROI1, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI1, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


ImageDimPlot(Tumor_ROI1,fov="fov", size=2.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(Tumor_ROI1,fov="fov.2", size=2.0, cols = "polychrome", axes=TRUE)

#CA
ImageFeaturePlot(Tumor_ROI1, fov="fov", features = c("Yap1"), max.cutoff = 0.5, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
    plot.background = element_rect(fill = "black", color = NA),
    panel.grid = element_blank()
  )
#A
ImageFeaturePlot(Tumor_ROI1, fov="fov.2", features = c("Yap1"), max.cutoff = 0.5, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )
#CA
ImageFeaturePlot(Tumor_ROI1, fov="fov", features = c("Hes1"), max.cutoff = 0.6, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )
#A
ImageFeaturePlot(Tumor_ROI1, fov="fov.2", features = c("Hes1"), max.cutoff = 0.6, size = 2, cols = c("mediumslateblue", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
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


#CA, A
ImageFeaturePlot(Tumor_ROI1, fov=c("fov","fov.2"), features = c("notch1"), max.cutoff =0.03, min.cutoff=0.0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )

#CA, A
ImageFeaturePlot(Tumor_ROI1, fov=c("fov","fov.2"), features = c("hippo1"), max.cutoff =0.07, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )



###################################################################################################################
########################## ROI subsets 2 ##########################################################################

#### KO tumor ROI2 ###

coords <- GetTissueCoordinates(CA_tumor)
head(coords)
tail (coords)
rownames(coords) <- colnames(CA_tumor)


max(coords$x)
max(coords$y)

ImageDimPlot(CA_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)


selected_cells <- rownames(coords[
  coords$x >= 6500 & coords$x <= 7500 &
    coords$y >= 1600 & coords$y <= 2600, ])
length(selected_cells)

CA_tumor@graphs <- list()
CA_tumor@neighbors <- list()
CA_tumor_ROI2 <- subset(CA_tumor, cells = selected_cells)

Images(CA_tumor_ROI2)
ImageDimPlot(CA_tumor,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageDimPlot(CA_tumor_ROI2,fov="fov", size=1.0, cols = "polychrome", axes=TRUE)
ImageFeaturePlot(CA_tumor_ROI2, fov="fov", features = c("Muc2", "Spdef"), max.cutoff = c(0.5,0.5), size = 1.5, cols = c("white", "red"))


#### No more ROI available for A tumor ###
######### ROI1 is used again for A tumor####


Tumor_ROI2 <- merge(A_tumor_ROI1, y=CA_tumor_ROI2, add.cell.ids=c("A_ROI1","CA_ROI2"))
table(Idents(Tumor_ROI2))

#hippo = list(c("Ajuba",	"Amot",	"Amotl1",	"Amotl2",	"Arrdc3",	"Cit",	"Dchs1",	"Dlg5",	"Fat4",	"Frmd6",	"Lats1",	"Lats2",	"Limd1",	"Map2k3",	"Mapk14",	"Mark3",	"Mob1a",	"Mob1b",	"Mob3b",	"Nek8",	"Nf2",	"Nuak2",	"Sav1",	"Schip1",	"Shank2",	"Sox11",	"Src",	"Stk3",	"Stk4",	"Tead1",	"Tead2",	"Tead3",	"Tead4",	"Tial1",	"Vgll4",	"Wtip",	"Wwc1",	"Wwc2",	"Wwtr1",	"Yap1"))
#secret = list(c("Rab11b",	"Rab27a",	"Rab33b",	"Tmem167",	"Tmem167b",	"Trappc11"))
#notch = list(c("Aak1",	"Acvrl1",	"Adam10",	"Adam17",	"Agxt",	"Akt1s1",	"Angpt4",	"Anxa4",	"Aph1a",	"Aph1b",	"Aph1c",	"App",	"Arrb1",	"Arrdc1",	"Ascl1",	"Atoh1",	"Bcl6",	"Bend6",	"Bmp2",	"Bmp2k",	"Bmp7",	"Cbfa2t2",	"Cc2d1a",	"Cc2d1b",	"Ccdc85c",	"Ccn3",	"Ccnc",	"Cd46",	"Cdh6",	"Cdk3",	"Cdk6",	"Cdkn1b",	"Cebpa",	"Cfap58",	"Cfd",	"Chac1",	"Cntn1",	"Cntn6",	"Dicer1",	"Dlk1",	"Dlk2",	"Dll1",	"Dll3",	"Dll4",	"Dlx1",	"Dlx2",	"Dner",	"Dtx1",	"Dtx2",	"Dtx3",	"Dtx3l",	"Dtx4",	"Egfl7",	"Enho",	"Epn1",	"Epn2",	"Etv2",	"Eya1",	"Fat4",	"Fbxw7",	"Fgf10",	"Foxa1",	"Foxa2",	"Foxc1",	"Foxc2",	"Galnt11",	"Gas2",	"Gata2",	"Gata5",	"Gdf2",	"Gdpd5",	"Got1",	"Grip2",	"Gsx2",	"Herc2",	"Hes1",	"Hes5",	"Hes7",	"Hey1",	"Hey2",	"Heyl",	"Hhex",	"Hif1an",	"Hnf1b",	"Hoxd3",	"Hp",	"Ift172",	"Ift74",	"Ift88",	"Il17a",	"Il2ra",	"Il6st",	"Itgb1bp1",	"Jag1",	"Jag2",	"Kcna5",	"Kctd10",	"Kit",	"Krt19",	"Lfng",	"Llgl1",	"Llgl2",	"Lrrk2",	"Maml1",	"Maml2",	"Maml3",	"Mecp2",	"Mesp1",	"Mesp2",	"Mfng",	"Mib1",	"Mib2",	"Mmp14",	"Ncstn",	"Nepro",	"Neurl1a",	"Neurl1b",	"Neurod4",	"Nfkbia",	"Nkap",	"Nle1",	"Nod2",	"Notch1",	"Notch2",	"Notch3",	"Notch4",	"Nr0b2",	"Nr1h4",	"Nrarp",	"Nrip2",	"Numb",	"Onecut1",	"Ovol2",	"Pdcd10",	"Perp",	"Pgam2",	"Pln",	"Pofut1",	"Poglut1",	"Postn",	"Prag1",	"Psen1",	"Psen2",	"Psenen",	"Ptp4a3",	"Rbm15",	"Rbpj",	"Rcan2",	"Reck",	"Rfng",	"Rian",	"Ripply1",	"Ripply2",	"Rita1",	"Robo1",	"Robo2",	"Rps19",	"S1pr3",	"S2bpcox16",	"Sel1l",	"Six1",	"Slc35c1",	"Slc35c2",	"Snai1",	"Snai2",	"Sorbs2",	"Sox2",	"Sox9",	"Spen",	"Src",	"Stat3",	"Susd5",	"Synj2bp",	"Tbx2",	"Tcim",	"Tgfb1",	"Tgfb2",	"Tgfbr2",	"Timp4",	"Tm9sf5",	"Tmem100",	"Traf7",	"Trp63",	"Tspan14",	"Tspan15",	"Tspan5",	"Tspear",	"Wdr12",	"Wnt1",	"Xirp1",	"Yap1",	"Zbtb7a",	"Zfp423",	"Zmiz1"))
#hif1a = list(c("Adm",	"Aldoa",	"Aldoc",	"Angpt2",	"Bnip3",	"Epo","Hk1",	"Hk2",	"Igf2",	"Ldha",	"Nos3",	"Pdk1","Pgk1",	"Tfrc",	"Tgfa",	"Vegfa"))


#APP_hippo <- AddModuleScore(object = Tumor_ROI2, features = hippo, name = "hippo", assay = "SCT", nbin=4)
#APP_secret <- AddModuleScore(object = Tumor_ROI2, features = secret, name = "secret", assay = "SCT", nbin=4)
#APP_notch <- AddModuleScore(object = Tumor_ROI2, features = notch, name = "notch", assay = "SCT", nbin=4)
#APP_hif1a <- AddModuleScore(object = Tumor_ROI2, features = hif1a, name = "hif1a", assay = "SCT", nbin=4)

#Tumor_ROI2$hippo1 <- APP_hippo$hippo1
#Tumor_ROI2$secret1 <- APP_secret$secret1
#Tumor_ROI2$notch1 <- APP_notch$notch1
#Tumor_ROI2$hif1a1 <- APP_hif1a$hif1a1
tail(Tumor_ROI2)

DotPlot(Tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")
DotPlot(Tumor_ROI2, features = c("hippo1","secret1","notch1"), cols = "RdYlBu",scale.min = 0,scale.max = 50, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "celltype")
DotPlot(Tumor_ROI2, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


DotPlot(Tumor_ROI2, features = c("Muc2","Spdef","Yap1","Ccn1","Ccn2","Axl","Hes1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")






head(Tumor_ROI2)


# CA A
ImageFeaturePlot(Tumor_ROI2, fov=c("fov","fov.2"), features = c("notch1"), max.cutoff =0.03, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )

# CA A
ImageFeaturePlot(Tumor_ROI2, fov=c("fov","fov.2"), features = c("hippo1"), max.cutoff =0.07, min.cutoff=0,  size = 1.5, cols = c("orchid1", "yellow"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )


# CA A
ImageFeaturePlot(Tumor_ROI1, fov=c("fov","fov.2"), features = c("hif1a1"), max.cutoff =0.07, min.cutoff=0,  size = 1.5, cols = c("white", "red"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )
ImageFeaturePlot(Tumor_ROI2, fov=c("fov","fov.2"), features = c("hif1a1"), max.cutoff =0.07, min.cutoff=0,  size = 1.5, cols = c("white", "red"), border.color= "black", border.size = 0.05, alpha = 0.85, axes=TRUE)+
  theme(panel.background = element_rect(fill = "black", color = NA),  # or "black", "grey90", etc.
        plot.background = element_rect(fill = "black", color = NA),
        panel.grid = element_blank()
  )



####### Merge A_ROI1, CA_ROI1, CA_ROI2 for plot #########

Tumor_ROI_all <- merge(A_tumor_ROI1, y=c(CA_tumor_ROI1, CA_tumor_ROI2), add.cell.ids=c("A_ROI1","CA_ROI1","CA_ROI2"))
head(Tumor_ROI_all)
Tumor_ROI_all$ROI <- str_extract(colnames(Tumor_ROI_all), "^[^_]+_[^_]+")
table(Tumor_ROI_all$ROI)

DotPlot(Tumor_ROI_all, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "ROI")
DotPlot(Tumor_ROI_all, features = c("hif1a1"), cols = "RdYlBu",scale.min = 0,scale.max = 10, col.max = 2.5, col.min = -1.5, group.by = "genotype", split.by = "genotype")


saveRDS(Tumor_ROI_all, "A_CA_ROIs_merged.rds")
