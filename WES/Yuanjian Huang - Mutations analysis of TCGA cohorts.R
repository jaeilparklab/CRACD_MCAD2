library(maftools)
library(mclust)
library(BSgenome.Hsapiens.UCSC.hg38)
library(NMF)
library(pheatmap)
library(data.table)
library(ASCAT)
library(RColorBrewer)
library(dplyr)

# Download TCGA COAD and READ MAF data from MC3 repository
coad = tcgaLoad(study = "COAD", source = c("MC3"), repo = c("gitee"))
read = tcgaLoad(study = "READ", source = c("MC3"), repo = c("gitee"))

# Check for duplicated samples between COAD and READ datasets
coad_samples <- getSampleSummary(coad)$Tumor_Sample_Barcode
read_samples <- getSampleSummary(read)$Tumor_Sample_Barcode
duplicated_samples <- intersect(coad_samples, read_samples)
length(duplicated_samples) > 0

# Merge COAD and READ MAF objects
coad_read_maf <- merge_mafs(maf = list(coad, read))
coad_read_maf
colnames(coad_read_maf@data)
write.mafSummary(maf = coad_read_maf, basename = 'coad_read_maf')

# Check clinical data columns and histological types
colnames(coad_read_maf@clinical.data)
unique(coad_read_maf@clinical.data$histological_type)

# Recode histological types for consistency
clin_data <- coad_read_maf@clinical.data
clin_data <- clin_data %>%
  mutate(
    histological_type = case_when(
      histological_type %in% c("Colon_Adenocarcinoma", "Rectal_Adenocarcinoma") ~ "Colorectal_Adenocarcinoma",
      histological_type %in% c("Colon_Mucinous_Adenocarcinoma", "Rectal_Mucinous_Adenocarcinoma") ~ "Colorectal_Mucinous_Adenocarcinoma",
      TRUE ~ histological_type
    )
  )
coad_read_maf@clinical.data <- clin_data
unique(coad_read_maf@clinical.data$histological_type)

# Filter MAF
filtered_maf = subsetMaf(maf = coad_read_maf, clinQuery = "histological_type %in% c('Colorectal_Adenocarcinoma', 'Colorectal_Mucinous_Adenocarcinoma')")
unique(filtered_maf@clinical.data$histological_type)
adeno_count <- sum(filtered_maf@clinical.data$histological_type == "Colorectal_Adenocarcinoma", na.rm = TRUE)
adeno_count
mucinous_count <- sum(filtered_maf@clinical.data$histological_type == "Colorectal_Mucinous_Adenocarcinoma", na.rm = TRUE)
mucinous_count
write.mafSummary(maf = filtered_maf, basename = 'filtered_maf')
plotmafSummary(maf = filtered_maf, rmOutlier = TRUE, dashboard = TRUE, titvRaw = FALSE, log_scale = FALSE, addStat = 'median', showBarcodes = FALSE)

# Define annotation colors for histological types
fabcolors <- c("blue", "red")
fabcolors
names(fabcolors) = c('Colorectal_Adenocarcinoma', 'Colorectal_Mucinous_Adenocarcinoma')
fabcolors
anno_cols = list(histological_type = fabcolors)
anno_cols
oncoplot(maf = filtered_maf, genes=c('CRACD', 'KIAA1211'), clinicalFeatures = 'histological_type', annotationColor = anno_cols, sortByAnnotation = TRUE)

# 1. Extract mutation data for KIAA1211 gene
kiaa1211_var <- filtered_maf@data %>%
  filter(Hugo_Symbol == "KIAA1211") %>%
  select(Cohort, Tumor_Sample_Barcode, Hugo_Symbol, Variant_Classification)

# 2. Extract sample ID and histological type from clinical data
clin_info <- filtered_maf@clinical.data %>%
  select(Tumor_Sample_Barcode, histological_type) %>%
  distinct()

# 3. Merge mutation data with clinical information
kiaa1211_full <- left_join(
  x = kiaa1211_var,
  y = clin_info,
  by = "Tumor_Sample_Barcode"
) %>%
  arrange(Tumor_Sample_Barcode)

# 4. Save full data
write.table(
  x = kiaa1211_full,
  file = "C:/Users/hyj/Documents/CRACD_mutations_TCGA.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# 5. Deduplicate by sample ID
kiaa1211_unique <- kiaa1211_full %>%
  group_by(Cohort, Tumor_Sample_Barcode, histological_type) %>%
  summarise(
    Hugo_Symbol = unique(Hugo_Symbol) %>% paste(collapse = "; "),
    Variant_Classification = unique(Variant_Classification) %>% paste(collapse = "; "),
    .groups = "drop"
  ) %>%
  arrange(Tumor_Sample_Barcode)

# 6. Save deduplicated data
write.table(
  x = kiaa1211_unique,
  file = "C:/Users/hyj/Documents/CRACD_mutations_TCGA_unique.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)

# 7. Print
cat("=== Data extraction and saving completed ===\n")
cat("1. Full data (with multiple mutations per sample): CRACD_mutations.txt\n")
cat("   Total rows (including duplicates): ", nrow(kiaa1211_full), "\n")
cat("2. Deduplicated data (one row per sample): CRACD_mutations_unique.txt\n")
cat("   Unique samples: ", nrow(kiaa1211_unique), "\n\n")



# Load in-house MAF files
C1003 = read.maf(maf = "J:/HYJ/Seq_data/2022_7_28_WES/Analysis/annotation/mutect2/C1033/C1033.mutect2.filtered_VEP.ann.nochr.maf", rmFlags = FALSE, removeDuplicatedVariants = TRUE, useAll = TRUE)
C1144 = read.maf(maf = "J:/HYJ/Seq_data/2022_7_28_WES/Analysis/annotation/mutect2/C1144/C1144.mutect2.filtered_VEP.ann.nochr.maf", rmFlags = FALSE, removeDuplicatedVariants = TRUE, useAll = TRUE)
C1239 = read.maf(maf = "J:/HYJ/Seq_data/2022_7_28_WES/Analysis/annotation/mutect2/C1239/C1239.mutect2.filtered_VEP.ann.nochr.maf", rmFlags = FALSE, removeDuplicatedVariants = TRUE, useAll = TRUE)
C1314 = read.maf(maf = "J:/HYJ/Seq_data/2022_7_28_WES/Analysis/annotation/mutect2/C1314/C1314.mutect2.filtered_VEP.ann.nochr.maf", rmFlags = FALSE, removeDuplicatedVariants = TRUE, useAll = TRUE)
F3053 = read.maf(maf = "J:/HYJ/Seq_data/2022_7_28_WES/Analysis/annotation/mutect2/F3053/F3053.mutect2.filtered_VEP.ann.nochr.maf", rmFlags = FALSE, removeDuplicatedVariants = TRUE, useAll = TRUE)
in_house <- merge_mafs(maf = list(C1003, C1144, C1239, C1314, F3053))
colnames(in_house@clinical.data)
in_house@clinical.data$histological_type <- "Colorectal_Mucinous_Adenocarcinoma"
colnames(in_house@clinical.data)
unique(in_house@clinical.data$histological_type)
write.mafSummary(maf = in_house, basename = 'in_house')

plotmafSummary(maf = in_house, rmOutlier = TRUE, dashboard = TRUE, titvRaw = FALSE, log_scale = FALSE, addStat = 'median', showBarcodes = FALSE)
oncoplot(maf = in_house, genes=c('CRACD', 'KIAA1211'))

# 1. Extract CRACD mutations from in-house MAF data
cracd_mutations <- in_house@data %>%
  filter(Hugo_Symbol == "CRACD") %>%
  select(
    Tumor_Sample_Barcode,
    Hugo_Symbol,
    Variant_Classification
  )

# 2. Add histological type annotation
cracd_full <- cracd_mutations %>%
  mutate(
    histological_type = "Colorectal_Mucinous_Adenocarcinoma"
  ) %>%
  arrange(Tumor_Sample_Barcode)

# 3. Save
write.table(
  x = cracd_full,
  file = "C:/Users/hyj/Documents/CRACD_mutations_in-house.txt",
  sep = "\t",
  row.names = FALSE,
  quote = FALSE,
  na = ""
)



# Merge TCGA and in-house datasets
combined <- merge_mafs(maf = list(coad, read, C1003, C1144, C1239, C1314, F3053))
plotmafSummary(maf = combined, rmOutlier = TRUE, dashboard = TRUE, titvRaw = FALSE, log_scale = FALSE, addStat = 'median', showBarcodes = FALSE)
oncoplot(maf = combined, genes=c('CRACD', 'KIAA1211'))



# Survival analysis for KIAA1211 in TCGA COAD cohort
mafSurvival(maf = coad, genes = 'KIAA1211', time = 'days_to_last_followup', Status = 'surv_event', isTCGA = TRUE)
