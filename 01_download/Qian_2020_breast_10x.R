#### Qian - breast ####

suppressPackageStartupMessages({
  library(Seurat)
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(EnsDb.Hsapiens.v86)
})

setwd("~/sc_breast/data/masters_project/01_download/qian/")

#### Download data ####

#### Read data ####
counts <- Read10X(data.dir = "./raw_data")

unique(substr(colnames(counts), 1, 10))

#### Prepare metadata ####
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0355-0/MediaObjects/41422_2020_355_MOESM13_ESM.pdf
# https://static-content.springer.com/esm/art%3A10.1038%2Fs41422-020-0355-0/MediaObjects/41422_2020_355_MOESM14_ESM.pdf
metadata <- as.tibble(fread("BC_metadata.csv"))

meta.paper = tibble(
  Qian_Patient_number = c(paste0(1:14)),
  Qian_PatientNumber = c(41, 42, 43, 44, 45, 46, 47, 48, 49, 50, 51, 52, 53, 54),
  ## Conditions to the integrated data:
  grade = c(rep(3, 13), 2),
  study = 'Qian et al',
  tech = '10x',
  histological_subtype = c(rep('Ductal',4), 'Apocrine', 'Ductal','Metaplastic', rep('Ductal',7)),
  age_range = c('40-45', '50-55', '50-55', '86-90', '60-65', '70-75', '80-85', '60-65', '40-45', '50-55', '46-50', '76-80', '60-65', '46-50'),
  molecular_subtype = c("HER2", rep("TNBC", 4), "HER2", "TNBC", "HER2", "ER/PR/HER2", rep("TNBC", 3), rep("ER+", 2))
)

metadata = inner_join(metadata %>% rename_with( ~ gsub("^", "Qian_", .x)), meta.paper)


metadata = metadata %>%
  ## Conditions to the integrated data:
  mutate(
    patient = paste0("Qian_", metadata$Qian_Patient_number),
    tumor_type = "breast",
    tech = "10x",
    tissue_origin = 'tumor',
    sample_type = "tumor",
    treatment_group = "treatment_naive",
    tumor_site = "primary",
    menopause = "unknown",
    gender = "F",
    BRCA_status = ifelse(patient == "Qian_3", "Positive", "Negative")
  ) %>%
  column_to_rownames("Qian_Cell")
metadata$condition <- ifelse(metadata$sample_type == "tumor", paste0("tumor_", metadata$molecular_subtype),meta$condition)



convert_and_make_seurat <- function(x) {
  ensembl <-
    mapIds(EnsDb.Hsapiens.v86, rownames(counts), 'GENEID', 'SYMBOL') %>% as.data.frame() 
  counts <- counts %>% as.data.frame()
  counts$gene <- ensembl$.
  counts <- counts %>%
    dplyr::filter(!is.na(gene))
  gene_id <-
    mapIds(EnsDb.Hsapiens.v86, counts$gene, 'SYMBOL', 'GENEID') %>% as.data.frame()
  rownames(counts) <- gene_id$.
  counts$gene <- NULL
  seurat <-
    CreateSeuratObject(counts = counts, meta.data = metadata)
}

seurat <- convert_and_make_seurat(counts)
head(seurat@meta.data)

print(seurat)
# An object of class Seurat
# Unable to map 343 of 33694 requested IDs
# 33351 features across 44024 samples within 1 assay

saveRDS(
  seurat,
  "~/sc_breast/data/masters_project/02_objects_pre_filter/qian_seurat.rds"
)

rm(list = ls())
