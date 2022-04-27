#!/usr/bin/env Rscript

# Esse script baixa os arquivos do Pal et al e adicona metadata

#load packages

library(Seurat)
library(dplyr)
library(readr)
library(magrittr)
library(tidyr)
library(tidyverse)
library(EnsDb.Hsapiens.v86)

# getGEOSuppFiles(GEO =
#   "GSE161529,"
#   makeDirectory = TRUE,
#   baseDir = getwd(),
#   fetch_files = TRUE,
#   filter_regex = NULL
# )
# untar("GSE161529_RAW.tar")

setwd("~/sc_breast/data/masters_project/01_download/pal")
path <- "~/sc_breast/data/masters_project/01_download/pal/GSE161529_RAW/"

metadata <- read_csv("Mestrado - Tabela Suplementar - Controle de Qualidade e Dados da Coorte - Página11 (2).csv")
rownames(metadata) <- metadata$`Sample Name`

colnames(metadata) <- tolower(colnames(metadata))
colnames(metadata) <- gsub(x = colnames(metadata), pattern = " ", replacement = "_")

metadata$patient <- ifelse(metadata$sample_name == 'B1-0023', '0023_p', metadata$patient)
metadata$patient <- ifelse(metadata$sample_name == 'ER-0064-T','0064_t', metadata$patient)
metadata$patient <- ifelse(metadata$sample_name == 'TN-0114-T2','0114_t', metadata$patient)

metadata$molecular_subtype <- ifelse(metadata$molecular_subtype == 'ER+/PR-/PR-', 'ER+', metadata$molecular_subtype)
metadata$molecular_subtype <- ifelse(metadata$molecular_subtype == 'ER+/PR+', 'ER/PR+', metadata$molecular_subtype)

metadata$sample_type <- ifelse(metadata$sample_type == 'Lymph node', 'lymphnode', metadata$sample_type)
metadata$sample_type <- tolower(metadata$sample_type)

metadata$histologycal_subtype <- ifelse(metadata$histologycal_subtype == "Invasive ductal carcinoma NST", "ductal",metadata$histologycal_subtype)
metadata$histologycal_subtype <- ifelse(metadata$histologycal_subtype == "Invasive lobular carcinoma", "lobular", metadata$histologycal_subtype)
metadata$histologycal_subtype <- ifelse(metadata$histologycal_subtype == "Unknown", "unknown",metadata$histologycal_subtype)

metadata$menopause <- tolower(metadata$menopause)
metadata$brca_status <- tolower(metadata$brca_status)
names(metadata)[names(metadata) == 'brca_status'] <- 'BRCA_status'

metadata$study = "Pal et al"
metadata$tech = '10x'
metadata$orig.ident <- metadata$sample_name
rownames(metadata) <- metadata$orig.ident
         
metadata <- subset(metadata, subset = metadata$gender == "F")
metadata$condition <- ifelse(metadata$sample_type == "tumor", paste0("tumor_", metadata$molecular_subtype),metadata$sample_type)


## list files
setwd(path)

files = list.files()
seurat.list <- lapply(files, function(dir) {
  mat <- Read10X(dir, gene.column = 1)
  mat = as(mat, "dgTMatrix")
  print("as dgT")
  mat <- mat %>% as.data.frame()
  ensembl <- mapIds(EnsDb.Hsapiens.v86, rownames(mat),'SYMBOL','GENEID') %>% as.data.frame() 
  print("convert ...")
  mat$gene <- ensembl$.
  mat <- mat %>%
    dplyr::filter(!is.na(gene))
  mat <- mat %>% group_by(gene) %>% summarise_each(funs(sum)) %>% as.data.frame()
  print("SUM")
  rownames(mat) <- mat$gene
  mat$gene <- NULL
  meta <- metadata[metadata$sample_name == dir,]
  meta = cbind(meta, barcodes = colnames(mat)) %>% column_to_rownames("barcodes")
  seurat <-
    CreateSeuratObject(
      counts = mat,
      project = paste0(basename(dir)),
      meta.data = meta,
      min.cells = 3,
      min.features = 100
    )
  print('done seurat')
  return(seurat)
})
names(seurat.list) = files

seurat_B10023 <- seurat.list$`B1-0023`
seurat.list$`B1-0023` <- NULL

seurat.merged <- merge(x = seurat_B10023, y = seurat.list)
print('merge seurat')

seurat.merged[["percent.mt"]] <- PercentageFeatureSet(seurat.merged, pattern = "^MT-")
seurat.merged <- subset(seurat.merged, subset = percent.mt < 15)

print(seurat.merged)


normal <- subset(seurat.merged, subset = condition == 'normal')
saveRDS(normal,
        "~/sc_breast/data/masters_project/02_objects_pre_filter/normal_pal.rds")

list <- c('tumor_ER/PR+', 'tumor_ER+')
er <- subset(seurat.merged, subset = condition %in% list)
saveRDS(er,
        "~/sc_breast/data/masters_project/02_objects_pre_filter/er_pal.rds")

list <- c('lymphnode', 'pre-neoplastic', 'tumor_HER2')
lph <- subset(seurat.merged, subset = condition %in% list)
saveRDS(lph,
        "~/sc_breast/data/masters_project/02_objects_pre_filter/lph_pal.rds")

list <- c('tumor_TNBC', 'tumor_TNBC_BRCA')
tnbc <- subset(seurat.merged, subset = condition %in% list)
saveRDS(tnbc,
        "~/sc_breast/data/masters_project/02_objects_pre_filter/tnbc_pal.rds")


rm(list = ls())
