#!/usr/bin/env Rscript

#Download https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-download

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratObject)
  library(dplyr)
  library(magrittr)
  library(tidyr)
  library(data.table)
  library(GEOquery)
  library(EnsDb.Hsapiens.v86)
})

#setwd
setwd("~/sc_breast/data/masters_project/01_download/wu/")

metadata <- read.csv("Whole_miniatlas_meta.csv") #sup1
metadata <- metadata[-1,]
rownames(metadata) <- metadata$NAME
names(metadata)[names(metadata) == 'Patient'] <- 'patient'

supp <- readxl::read_excel("41588_2021_911_MOESM4_ESM.xlsx", col_names = FALSE)  #sup2
supp <- supp[-(1:3),]
colnames(supp) <- supp[1,]
supp <- supp[-1,]

names(supp)[names(supp) == 'Case ID'] <- 'patient'
supp$patient <- (gsub("^", "CID", supp$patient))
supp$patient <- {gsub("CID4290", "CID4290A", gsub("CID4404-1","CID44041", gsub("CID4497-1", "CID44971", gsub("CID4499-1", "CID44991", gsub("CID4517-1", "CID45171", gsub("CID4530", "CID4530N",supp$patient))))))}
supp$molecular_subtype = c('ER/PR/HER2', rep('HER2',2), 'ER/PR+', 'TNBC', 'ER/PR+', 'ER+', 'ER/PR+', 'ER/HER2+', rep('ER/PR+',3), 'TNBC', 'ER+','ER/PR+', 'TNBC',  'ER/PR+', rep("TNBC",2), 'TNBC_BRCA', rep('TNBC',2),'HER2', 'TNBC', rep('ER/PR+',2))
supp <- supp[ ,-(5:15)]
colnames(supp) <- tolower(colnames(supp))
         
meta <- inner_join(metadata, supp, "patient") #merge sup1 and sup3
rownames(meta) <- meta$NAME

meta = meta %>%
  ## Conditions to the integrated data:
  mutate(
    patient = paste0("Wu_", meta$patient),
    study = "Wu et al",
    tech = "10x",
    tissue_origin = 'tumor',
    sample_type = "tumor",
    treatment_group = "treatment_naive",
    histologycal_subtype <- ifelse(meta$patient == 'Wu_CID4471'|meta$patient == 'Wu_CID4535', "Lobular", "Ductal"),
    tumor_site = "primary",
    menopause <- "unknown",
    BRCA_status <- ifelse(meta$patient == "Wu_CID44991", "positive", "unknown")
  )
meta$condition <- ifelse(meta$sample_type == "tumor", paste0("tumor_", meta$molecular_subtype),meta$condition)

#### Download data ----
# run once
# https://singlecell.broadinstitute.org/single_cell/study/SCP1039/a-single-cell-and-spatially-resolved-atlas-of-human-breast-cancers#study-download

# getGEOSuppFiles(
#   GEO =
#     "GSE176078",
#   makeDirectory = TRUE,
#   baseDir = getwd(),
#   fetch_files = TRUE,
#   filter_regex = NULL
# )
# setwd("GSE176078/")
# untar("GSE176078_RAW.tar")
# untar("GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz")

setwd("GSE176078/")

count <- Read10X("Wu_etal_2021_BRCA_scRNASeq/", gene.column = 1)

convert_and_make_seurat <- function(x) {
  ensembl <-
    mapIds(EnsDb.Hsapiens.v86, rownames(count), 'GENEID', 'SYMBOL') %>% as.data.frame() 
  dgT_matrix = as(count, "dgTMatrix")
  count <- dgT_matrix %>% as.data.frame()
  count$gene <- ensembl$.
  count <- count %>%
    dplyr::filter(!is.na(gene))
  gene_id <-
    mapIds(EnsDb.Hsapiens.v86, count$gene, 'SYMBOL', 'GENEID') %>% as.data.frame()
  rownames(count) <- gene_id$.
  count$gene <- NULL
  seurat <-
    CreateSeuratObject(counts = count, meta.data = meta)
    #remove treatment patients
  '%ni%' <- Negate("%in%")
  list <- c('Wu_CID3963', "Wu_CID4066", "Wu_CID4398", "Wu_CID4513", "Wu_CID4523")
  seurat <- subset(seurat, subset = patient %ni% list)
}

seurat <- convert_and_make_seurat(count)
head(seurat@meta.data)

print(seurat)
# An object of class Seurat
# Unable to map 301 of 29733 requested IDs.
# 29432 features across 79404 samples within 1 assay 

saveRDS(seurat,
        "~/sc_breast/data/masters_project/02_objects_pre_filter/wu.rds")

rm(list = ls())
