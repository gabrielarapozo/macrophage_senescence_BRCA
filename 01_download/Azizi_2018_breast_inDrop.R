#### Azizi 2018 ----

#download Azizi et al data, add metadata, and harmonizes the gene IDs

# load packages
suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(Seurat)
  library(GEOquery)
  library(EnsDb.Hsapiens.v86)
})

#### Download data ----
# run once

# setwd("~/sc_breast/data/masters_project/01_download/azizi/")
# getGEOSuppFiles("GSE114725",  makeDirectory = TRUE,  baseDir = getwd(),  fetch_files = TRUE,  filter_regex = NULL)
  setwd("~/sc_breast/data/masters_project/01_download/azizi/GSE114725/")
# gunzip("GSE114725_rna_raw.csv.gz")

matrix <- fread("GSE114725_rna_raw.csv")

#### prepare count matrix ----
counts = matrix %>% mutate(barcode=paste0("barcode_",cellid)) %>%
  column_to_rownames("barcode") %>%
  dplyr::select(6:ncol(matrix)) %>% t()

#### Prepare metadata Azizi ----
# Figure S1: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6348010/#S14
metadata = matrix %>% dplyr::select(1:5) %>%
  mutate(barcode=paste0("barcode_",cellid)) %>% dplyr::select(c(1:3,6))

meta.paper = tibble(patient=c(paste0("BC",1:8)), histologycal_subtype=c("Ductal"), grade=c(1,2,3,1,3,2,3,2),
                    n_metastasis=c(0,1,0,1,0,0,0,0),
                    molecular_subtype=c("ER/PR+","ER+","TNBC","ER/PR+","TNBC","ER+","HER2","TNBC"),
                    age = c(38, 60, 43, 52, 78, 58, 65,72),
                    menopause = c('Pre', 'Post', 'Pre', 'Pre', 'Post', 'Post', 'Post', 'Post'),
                    BRCA_status = c("Negative", "Unknown", "Negative", "Negative", "Unknown", "Unknown", "Unknown", "Unknown"),
                    gender = "F")

metadata = inner_join(metadata, meta.paper)

metadata = metadata %>%
  mutate(study= 'Azizi et al',
         tech = 'InDrop',
         tumor_type ="breast",
         tech="inDrop",
         sample_type = tolower(metadata$tissue),
         study="Azizi et al",
         tumor_site=ifelse(tissue=="TUMOR","primary",NA)) %>% 
      column_to_rownames("barcode")
metadata$replicate <- NULL
metadata$tissue <- NULL
metadata$patient <- gsub("BC", "Azizi_", metadata$patient)

metadata$condition = ifelse(metadata$sample_type == "blood", "blood", NA)
metadata$condition = ifelse(metadata$sample_type == "lymphnode", "lymphnode", metadata$condition)
metadata$condition = ifelse(metadata$sample_type == "normal", "normal", metadata$condition)
metadata$condition = ifelse(metadata$sample_type == "tumor", "tumor", metadata$condition)
metadata$condition = ifelse(metadata$condition == "tumor", paste0("tumor_", metadata$molecular_subtype),metadata$condition)

head(metadata)

#### harmonize gene symbol and create object ----

convert_and_make_seurat <- function(x){
  ensembl <- mapIds(EnsDb.Hsapiens.v86, rownames(counts),'GENEID','SYMBOL') %>% as.data.frame() 
  counts <- counts %>% as.data.frame() 
  counts$gene <- ensembl$.
  counts <- counts %>%
  dplyr::filter(!is.na(gene))
  gene_id <- mapIds(EnsDb.Hsapiens.v86, counts$gene,'SYMBOL','GENEID') %>% as.data.frame()
  rownames(counts) <- gene_id$.
  counts$gene <- NULL
  seurat <- CreateSeuratObject(counts = counts, meta.data = metadata)
}

seurat <- convert_and_make_seurat(counts)
head(seurat@meta.data)

print(seurat)
# An object of class Seurat 
#Unable to map 274 of 14875 requested IDs.
# 14601 features across 47016 samples within 1 assay 

saveRDS(seurat, "~/sc_breast/data/masters_project/02_objects_pre_filter/azizi_seurat.rds")

rm(list = ls()) 
