# Tabula Sapiens

suppressPackageStartupMessages({
  library(data.table)
  library(readxl)
  library(Seurat)
  library(reticulate)
  library(dplyr)
  library(EnsDb.Hsapiens.v86)
})

sc <- import('scanpy')

adata <- sc$read("/home/gabrielarapozo/sc_breast/data/masters_project/01_download/tabula/tabula.h5ad")

adata_to_seurat <- function(x){
  exprs <- t(adata$X) %>% as.matrix()
  colnames(exprs) <- adata$obs_names$to_list()
  rownames(exprs) <- adata$var_names$to_list()
  seurat <- CreateSeuratObject(exprs)
  seurat <- AddMetaData(seurat, adata$obs)
}

seurat <- adata_to_seurat(adata)
genes <- readr::read_delim("/home/gabrielarapozo/sc_breast/data/masters_project/01_download/tabula/tabula.csv")
genes$ensemblid <- substr(genes$ensemblid,1,15)

mat <- seurat@assays$RNA@counts %>% as.data.frame()
ensembl <-
mapIds(EnsDb.Hsapiens.v86, genes$ensemblid, 'SYMBOL', 'GENEID') %>% as.data.frame() #Unable to map 301 of 29733 requested IDs.
mat$gene <- ensembl$.
mat <- mat %>%
  dplyr::filter(!is.na(gene))
mat <- mat %>% group_by(gene) %>% summarise_each(funs(sum)) %>% as.data.frame()
print("SUM")
rownames(mat) <- mat$gene
mat$gene <- NULL

seurat <- CreateSeuratObject(mat, meta.data = seurat@meta.data)

saveRDS(seurat, "/home/gabrielarapozo/sc_breast/data/masters_project/02_objects_pre_filter/tabula.rds")
