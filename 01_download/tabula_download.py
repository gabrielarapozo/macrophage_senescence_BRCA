#!/usr/bin/env python
# coding: utf-8

# # Tabula Sapiens - mammary

import scanpy as sc

#https://figshare.com/articles/dataset/Tabula_Sapiens_release_1_0/14267219?file=28902069
adata = sc.read("/home/gabrielarapozo/sc_breast/data/masters_project/01_download/tabula/TS_Mammary.h5ad")

df = adata.var
df.to_csv("/home/gabrielarapozo/sc_breast/data/masters_project/01_download/tabula/tabula.csv")

adata.write("/home/gabrielarapozo/sc_breast/data/masters_project/01_download/tabula/tabula.h5ad")
