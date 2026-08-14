# Find genes shared in all four studies of Skin analysis

pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder, UpSetR)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())

seu_list <- list(
  Rojahn_2020   = readRDS('./data/GSE153760_seuratobj.rds'),
  He_2020       = readRDS('./data/GSE147424_seuratobj.rds'),
  Rindler_2021  = readRDS('./data/GSE162054_seuratobj.rds'),
  Alkon_2022    = readRDS('./data/GSE180885_seuratobj.rds'),
  Bangert_2021  = readRDS('./data/GSE158432_seuratobj.rds')
)

sapply(seu_list, ncol)  # number of cells in each study
# Rojahn_2020      He_2020 Rindler_2021   Alkon_2022 Bangert_2021 
# 57607        39600         2637        25034        19379 

sapply(seu_list, nrow)  # number of genes in each study
# Rojahn_2020      He_2020 Rindler_2021   Alkon_2022 Bangert_2021 
# 24221        18433        18403        22108        22176 

seu_list$Rojahn_2020@meta.data |> 
  dplyr::count(orig.ident) |> 
  arrange(orig.ident) 
# Overlap of genes expressed in at least 1% of the cells -------------------
lapply(seu_list, rownames) |> unlist() |> unique() |> length()  # 26,933

expressed <- lapply(seu_list, 
                    function(obj) which( rowMeans( GetAssayData(obj, layer = "data") > 0 ) > 0.01 ) |> names())

sapply(expressed, length)
# Rojahn_2020      He_2020 Rindler_2021   Alkon_2022 Bangert_2021 
# 12396        11572        13378        12149        13359 

expressed |> unlist() |> unique() |> length() #15,055
fromList(expressed) |> upset(order.by = "freq")

tb <- expressed |> unlist() |> table()
keep <- which(tb >= 3)  |> names() #  9982 + 2100 = (9982 + 1646 +258 +183 +13)

seu_list_filt <- lapply(seu_list, function(obj) subset(obj, features = keep))

