# Find genes shared in all four studies of Skin analysis

pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())

seu1 <- readRDS('./data/GSE153760_seuratobj.rds')
seu2 <- readRDS('./data/GSE147424_seuratobj.rds')
seu3 <- readRDS('./data/GSE162054_seuratobj.rds')
seu4 <- readRDS('./data/GSE180885_seuratobj.rds')

seurat_list <- list(seu1, seu2, seu3, seu4)
shared_genes <- Reduce(intersect, lapply(seurat_list, rownames))
write.table(shared_genes, "./data/shared_genes_across_studies.txt",
            row.names = FALSE, col.names = FALSE, quote = FALSE)
#subset of seurat object for shared genes
seu1_sub <- subset(seu1, features = shared_genes)
seu2_sub <- subset(seu2, features = shared_genes)
seu3_sub <- subset(seu3, features = shared_genes)
seu4_sub <- subset(seu4, features = shared_genes)

length(shared_genes)
class(shared_genes)
