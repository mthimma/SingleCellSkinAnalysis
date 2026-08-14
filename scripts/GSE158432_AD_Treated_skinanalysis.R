# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())

# Download data, rename and organize into folders -------------------------

# Data download from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE158432&format=file

fns <- list.files("data/GSE158432/", full.names = TRUE, pattern = ".gz")

for(fn in fns){
  
  new.fn <- sapply( strsplit(fn, split = "_"), function(x) paste0("data/GSE158432/", x[2], "/", x[3]) )  
  
  dir.create( dirname(new.fn), recursive = TRUE, showWarnings = FALSE )
  
  file.rename(fn, new.fn)
  
}
rm(fn, fns)

# Read in -----------------------------------------------------------------
dirs <- list.files("data/GSE158432/", full.names = TRUE)
names(dirs) <- basename(dirs)

counts <- Read10X(dirs, strip.suffix = TRUE)
dim(counts)  # 33538 26635

raw <- CreateSeuratObject(counts       = counts,
                          min.cells    = 3,
                          min.features = 200)

raw         ## 22176 features across 19379 samples within 1 assay 
rm(counts)

## Add in extra meta study ----
raw$study <- "AD_TmtwithIL-4Rαblockerdupilumab"
raw$condition <- ifelse( grepl("^AD", raw$orig.ident), "AD", "Healthy" )

table(raw$orig.ident)
# AD10 AD11 AD12 AD13 AD14 AD15 AD16 AD17 AD18 AD19 
# 1054 1499 1570 4562 2976  489 1866 2478 1942  943 

## Add in mitochondrial read ----
grep("^MT-", rownames(raw), value = TRUE)  # mitochondrial genes

raw <- PercentageFeatureSet(raw,
                            pattern  = "^Mt-|^mt-|^MT-",
                            col.name = "percent.mt")
png("./results/GSE158432_ViolinPlotofFeatures.png", width=480*2, height = 480*2 )
VlnPlot(raw,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size  = 0, ncol     = 1) &
  labs(x = NULL)
dev.off()
# Standard pipeline -------------------------------------------------------

norm <- raw %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) |> 
  FindNeighbors(verbose = FALSE)
norm@assays$RNA@layers$scale.data <- NULL
## Determine the dimensionality of the data
ElbowPlot(norm, ndims = 50, reduction = "pca")

nPCs <- 30

norm <- RunUMAP(norm,
                dims           = 1:nPCs,
                reduction      = "pca",   #  input embedding name
                reduction.name = "umap")  # output embedding name


rm(raw, nPCs); gc()

# Color the UMAP by sample and status
png("./results/GSE158432_umapcluster_bysamplecondn.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("orig.ident", "condition"),
        reduction = "umap", label=TRUE)
dev.off()

# Color the UMAP by the nFeature_RNA and percent.mt
FeaturePlot(norm, c("nFeature_RNA", "percent.mt"), reduction = "umap") &
  theme_bw() &
  NoGrid()

norm <- FindClusters(norm, resolution = 0.1) #c(0.05, 0.1, 0.2, 0.3))
png("./results/GSE158432_umapcluster.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("RNA_snn_res.0.1", "condition", "seurat_clusters"), label = TRUE)
dev.off()

norm@meta.data |> 
  tabyl(RNA_snn_res.0.1, condition) |> 
  knitr::kable()

# |RNA_snn_res.0.1 |   AD|
#   |:---------------|----:|
#   |0               | 7266| T_cells
#   |1               | 4318| Keratinocytes
#   |2               | 2368| Dendritic Cells
#   |3               | 2200| Keratinocytes, Fibroblasts
#   |4               | 1830| Fibroblasts
#   |5               |  735| T_cells
#   |6               |  366|
#   |7               |  150|
#   |8               |  146| Dendritic Cells

#Annotate the clusters
markers <- list(
  
  Keratinocytes = c(
    "KRT14", "KRT5", "KRT1", "KRT10",
    "KRT16", "KRT17", "FLG", "LOR"
  ),
  
  Fibroblasts = c(
    "COL1A1", "COL1A2", "DCN", "LUM", "PDGFRA"
  ),
  
  Endothelial = c("PECAM1", "VWF", "EMCN", "CDH5", "KDR"),
  
  Lymphatic_Endothelial = c(
    "LYVE1", "PDPN", "FLT4", "CCL21"
  ),
  
  Pericytes = c(
    "RGS5", "CSPG4", "PDGFRB", "MCAM"
  ),
  
  Smooth_Muscle = c(
    "ACTA2", "TAGLN", "MYH11"
  ),
  
  Melanocytes = c(
    "PMEL", "MLANA", "TYR", "DCT", "SOX10"
  ),
  
  Schwann_Cells = c(
    "S100B", "PLP1", "MPZ", "SOX10"
  ),
  
  T_Cells = c(
    "CD3D", "CD3E", "TRBC1", "IL7R"
  ),
  
  NK_Cells = c(
    "NKG7", "GNLY", "PRF1", "KLRD1"
  ),
  
  B_Cells = c(
    "MS4A1", "CD79A", "CD74"
  ),
  
  Plasma_Cells = c(
    "MZB1", "JCHAIN", "XBP1"
  ),
  
  Monocytes = c(
    "LST1", "FCN1", "S100A8", "S100A9"
  ),
  
  Macrophages = c(
    "CD68", "C1QA", "C1QB", "APOE"
  ),
  
  Dendritic_Cells = c(
    "FCER1A", "CD1C", "CLEC10A"
  ),
  
  Langerhans_Cells = c(
    "CD207", "CD1A", "CD74"
  ),
  
  Mast_Cells = c(
    "KIT", "TPSAB1", "CPA3", "HDC"
  )
)

# markers <- readxl::read_excel("data/Lineage_Manager_summary_keep_duplicates.xlsx") |> 
#   separate_longer_delim(Markers, delim = ", ") |> 
#   select(major = Sheet_Name, Markers) |> 
#   distinct()
# 
# markers <- split(markers$Markers, markers$major)

norm <- AddModuleScore_UCell(norm, markers)

cts <- grep("_UCell", colnames(norm@meta.data), v = TRUE)
cts

png("./results/GSE158432_FeaturePlotOfMarker.png", width=480*2, height = 480*2 )
FeaturePlot(norm, features = cts)
dev.off()

#Get number of cells for each cluster grouped by sample and condition
# Assign each cell to the highest-scoring UCell cell type
norm$celltype <- sub("_UCell$", "", cts[max.col(norm@meta.data[, cts], ties.method = "first")])

tab <- norm@meta.data |> 
  dplyr::count(celltype, seurat_clusters, orig.ident, condition) |> 
  unite(group, seurat_clusters, orig.ident, condition, sep = "_") |> 
  pivot_wider(
    names_from = group,
    values_from = n,
    values_fill = 0
  ) 

# Convert numeric part to matrix
mat <- as.matrix(tab[, -1])
rownames(mat) <- tab$celltype
# Add row and column totals
mat <- addmargins(mat)


norm@meta.data |> 
  tabyl(celltype, condition) |> 
  knitr::kable()

# |celltype         |   AD|
#   |:----------------|----:|
#   |B_Cells          |  298|
#   |Dendritic_Cells  | 1442|
#   |Keratinocytes    | 6600|
#   |Langerhans_Cells |  732|
#   |Macrophages      |   94|
#   |Mast_Cells       |   30|
#   |Melanocytes      | 1784|
#   |Monocytes        |  208|
#   |NK_Cells         |  616|
#   |Pericytes        |   26|
#   |Plasma_Cells     |  138|
#   |Schwann_Cells    |   23|
#   |Smooth_Muscle    |   16|
#   |T_Cells          | 7372|
# 
# Write to file
write.table( mat, file = "results/GSE158432_NumofCellsByCelltypeSampleCondition.tsv", sep = "\t",
             quote = FALSE, row.names = TRUE, col.names = TRUE)

saveRDS(norm, file = "data/GSE158432_seuratobj.rds", compress = T)


