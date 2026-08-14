# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())

# Download data, rename and organize into folders -------------------------

# Data download from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE153760&format=file

fns <- list.files("data/GSE153760//rawdata/", full.names = TRUE, pattern = ".gz")

for(fn in fns){
  
  new.fn <- sapply( strsplit(fn, split = "_"), function(x) paste0("data/GSE153760/rawdata/", x[2], "/", x[3]) )  
  
  dir.create( dirname(new.fn), recursive = TRUE, showWarnings = FALSE )
  
  file.rename(fn, new.fn)
  
}
rm(fn, fns)

# Read in -----------------------------------------------------------------
dirs <- list.files("data/GSE153760/rawdata/", full.names = TRUE)
names(dirs) <- basename(dirs)

counts <- Read10X(dirs, strip.suffix = TRUE)
dim(counts)  # 33538 60584

raw <- CreateSeuratObject(counts       = counts,
                          min.cells    = 3,
                          min.features = 200)

raw         ## 24221 features across 57607 samples within 1 assay 
rm(counts)

## Add in extra meta study ----
raw$study <- "ADVsControl_Rojahn_2020"
raw$condition <- ifelse( grepl("^AD", raw$orig.ident), "AD", "Healthy" )

table(raw$orig.ident)

# AD1  AD2  AD3  AD4  AD5  AD6  AD7  AD8  HC1  HC2  HC3  HC4  HC5  HC6  HC7 
# 4846 3803 1819 3222 5628 5104 6246 5555   87  672 2662 2103  668 5455 9737


## Add in mitochondrial read ----
grep("^MT-", rownames(raw), value = TRUE)  # mitochondrial genes

raw <- PercentageFeatureSet(raw,
                            pattern  = "^Mt-|^mt-|^MT-",
                            col.name = "percent.mt")
png("./results/GSE153760_ViolinPlotofFeatures.png", width=480*2, height = 480*2 )
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
png("./results/GSE153760_umapcluster_bysamplecondn.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("orig.ident", "condition", "seurat_clusters"),
        reduction = "umap", label = TRUE)
dev.off()

# Color the UMAP by the nFeature_RNA and percent.mt
FeaturePlot(norm, c("nFeature_RNA", "percent.mt"), reduction = "umap") &
  theme_bw() &
  NoGrid()

norm <- FindClusters(norm, resolution = 0.1) #c(0.05, 0.1, 0.2, 0.3))
png("./results/GSE153760_umapcluster.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("RNA_snn_res.0.1", "condition"), label = TRUE)
dev.off()

norm@meta.data |> 
  tabyl(RNA_snn_res.0.1, condition) |> 
  knitr::kable()

# |RNA_snn_res.0.1 |    AD| Healthy|
#   |:---------------|-----:|-------:|
#   |0               |  9835|    3058| Keratinocytes
#   |1               | 10276|    2524| Dendritic Cells
#   |2               |  6584|    4193| T_cells
#   |3               |  3414|    5138| Keratinocytes
#   |4               |   886|    2811| Fibroblasts
#   |5               |  1480|     776| Keratinocytes
#   |6               |   951|     916|
#   |7               |   730|    1012|
#   |8               |  1367|     309| Mast cells
#   |9               |   700|     647| Fibroblasts


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

png("./results/GSE153760_FeaturePlotOfMarker.png", width=480*2, height = 480*2 )
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

# |celltype              |    AD| Healthy|
#   |:---------------------|-----:|-------:|
#   |B_Cells               |  1254|     801|
#   |Dendritic_Cells       |  4912|    1206|
#   |Endothelial           |   414|     613|
#   |Fibroblasts           |   866|    2780|
#   |Keratinocytes         | 15215|    9066|
#   |Langerhans_Cells      |  1099|     377|
#   |Lymphatic_Endothelial |    20|     147|
#   |Macrophages           |  1715|     618|
#   |Mast_Cells            |  1147|     315|
#   |Melanocytes           |   941|     915|
#   |Monocytes             |  2495|     124|
#   |NK_Cells              |   524|     310|
#   |Pericytes             |    79|      61|
#   |Plasma_Cells          |   221|     264|
#   |Schwann_Cells         |     8|      15|
#   |Smooth_Muscle         |   603|     648|
#   |T_Cells               |  4710|    3124|
#   
#xxxxxxx
# Write to file
write.table( mat, file = "results/GSE153760_NumofCellsByCelltypeSampleCondition.tsv", sep = "\t",
             quote = FALSE, row.names = TRUE, col.names = TRUE)

saveRDS(norm, file = "data/GSE153760_seuratobj.rds", compress = T)

