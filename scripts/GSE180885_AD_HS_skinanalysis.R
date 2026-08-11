# Setup -------------------------------------------------------------------

pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())

# Download data, rename and organize into folders -------------------------

# Data download from https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE180885&format=file

fns <- list.files("data/GSE180885/", full.names = TRUE, pattern = ".gz")

for(fn in fns){

  new.fn <- sapply( strsplit(fn, split = "_"), function(x) paste0("../data/GSE180885/", x[2], "/", x[3]) )  
  
  dir.create( dirname(new.fn), recursive = TRUE, showWarnings = FALSE )
  
  file.rename(fn, new.fn)
  
}
rm(fn, fns)

# Read in -----------------------------------------------------------------
dirs <- list.files("data/GSE180885/", full.names = TRUE)
names(dirs) <- basename(dirs)

counts <- Read10X(dirs, strip.suffix = TRUE)
dim(counts)  # 33538 26278

raw <- CreateSeuratObject(counts       = counts,
                          min.cells    = 3,
                          min.features = 200)

raw         # 19195 features across 25034 samples within 1 assay 
rm(counts)

## Add in extra meta study ----
raw$study <- "Alkon_JACI2022"
raw$condition <- ifelse( grepl("^AD", raw$orig.ident), "AD", "Healthy" )

table(raw$orig.ident)

## Add in mitochondrial read ----
grep("^MT-", rownames(raw), value = TRUE)  # mitochondrial genes

raw <- PercentageFeatureSet(raw,
                            pattern  = "^Mt-|^mt-|^MT-",
                            col.name = "percent.mt")

VlnPlot(raw,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
        pt.size  = 0, ncol     = 1) &
  labs(x = NULL)

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
png("./results/GSE180885_umapcluster_bysamplecondn.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("orig.ident", "condition"),
        reduction = "umap")
dev.off()

# Color the UMAP by the nFeature_RNA and percent.mt
FeaturePlot(norm, c("nFeature_RNA", "percent.mt"), reduction = "umap") &
  theme_bw() &
  NoGrid()

norm <- FindClusters(norm, resolution = c(0.05, 0.1, 0.2, 0.3))
png("./results/GSE180885_umapcluster.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("RNA_snn_res.0.05", "condition"), label = TRUE)
dev.off()

norm@meta.data |> 
  tabyl(RNA_snn_res.0.05, condition) |> 
  knitr::kable()

#   |RNA_snn_res.0.05 |   AD | Healthy|
#   |:----------------|-----:|-------:|
#   |0                | 10445|      95| KC
#   |1                |  1189|    5817|T/NK
#   |2                |   646|    3020|T/NK
#   |3                |     3|    1905|KC
#   |4                |   784|       3|Fibrobloast
#   |5                |   103|     326|
#   |6                |   334|       3|
#   |7                |   284|       1|Fibrobloast
#   |8                |    71|       5|

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

png("./results/GSE180885_FeaturePlotOfMarker.png", width=480*2, height = 480*2 )
FeaturePlot(norm, features = cts)
dev.off()



