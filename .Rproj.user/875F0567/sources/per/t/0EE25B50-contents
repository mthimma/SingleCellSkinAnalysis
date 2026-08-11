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
dim(counts)  # 33538 87219

raw <- CreateSeuratObject(counts       = counts,
                          min.cells    = 3,
                          min.features = 200)

raw         ## 24984 features across 76986 samples within 1 assay 
rm(counts)

## Add in extra meta study ----
raw$study <- "ADVsControl_TmtwithIL-4Rαblockerdupilumab"
raw$condition <- ifelse( grepl("^AD", raw$orig.ident), "AD", "Healthy" )

table(raw$orig.ident)
#AD1 AD10 AD11 AD12 AD13 AD14 AD15 AD16 AD17 AD18 AD19  AD2  AD3  AD4  AD5  AD6  AD7  AD8  HC1  HC2  HC3  HC4  HC5 
# 4846 1054 1499 1570 4562 2976  489 1866 2478 1942  943 3803 1819 3222 5628 5104 6246 5555   87  672 2662 2103  668 
# HC6  HC7 
# 5455 9737 

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
DimPlot(norm, group.by = c("orig.ident", "condition"),
        reduction = "umap")
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
#   |0               | 14918|    4136| T/NK
#   |1               | 12762|    2522|
#   |2               | 10042|    4373| Keratinocytes
#   |3               |  8716|    3848| Keratinocytes
#   |4               |  2793|     915|
#   |5               |   878|    2801| Fibro
#   |6               |  2512|     809| Keratinocytes
#   |7               |  1537|     313|
#   |8               |   733|    1015|
#   |9               |   711|     652| Fibro

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



