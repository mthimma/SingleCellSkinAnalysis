# Setup -------------------------------------------------------------------
pacman::p_load(tidyverse, janitor, scDblFinder, patchwork,
               Seurat, UCell, scDblFinder)
setwd("/Users/manjulat/Projects/Araichi/SkinDataAnalysis")
rm(list=ls())


# Download data, rename and organize into folders -------------------------

# Data download from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE147424

mapping_file <- "data/GSE147424/samples_mapping.txt"
map_df <- read.delim(mapping_file, skip =2, header = FALSE, stringsAsFactors = FALSE) |> 
  setNames( c("gsm", "sample_condition") ) |> 
  separate(sample_condition, c("sample", "condition")) |> 
  column_to_rownames("gsm")

read_count_matrix <- function(file) {
  x <- data.table::fread(file, data.table = FALSE) |> 
    column_to_rownames("V1") |> 
    as.matrix()
}

fns <- list.files("data/GSE147424/", full.names = TRUE, pattern = ".txt.gz$")

seurat_list <- map(fns, function(f) {
  
  fname <- basename(f)
  gsm_id <- str_extract(basename(f), "^GSM\\d+")
  cat(gsm_id, "\n")
  
  if (!gsm_id %in% rownames(map_df)) {
    stop("GSM ID not found in mapping file: ", gsm_id)
  }
  
  obj <- CreateSeuratObject(
    counts = read_count_matrix(f),
    project = "GSE147424",
    assay = "RNA",
    min.cells = 3,
    min.features = 200
  )
  
  # Add metadata
  obj$sample    <- map_df[gsm_id, "sample"]
  obj$condition <- map_df[gsm_id, "condition"]
  obj$gsm_id    <- gsm_id

  return(obj)
})

sapply(seurat_list, dim)
combined_seurat <- Reduce(function(x, y) merge(x, y), seurat_list)
combined_seurat <- JoinLayers(combined_seurat)
rm(seurat_list, fns, read_count_matrix, map_df, mapping_file); gc()

combined_seurat # 18,433 genes x 39600 cells

combined_seurat@meta.data |>
  dplyr::count(sample, condition) |> 
  arrange(condition, desc(n))


## Add in mitochondrial read ----
grep("^MT-", rownames(combined_seurat), value = TRUE)  # mitochondrial genes
# No mitochondrial genes found
png("./results/GSE147424_nFeatnCount_Vlnplot.png", width=480*2, height = 480*2 )
VlnPlot(combined_seurat,
        features = c("nFeature_RNA", "nCount_RNA"),
        pt.size  = 0, ncol     = 1) &
  labs(x = NULL)
dev.off()

# Standard pipeline -------------------------------------------------------

norm <- combined_seurat %>%
  NormalizeData(verbose = FALSE) %>%
  FindVariableFeatures(verbose = FALSE) %>%
  ScaleData(verbose = FALSE) %>%
  RunPCA(verbose = FALSE) |> 
  FindNeighbors(verbose = FALSE)

norm@assays$RNA@layers$scale.data <- NULL
norm

rm(combined_seurat)
gc()

## Determine the dimensionality of the data
png("./results/GSE147424_Elbowplot.png", width=480*2, height = 480*2 )
ElbowPlot(norm, ndims = 50, reduction = "pca")
dev.off()

nPCs <- 40

norm <- RunUMAP(norm,
                dims           = 1:nPCs,
                reduction      = "pca",   #  input embedding name
                reduction.name = "umap")  # output embedding name

# Color the UMAP by sample and status
png("./results/GSE147424_umapcluster.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("sample", "condition"), reduction = "umap")
dev.off()
# Color the UMAP by the nFeature_RNA
FeaturePlot(norm, c("nFeature_RNA"), reduction = "umap") &
  theme_bw() &
  NoGrid()

norm <- FindClusters(norm, resolution = 0.2)

png("./results/GSE147424_umapclusterbycondition.png", width=480*2, height = 480*2 )
DimPlot(norm, group.by = c("RNA_snn_res.0.2", "condition"), label = TRUE)
dev.off()

norm@meta.data |> 
  tabyl(RNA_snn_res.0.2, condition) |> 
  knitr::kable()

#   |RNA_snn_res.0.2 |    H|   LS|   NL|
#   |:---------------|----:|----:|----:|
#   |0               | 5577| 1774| 2646|
#   |1               | 4615| 2575| 1803| Fibroblast
#   |2               | 4051|  518|  117| KC
#   |3               | 3821|  550|  159| KC
#   |4               | 1691| 1788|  427| T/NK
#   |5               | 2065|  579|  655|
#   |6               |  557| 1302|  179|
#   |7               |  574|  292|   82|
#   |8               |  383|  153|   59|
#   |9               |  474|   46|   27|
#   |10              |   19|   37|    5|

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

norm@meta.data <- norm@meta.data |> 
  select(-ends_with("UCell"))

norm <- AddModuleScore_UCell(norm, markers)

cts <- grep("_UCell", colnames(norm@meta.data), v = TRUE)
cts

png("./results/GSE147424_FeaturePlotOfMarker.png", width=480*2, height = 480*2 )
FeaturePlot(norm, features = cts)
dev.off()


