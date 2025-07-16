library(Seurat)
library(ggplot2)


ref_obj <- readRDS("E:/All_users/Giacomo/straital_graft/data/Adult_astrocytes.rds")#
#ref_obj <- readRDS("/home/giacomo/Documents/Analysis/Striatal_graft/local/data/Adult_astrocytes.rds")


ref_obj <- CreateSeuratObject(counts=ref_obj@assays$RNA$data, meta.data = ref_obj@meta.data)


ref_obj <- subset(ref_obj, subset = development_stage == "29-year-old stage")


ref_obj$dataset <- "Siletti"

DefaultAssay(ref_obj) <- "RNA"
ref_obj <- JoinLayers(ref_obj)


query_obj <- readRDS("E:/All_users/Giacomo/straital_graft/data/Astrocyte_subset_reannotated.rds")

query_obj$dataset <- query_obj$orig.ident

DefaultAssay(query_obj) <- "RNA"

query_obj <- JoinLayers(query_obj)


merged_obj <- merge(ref_obj, query_obj)

merged_obj <- JoinLayers(merged_obj)
merged_obj <- CreateSeuratObject(counts=merged_obj@assays$RNA$counts, meta.data = merged_obj@meta.data)

merged_obj[["RNA"]] <- split(merged_obj[["RNA"]], f = merged_obj$dataset)


merged_obj <- NormalizeData(merged_obj)
merged_obj <- FindVariableFeatures(merged_obj)
merged_obj <- ScaleData(merged_obj)
merged_obj <- RunPCA(merged_obj)
ElbowPlot(merged_obj, ndims = 50)
merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
# merged_obj <- FindNeighbors(merged_obj, dims = 1:30)
#merged_obj <- FindClusters(merged_obj)


merged_obj <- RunUMAP(merged_obj, dims = 1:30)
# merged_obj <- RunUMAP(merged_obj, dims = 1:30)


DimPlot(merged_obj, group.by = c("dataset"))

merged_obj <- IntegrateLayers(object = merged_obj, method = CCAIntegration, orig.reduction = "pca",
                             new.reduction = "integrated.cca", verbose = FALSE)
merged_obj <- FindNeighbors(merged_obj, reduction = "integrated.cca", dims = 1:30)

merged_obj <- RunUMAP(merged_obj, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

DimPlot(merged_obj, group.by = c("dataset"), label= T, reduction="umap.cca")

DimPlot(merged_obj, group.by = c("tissue"), label= T, reduction="umap.cca")

merged_obj$annotation <- merged_obj$tissue
merged_obj@meta.data[which(!(merged_obj@meta.data$dataset == "Siletti")),]$annotation <- "graft-derived astrocytes"




brain.anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, normalization.method = "LogNormalize",
                                     reference.reduction = "pca")
predictions <- TransferData(anchorset = brain.anchors, refdata = ref_obj$Subregion, 
                            dims = 1:30)

query_obj <- AddMetaData(query_obj,metadata=predictions)

library(pheatmap)  
group_col <- "integrated_clusters_k0.2"

pattern <- "prediction.score."

# The statistic to summarize each group (mean or median)
stat_fun <- mean  # or median


# Extract the object metadata
meta_data <- query_obj@meta.data
meta_data$prediction.score.max <- NULL

# Identify which columns contain your identity scores
score_cols <- grep(pattern, colnames(meta_data), value = TRUE)


# Subset the metadata to include only the grouping column + the score columns
sub_meta <- meta_data[, c(group_col, score_cols)]


# Summarize score values by group
agg_scores <- sub_meta %>%
  group_by(.data[[group_col]]) %>%
  summarize(across(all_of(score_cols), stat_fun)) %>%
  as.data.frame()

# Move the grouping column into row names
rownames(agg_scores) <- agg_scores[[group_col]]
agg_scores[[group_col]] <- NULL



my_rows_order <- c("0","1","2","3","4","5")

my_cols_order <- c("cerebral.nuclei","cerebral.cortex","hippocampal.formation","Deep.layer.near.projecting","Hippocampal.CA1.3","Hippocampal.CA4","Hippocampal.dentate.gyrus",
                   "Thalamic.excitatory","Mammillary.body","Midbrain.derived.inhibitory","Cerebellar.inhibitory","Upper.rhombic.lip","Lower.rhombic.lip")


# Reorder agg_scores so that rows appear in 'my_row_order'
agg_scores_ordered <- agg_scores[my_rows_order,my_cols_order]

pheatmap(
  agg_scores,
  scale = "row",               
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  filename = "E:/All_users/Giacomo/straital_graft/plots/Subregion_heatmap_astro_clusters.pdf",
  width=10
)
=======
DimPlot(merged_obj, group.by = c("annotation"), label= T, reduction="umap.cca", pt.size=1, raster=T)

