library(Seurat)
library(ggplot2)
library(schard)

ref_obj <- h5ad2seurat('E:/All_users/Giacomo/straital_graft/data/Adult_Brain_Neurons_29yo.h5ad')

ref_obj@assays$RNA$counts <- ref_obj@assays$RNA$data
ref_obj <- CreateSeuratObject(counts=ref_obj@assays$RNA$counts, meta.data = ref_obj@meta.data)
ref_obj <- FindVariableFeatures(ref_obj)
ref_obj <- ScaleData(ref_obj)
ref_obj <- RunPCA(ref_obj)
ElbowPlot(ref_obj, ndims = 50)

ref_obj <- FindNeighbors(ref_obj, dims = 1:30)
ref_obj <- RunUMAP(ref_obj, dims = 1:30, reduction = "pca")


ref_obj <- RunUMAP(ref_obj, dims = 1:30, reduction = "pca")

query_obj <- subset(query_obj, subset= annotation_markers2 %in% c("CGE interneurons","Emerging MSNs", "MGE-like interneurons", "MGE-LAMP5 interneurons", "Pre-D1-MSN-like", "Pre-D2-MSN-like"))

DefaultAssay(query_obj) <- "RNA"
query_obj[["RNA"]] <- split(query_obj[["RNA"]], f = query_obj$orig.ident)

query_obj <- NormalizeData(query_obj)
query_obj <- FindVariableFeatures(query_obj)
query_obj <- ScaleData(query_obj)
query_obj <- RunPCA(query_obj)
ElbowPlot(query_obj, ndims = 50)
query_obj <- FindNeighbors(query_obj, dims = 1:30)
query_obj <- FindClusters(query_obj)


query_obj <- RunUMAP(query_obj, dims = 1:30)


DimPlot(query_obj, group.by = c("orig.ident"))

query_obj <- IntegrateLayers(object = query_obj, method = CCAIntegration, orig.reduction = "pca",
                             new.reduction = "integrated.cca", verbose = FALSE)
query_obj <- FindNeighbors(query_obj, reduction = "integrated.cca", dims = 1:30)

query_obj <- RunUMAP(query_obj, reduction = "integrated.cca", dims = 1:30)


brain.anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, normalization.method = "LogNormalize",
                                     reference.reduction = "pca")
predictions <- TransferData(anchorset = brain.anchors, refdata = ref_obj$Supercluster_term, 
                            dims = 1:30)

query_obj <- AddMetaData(query_obj,metadata=predictions)


library(pheatmap)  
group_col <- "annotation_markers2"  

# Exact names for identity score columns

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



my_cols_order <- c("Medium.spiny.neuron","CGE.interneuron","MGE.interneuron","LAMP5.LHX6.and.Chandelier","Eccentric.medium.spiny.neuron",
                   "Amygdala.excitatory","Upper.layer.intratelencephalic","Deep.layer.corticothalamic.and.6b","Deep.layer.intratelencephalic",
                   "Deep.layer.near.projecting","Hippocampal.CA1.3","Hippocampal.CA4","Hippocampal.dentate.gyrus",
                   "Thalamic.excitatory","Mammillary.body","Midbrain.derived.inhibitory","Cerebellar.inhibitory","Upper.rhombic.lip","Lower.rhombic.lip")

my_rows_order <- c("Pre-D1-MSN-like", "Pre-D2-MSN-like", "Emerging MSNs", "MGE-like interneurons", "MGE-LAMP5 interneurons", "CGE interneurons")

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
  filename = "E:/All_users/Giacomo/straital_graft/plots/Supercluster_term_heatmap_no_splatter_neuro_annotation.pdf",
  width=10
)
