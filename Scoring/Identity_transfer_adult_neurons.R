library(Seurat)
library(ggplot2)

ref_obj <- readRDS("/home/giacomo/Documents/Analysis/Striatal_graft/local/data/Adult_astrocytes.rds")

ref_obj@assays$RNA$counts <- ref_obj@assays$RNA$data
ref_obj <- CreateSeuratObject(counts=ref_obj@assays$RNA$counts, meta.data = ref_obj@meta.data)
ref_obj <- NormalizeData(ref_obj)
ref_obj <- FindVariableFeatures(ref_obj)
ref_obj <- ScaleData(ref_obj)
ref_obj <- RunPCA(ref_obj)
ElbowPlot(ref_obj, ndims = 50)

ref_obj <- FindNeighbors(ref_obj, dims = 1:30)
ref_obj <- RunUMAP(ref_obj, dims = 1:30, reduction = "pca")


ref_obj <- RunUMAP(ref_obj, dims = 1:30, reduction = "pca")


query_obj <- readRDS("/home/giacomo/Documents/Analysis/Striatal_graft/dataset/Astrocyte_subset_reannotated.rds")

DefaultAssay(query_obj) <- "RNA"
query_obj <- JoinLayers(query_obj)
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

#query_obj <- FindClusters(query_obj)

#query_obj <- RunUMAP(query_obj, reduction = "integrated.cca", dims = 1:50)
query_obj <- RunUMAP(query_obj, reduction = "integrated.cca", dims = 1:30)
#query_obj <- RunUMAP(query_obj, reduction = "integrated.cca", dims = 1:50)
# query_obj <- RunUMAP(query_obj, reduction = "integrated.cca", dims = 1:30)

DimPlot(query_obj, group.by = c("orig.ident"))



# brain.anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, normalization.method = "LogNormalize",
#                                      reference.reduction = "pca")
query_obj <- JoinLayers(query_obj)
query_obj <- MapQuery(anchorset = brain.anchors, reference = ref_obj, query = query_obj,
                      reference.reduction = "pca", reduction.model = "umap")

query_obj$annotation2 <- "Astrocytes"
p1 <- DimPlot(ref_obj, reduction = "umap", group.by = "tissue", label = TRUE, label.size = 3,
              repel = TRUE) + NoLegend() + ggtitle("Adult Human annotations")
p2 <- DimPlot(query_obj, reduction = "ref.umap", group.by = "annotation_markers2", label = TRUE,
              label.size = 3, repel = TRUE) + NoLegend() + ggtitle("Graft Annotation")
p3 <-p1 + p2
                      

brain.anchors <- FindTransferAnchors(reference = ref_obj, query = query_obj, dims = 1:30, normalization.method = "LogNormalize",
                                     reference.reduction = "pca")
predictions <- TransferData(anchorset = brain.anchors, refdata = ref_obj$supercluster_term, #Once you found the anchors if you want to change metadata to inspect is sufficient to change refdata argument 
                            dims = 1:30)

query_obj <- AddMetaData(query_obj,metadata=predictions)

# query_obj$astro_vs_neuro <- NA
# query_obj@meta.data[which(query_obj@meta.data$labels == "Astrocytes"),]$astro_vs_neuro <- "Astrocytes"
# query_obj@meta.data[which(!(query_obj@meta.data$labels == "Astrocytes")),]$astro_vs_neuro <- "Neurons"

##############################
## Load required libraries  ##
##############################
library(Seurat)
library(dplyr)
library(pheatmap)  # or use ComplexHeatmap, etc.

###########################################
## User-defined parameters and settings  ##
###########################################

#The metadata column that defines groups/subdivisions
group_col <- "annotation_markers2"  # <-- Replace with your metadata column name

# A pattern or exact names for identity score columns

pattern <- "prediction.score."

#The statistic to summarize each group (mean or median)
stat_fun <- mean  # or median


# Extract the object metadata
meta_data <- query_obj@meta.data
meta_data$prediction.score.max <- NULL

# Identify which columns contain your identity scores
# Here we find column names matching 'pattern'
score_cols <- grep(pattern, colnames(meta_data), value = TRUE)

# (Optional) Check that we found the expected columns
print(score_cols)

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


                               
library(pheatmap) 
# order columns and rows and then generate the plot


#my_cols_order <- c("cerebral.nuclei","cerebral.cortex","hippocampal.formation","Deep.layer.near.projecting","Hippocampal.CA1.3","Hippocampal.CA4","Hippocampal.dentate.gyrus",
#"Thalamic.excitatory","Mammillary.body","Midbrain.derived.inhibitory","Cerebellar.inhibitory","Upper.rhombic.lip","Lower.rhombic.lip")

#my_rows_order <- c("MSN", "D1-MSN progenitors", "Emerging MSNs", "CGE interneurons", "MGE-LAMP5 interneurons", "MGE-like interneurons")
#my_rows_order <- c("0","1","2","3","4","5")

my_rows_order <- c("Pre-D1-MSN-like", "Pre-D2-MSN-like", "Emerging MSNs", "MGE-like interneurons", "MGE-LAMP5 interneurons", "CGE interneurons")
#my_rows_order <- c("0","1","2","3","4","5")
my_cols_order <- c("Medium.spiny.neuron","CGE.interneuron","MGE.interneuron","LAMP5.LHX6.and.Chandelier","Eccentric.medium.spiny.neuron",
                  "Amygdala.excitatory","Upper.layer.intratelencephalic","Deep.layer.corticothalamic.and.6b","Deep.layer.intratelencephalic",
                  "Deep.layer.near.projecting","Hippocampal.CA1.3","Hippocampal.CA4","Hippocampal.dentate.gyrus",
                  "Thalamic.excitatory","Mammillary.body","Midbrain.derived.inhibitory","Cerebellar.inhibitory","Upper.rhombic.lip","Lower.rhombic.lip")
colnames <- colnames(agg_scores)
colnames <- gsub(pattern = "prediction.score.", replacement = "", x= colnames)
colnames(agg_scores) <- colnames

# Reorder agg_scores so that rows appear in 'my_row_order'
agg_scores_ordered <- agg_scores[my_rows_order,my_cols_order]

pheatmap(
  agg_scores_ordered,
  scale = "row",               # scales each row to mean=0, sd=1 (optional)
  cluster_rows = TRUE,
  cluster_cols = TRUE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  border_color = NA,
  #filename= "/home/giacomo/Documents/Analysis/Striatal_graft/dataset/Plots/tissue_only_astrocyte_adult_labels_ref_subset_ordered.pdf",
  filename = "E:/All_users/Giacomo/straital_graft/plots/supercluster_term_only_neurons_adult_labels_ref_subset_ordered2.2.pdf",
  width=10,
  fontsize = 17
)

# Turn metadata into a data.frame
df <- query_obj@meta.data
# If your column isn’t already numeric, coerce it
df$score <- as.numeric(df$prediction.score.CGE.interneuron)

# Histogram
p1 <- ggplot(df, aes(x = score)) +
  geom_histogram(binwidth = 0.01, fill = "steelblue", color = "white") +
  theme_classic() +
  labs(
    title = "Histogram of Label-Transfer Scores",
    x     = "Prediction score",
    y     = "Cell count"
  )
ggsave("E:/All_users/Giacomo/straital_graft/plots/CGE_distrib.png", p1)

# Density
ggplot(df, aes(x = score)) +
  geom_density(fill = "lightgreen", alpha = 0.6) +
  theme_classic() +
  labs(
    title = "Density of Label-Transfer Scores",
    x     = "Prediction score",
    y     = "Density"
  )
