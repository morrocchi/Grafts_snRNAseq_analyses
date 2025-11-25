# Set up the environment ----
library(ggplot2) 
library(tidyverse) 
library(dplyr) 
library(RColorBrewer)
library(Seurat)
library(cowplot)
library(SeuratData)
library(org.Mm.eg.db, verbose = FALSE)
library(dittoSeq)
library(readxl)

options(future.globals.maxSize = 10000000000)


setwd("~/Documents/Analysis/Striatal_graft/dataset")


human_hMSN_6M <- readRDS(file = "/home/giacomo/Documents/Analysis/Striatal_graft/local/data/Final_ENSEMBL_GENES_grafts_Conforti_markers_annotation.rds")

DefaultAssay(human_hMSN_6M) <- "RNA"

human_hMSN_6M<-JoinLayers(human_hMSN_6M)

human_hMSN_6M<-CreateSeuratObject(counts = human_hMSN_6M@assays$RNA$counts ,meta.data = human_hMSN_6M@meta.data)

human_hMSN_6M[["RNA"]] <- split(human_hMSN_6M[["RNA"]], f = human_hMSN_6M$orig.ident)

human_hMSN_6M<-SCTransform(human_hMSN_6M,vars.to.regress=c("nFeature_RNA"))
human_hMSN_6M <- RunPCA(human_hMSN_6M, dims = 1:30)
human_hMSN_6M <- FindNeighbors(human_hMSN_6M, dims=1:30)
# integrate the datasets
human_hMSN_6M <- IntegrateLayers(human_hMSN_6M, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                 verbose = F, normalization.method = "SCT")

# cluster the integrated data
human_hMSN_6M <- FindNeighbors(human_hMSN_6M, reduction = "integrated.cca", dims = 1:30)
human_hMSN_6M <- RunUMAP(human_hMSN_6M, reduction.name = "umap.cca", dims=1:30)
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 1, cluster.name = "integrated_clusters_k1")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 0.5, cluster.name = "integrated_clusters_k0.5")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 0.45, cluster.name = "integrated_clusters_k0.45")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 0.35, cluster.name = "integrated_clusters_k0.35")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 0.3, cluster.name = "integrated_clusters_k0.3")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 1.5, cluster.name = "integrated_clusters_k1.5")
human_hMSN_6M <- FindClusters(human_hMSN_6M, resolution = 0.1, cluster.name = "integrated_clusters_k0.1")

human_hMSN_6M <- RunUMAP(human_hMSN_6M, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")


# Subset only cells that are likely to be astrocytes based on the markers of Sofroniew ----
Sofroniew_2022_astro_genes <- read_excel("/home/giacomo/Documents/Analysis/Striatal_graft/local/Sofroniew 2022_HealtyReactiveEMTAstro.xlsx", 
                                         sheet = "AstroGenePanel")

Sofroniew_2022_astro_genes <- toupper(Sofroniew_2022_astro_genes$`Gene Astro`)
#Sofroniew_2022_astro_genes <- na.omit(Sofroniew_2022_astro_genes)

# Convert SYMBOLS to ENSEMBL ids ----

library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

df <- data.frame(gene_symbol = Sofroniew_2022_astro_genes)

conversion_result <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = df$gene_symbol,
  mart       = mart
)

df_converted <- merge(
  x  = df,
  y  = conversion_result,
  by.x = "gene_symbol",
  by.y = "hgnc_symbol",
  all.x = TRUE
)



not_found <- Sofroniew_2022_astro_genes[which(!(Sofroniew_2022_astro_genes %in% conversion_result$hgnc_symbol))]
ENS <- c("ENSG00000099194","ENSG00000197238","ENSG00000143514","ENSG00000122378","ENSG00000111271",
         "ENSG00000132613","ENSG00000142871","ENSG00000180573","ENSG00000145476","ENSG00000172640",
         "ENSG00000133048","ENSG00000241186","ENSG00000109062","ENSG00000176076","ENSG00000186204",
         "ENSG00000164976","ENSG00000134716","ENSG00000187773","ENSG00000125148","ENSG00000171903",
         "ENSG00000099194","ENSG00000074800","ENSG00000182329","ENSG00000205362","ENSG00000227742",
         "ENSG00000227742","ENSG00000133816","ENSG00000169239","ENSG00000205702")

Sofroniew_2022_astro_genes <- df_converted[!(is.na(df_converted$ensembl_gene_id)),]$ensembl_gene_id

Sofroniew_2022_astro_genes <- unique(c(Sofroniew_2022_astro_genes,ENS))




genes.list <- list(Sofroniew_2022_astro_genes)
enrich.name <- "Astrocyte_score"
human_hMSN_6M <- AddModuleScore(human_hMSN_6M,
                                features = genes.list,
                                pool = NULL,
                                n.bin = 5,
                                seed = 1,
                                ctrl = length(genes.list),
                                k = FALSE,
                                name = enrich.name,
                                random.seed = 1)



#save900x800
FeaturePlot(human_hMSN_6M, features = "Astrocyte_score1", reduction = "umap.cca", pt.size = 1, max.cutoff="q99",min.cutoff="q1", order=T) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + theme(axis.title.x = element_blank(),
                                                                                                                                    axis.text.x = element_blank(),
                                                                                                                                    axis.ticks.x = element_blank(),
                                                                                                                                    axis.title.y = element_blank(),
                                                                                                                                    axis.text.y = element_blank(),
                                                                                                                                    axis.ticks.y = element_blank(),
                                                                                                                                    axis.line = element_blank())


Idents(human_hMSN_6M) <- human_hMSN_6M$integrated_clusters_k0.5
Astrocyte_subset <- subset(human_hMSN_6M, idents=c("0", "1", "2", "4", "7"))


DefaultAssay(Astrocyte_subset) <- "RNA"

Astrocyte_subset<-JoinLayers(Astrocyte_subset)

Astrocyte_subset[["RNA"]] <- split(Astrocyte_subset[["RNA"]], f = Astrocyte_subset$orig.ident)

Astrocyte_subset<-SCTransform(Astrocyte_subset,vars.to.regress=c("nFeature_RNA"))
Astrocyte_subset <- RunPCA(Astrocyte_subset, dims = 1:30)
Astrocyte_subset <- FindNeighbors(Astrocyte_subset, reduction = "pca", dims = 1:30)

# integrate the datasets
Astrocyte_subset <- IntegrateLayers(Astrocyte_subset, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca",
                                    verbose = F, normalization.method = "SCT")



# cluster the integrated data
Astrocyte_subset <- FindNeighbors(Astrocyte_subset, reduction = "integrated.cca", dims = 1:30)
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 1, cluster.name = "integrated_clusters_k1")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 0.5, cluster.name = "integrated_clusters_k0.5")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 0.4, cluster.name = "integrated_clusters_k0.4")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 1.5, cluster.name = "integrated_clusters_k1.5")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 0.1, cluster.name = "integrated_clusters_k0.1")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 0.2, cluster.name = "integrated_clusters_k0.2")
Astrocyte_subset <- FindClusters(Astrocyte_subset, resolution = 0.3, cluster.name = "integrated_clusters_k0.3")

Astrocyte_subset <- RunUMAP(Astrocyte_subset, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")

dittoDimPlot(Astrocyte_subset, reduction.use = "umap.cca", var = "integrated_clusters_k0.2", main = NULL, size = 1.5,
             do.label = T, labels.size = 6, labels.highlight = F, opacity = 1, xlab = "umap1", ylab = "umap2") + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.line = element_blank())

dittoDimPlot(Astrocyte_subset, reduction.use = "umap.cca", var = "SCT_snn_res.0.45", main = NULL, size = 1.5,
             do.label = T, labels.size = 6, labels.highlight = F, opacity = 1, xlab = "umap1", ylab = "umap2") + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.line = element_blank())

saveRDS(Astrocyte_subset, file = "Astrocyte_subset_reannotated.rds")

dittoDimPlot(Astrocyte_subset, reduction.use = "harmony", var = "SCT_snn_res.0.45", main = NULL, size = 1.5,
             do.label = T, labels.size = 6, labels.highlight = F, opacity = 1, xlab = "umap1", ylab = "umap2") + 
  theme_classic() + theme(legend.position = "none", axis.title.x = element_blank(),
                          axis.text.x = element_blank(),
                          axis.ticks.x = element_blank(),
                          axis.title.y = element_blank(),
                          axis.text.y = element_blank(),
                          axis.ticks.y = element_blank(),
                          axis.line = element_blank())


Sofroniew_2022_ReactivityHealthy_genes <- read_excel("/home/giacomo/Documents/Analysis/Striatal_graft/local/Sofroniew 2022_HealtyReactiveEMTAstro.xlsx", 
                                                     sheet = "HealthyReactiveEMT")

Sofroniew_2022_reactivity_genes <- Sofroniew_2022_ReactivityHealthy_genes[,2]
Sofroniew_2022_reactivity_genes <- toupper(Sofroniew_2022_reactivity_genes$`Reactivity list`)
Sofroniew_2022_reactivity_genes <- na.omit(Sofroniew_2022_reactivity_genes)

Sofroniew_2022_EMT_genes <- Sofroniew_2022_ReactivityHealthy_genes[,3]
Sofroniew_2022_EMT_genes <- toupper(Sofroniew_2022_EMT_genes$`EMT list`)
Sofroniew_2022_EMT_genes <- na.omit(Sofroniew_2022_EMT_genes)


library(biomaRt)

mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

df_reactivity <- data.frame(gene_symbol = Sofroniew_2022_reactivity_genes)

conversion_result_reactivity <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = df_reactivity$gene_symbol,
  mart       = mart
)

df_converted_reactivity <- merge(
  x  = df_reactivity,
  y  = conversion_result_reactivity,
  by.x = "gene_symbol",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

reactivity_genes <- df_converted_reactivity[!(is.na(df_converted_reactivity$ensembl_gene_id)),"ensembl_gene_id"]

genes.list <- list(reactivity_genes)
enrich.name <- "Reactivity_score"
Astrocyte_subset <- AddModuleScore(Astrocyte_subset,
                                   features = genes.list,
                                   pool = NULL,
                                   n.bin = 5,
                                   seed = 1,
                                   ctrl = length(genes.list),
                                   k = FALSE,
                                   name = enrich.name,
                                   random.seed = 1)


FeaturePlot(Astrocyte_subset, features = "Reactivity_score1", reduction = "umap.cca", pt.size = 1.5, max.cutoff="q99",min.cutoff="q1", order = T) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + theme(axis.title.x = element_blank(),
                                                                                                                                    axis.text.x = element_blank(),
                                                                                                                                    axis.ticks.x = element_blank(),
                                                                                                                                    axis.title.y = element_blank(),
                                                                                                                                    axis.text.y = element_blank(),
                                                                                                                                    axis.ticks.y = element_blank(),
                                                                                                                                    axis.line = element_blank())


df_EMT <- data.frame(gene_symbol = Sofroniew_2022_EMT_genes)

conversion_result_EMT <- getBM(
  attributes = c("hgnc_symbol", "ensembl_gene_id"),
  filters    = "hgnc_symbol",
  values     = df_EMT$gene_symbol,
  mart       = mart
)

df_converted_EMT <- merge(
  x  = df_EMT,
  y  = conversion_result_EMT,
  by.x = "gene_symbol",
  by.y = "hgnc_symbol",
  all.x = TRUE
)

EMT_genes <- df_converted_EMT[!(is.na(df_converted_EMT$ensembl_gene_id)),"ensembl_gene_id"]

genes.list <- list(EMT_genes)
enrich.name <- "EMT_score"
Astrocyte_subset <- AddModuleScore(Astrocyte_subset,
                                   features = genes.list,
                                   pool = NULL,
                                   n.bin = 5,
                                   seed = 1,
                                   ctrl = length(genes.list),
                                   k = FALSE,
                                   name = enrich.name,
                                   random.seed = 1)


FeaturePlot(Astrocyte_subset, features = "EMT_score1", reduction = "umap.cca", pt.size = 1.5, max.cutoff="q99",min.cutoff="q1", order = T) +
  theme(plot.margin = unit(c(1, 1, 1, 1), "cm"))+ scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "RdBu")))  + theme(axis.title.x = element_blank(),
                                                                                                                                    axis.text.x = element_blank(),
                                                                                                                                    axis.ticks.x = element_blank(),
                                                                                                                                    axis.title.y = element_blank(),
                                                                                                                                    axis.text.y = element_blank(),
                                                                                                                                    axis.ticks.y = element_blank(),
                                                                                                                                    axis.line = element_blank())



Astrocyte_subset <- readRDS(file = "Astrocyte_subset_reannotated.rds")

## Integration with Visium ----
library(SeuratData)
brain <- LoadData("stxBrain", type = "anterior1")
brain <- SCTransform(brain, assay = "Spatial", verbose = FALSE)
brain <- RunPCA(brain, assay = "SCT", verbose = FALSE)
brain <- FindNeighbors(brain, reduction = "pca", dims = 1:30)
brain <- FindClusters(brain, verbose = FALSE)
brain <- RunUMAP(brain, reduction = "pca", dims = 1:30)

brain@images$anterior1@scale.factors$spot <- 89.4725
SpatialFeaturePlot(brain, features = "nCount_Spatial", interactive = FALSE, pt.size.factor = 2) + theme(legend.position = "right")


DefaultAssay(brain) <- "SCT"

Astrocyte_subset <- readRDS(file = "Astrocyte_subset.rds")

DefaultAssay(Astrocyte_subset) <- "RNA"
Astrocyte_subset <- JoinLayers(Astrocyte_subset)
Astrocyte_subset <- SCTransform(Astrocyte_subset)
Astrocyte_subset <- FindVariableFeatures(Astrocyte_subset, selection.method = "vst", nfeatures = 2000)
DefaultAssay(Astrocyte_subset) <- "SCT"

gene_names <- rownames(Astrocyte_subset)
gene_names_updated <- tolower(gene_names)
gene_names_updated <- sapply(gene_names_updated, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
gene_names_updated <- as.character(gene_names_updated)
Astrocyte_subset@assays$SCT@counts@Dimnames[[1]] <- gene_names_updated
Astrocyte_subset@assays$SCT@data@Dimnames[[1]] <- gene_names_updated

scale_data <- tolower(rownames(Astrocyte_subset@assays$SCT@scale.data))
scale_data <- sapply(scale_data, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
rownames(Astrocyte_subset@assays$SCT@scale.data) <- scale_data
rownames(Astrocyte_subset@assays$SCT@meta.features) <- gene_names_updated

var_features <- Astrocyte_subset@assays$SCT@var.features
var_features <- tolower(var_features)
var_features <- sapply(var_features, function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
})
Astrocyte_subset@assays$SCT@var.features <- var_features

saveRDS(Astrocyte_subset, file = "Astrocyte_subset_GenesToLower4Visium.rds")

head(rownames(Astrocyte_subset))

#integrate the datasets
anchors <- FindTransferAnchors(reference = Astrocyte_subset, query = brain, normalization.method = "SCT", 
                               recompute.residuals = F , npcs = 5, dims = 1:5)

predictions.assay <- TransferData(anchorset = anchors, refdata = Astrocyte_subset$integrated_clusters_k0.3, prediction.assay = TRUE,
                                  weight.reduction = brain[["pca"]], dims = 1:5)


brain[["predictions"]] <- predictions.assay

DefaultAssay(brain) <- "predictions"

SpatialFeaturePlot(brain, features = "0", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "1", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "2", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "3", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "4", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "5", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800
SpatialFeaturePlot(brain, features = "6", pt.size.factor = 3, ncol = 1, crop = TRUE) #save 800x800

