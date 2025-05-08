library(Seurat)
library(tidyverse)
library(ggplot2)
library(scRepertoire)
library(devtools)
library(presto)
library(BiocManager)
library(data.table)
library(rhdf5)
library(SingleR)
library(celldex)
library(monocle3)
library(monocle)
library(remotes)
library(R.utils)
library(SeuratWrappers)
library(immunarch)
library(scMayoMap)

#install celldex and get MouseRNAseqData dataset to report later
#look for tcell receptor and general receptor seq tool
#try to improve pseudotime analysis
#adjust QC parameters

D03 <- Read10X(data.dir = "../Desktop/sratoolkit/bin/cellranger_results/day_3/")
D07 <- Read10X(data.dir = "../Desktop/sratoolkit/bin/cellranger_results/day_7/")
D14 <- Read10X(data.dir = "../Desktop/sratoolkit/bin/cellranger_results/day_14/")

D03 <- CreateSeuratObject(counts = D03, project = "D03")
D07 <- CreateSeuratObject(counts = D07, project = "D07")
D14 <- CreateSeuratObject(counts = D14, project = "D14")

GSE <- merge(x = D03, y = c(D07, D14))
GSE <- NormalizeData(GSE)

GSE[["percent.mt"]] <- PercentageFeatureSet(GSE, pattern = "^MT-")
FeatureScatter(GSE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
GSE <- subset(GSE, subset = nCount_RNA < 15000 & nFeature_RNA > 300 & nFeature_RNA < 3000 & percent.mt < 5)
#GSE <- subset(GSE, subset = percent.mt < 5)

GSE <- FindVariableFeatures(GSE, nfeatures = 4000)
GSE.varplot <- VariableFeaturePlot(GSE)
top10 <- head(VariableFeatures(GSE), n = 10)
LabelPoints(GSE.varplot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

GSE <- ScaleData(GSE)

GSE <- RunPCA(GSE)
ElbowPlot(GSE, reduction = "pca", ndims = 40)

GSE <- IntegrateLayers(object = GSE, method = CCAIntegration,
                       orig.reduction = "pca", new.reduction = "cca",
                       verbose = FALSE)
GSE <- JoinLayers(GSE)

GSE <- FindNeighbors(GSE, dims = 1:34, reduction = "cca")
GSE <- FindClusters(GSE, resolution = 0.6, cluster.name = "cca_clusters",
                    algorithm = 4)

GSE <- RunUMAP(GSE, dims = 1:34, reduction.name = "cca")
DimPlot(GSE, reduction = "cca", group.by = "orig.ident")
DimPlot(GSE, reduction = "cca", group.by = "cca_clusters")

view(GSE.markers_test)
GSE.markers_test <- FindAllMarkers(GSE)
GSE.markers_test <- FindAllMarkers(GSE, test.use = "roc")
view(GSE.markers_test)

GSE.markers <- FindAllMarkers(GSE, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)
GSE.markers_power <- group_by(GSE.markers, cluster) %>% arrange(power) %>% slice(which.max(power))
view(GSE.markers_power)



view(GSE)
GSE_data <- GetAssayData(object = GSE, assay = "RNA", layer = "scale.data")
fwrite(GSE_data, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_data.csv")
write.csv(GSE_data, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_data.csv")
view(GSE.markers)
GSE_markers <- GetAssayData(object = GSE, assay = "RNA", layer = "scale.data")
fwrite(GSE.markers_power, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_markers.csv")
write.csv(GSE, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_markers")




#SingleR
listReferences()
reference_dataset <- MouseRNAseqData()
reference_dataset <- NormalizeData(reference_dataset)
MouseRNAseqData()
BiocManager::install("celldex")

singleR_results <- SingleR(test = GSE@assays$RNA$counts, ref = reference_dataset,
                           labels = reference_dataset$label.fine)
view(singleR_results)
table(singleR_results$labels)
rownames(singleR_results) %in% colnames(GSE)  #double checking that everything is aligned


mouse_immune <- ImmGenData()
mouse_immune <- NormalizeData(reference_dataset)
mouse_immune <- LogNormalize(reference_dataset)

immgen <- SingleR(test = GSE@assays$RNA$counts, ref = mouse_immune,
                         labels = mouse_immune$label.main)
immgen_fine <- SingleR(test = GSE@assays$RNA$counts, ref = mouse_immune,
                           labels = mouse_immune$label.fine)
table(immgen$labels)
table(immgen_fine$pruned.labels)
rownames(immgen) %in% colnames(GSE)  #double checking that everything is aligned

GSE[["cell_type"]] = immgen$labels
GSE[["cell_type refined"]] <- immgen_fine$pruned.labels
view(GSE)

DimPlot(GSE, reduction = "cca", label = TRUE, group.by = "cell_type", repel = TRUE)




GSE <- JoinLayers(GSE)
cds <- SeuratWrappers::as.cell_data_set(GSE)
cds <- cluster_cells(cds)

plot_cells(cds, show_trajectory_graph = FALSE, color_cells_by = "partition")
cds <- learn_graph(cds, use_partition = FALSE)
cds <- order_cells(cds)

plot_cells(cds, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)


DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident")


#T cells
GSE_tcell <- subset(GSE, cell_type == "T cells")
DimPlot(GSE_tcell, reduction = "cca", split.by = "seurat_clusters")
DimPlot(GSE_tcell, reduction = "cca", group.by = "seurat_clusters")
DimPlot(GSE_tcell, reduction = "cca", group.by = "orig.ident")


ncol(subset(GSE_tcell, orig.ident == "D03"))
ncol(subset(GSE_tcell, orig.ident == "D07"))
ncol(subset(GSE_tcell, orig.ident == "D14"))

asdfasdf
GSE_tcell <- FindNeighbors(GSE_tcell, dims = 1:2, reduction = "cca")
GSE_tcell <- FindClusters(GSE_tcell, resolution = 0.6, cluster.name = "cca_clusters",
                      algorithm = 4)
DimPlot(GSE_tcell, reduction = "cca")
GSE_tcell.markers <- FindAllMarkers(GSE_tcell, test.use = "roc")
view(GSE_tcell.markers)

FeaturePlot(GSE_tcell, reduction = "cca", features = c("Trbc1", "Trbc2"), blend = TRUE)
?FeaturePlot
GSE_tcell <- JoinLayers(GSE_tcell)
cds_tcell <- SeuratWrappers::as.cell_data_set(GSE_tcell)
cds_tcell <- cluster_cells(cds_tcell)

plot_cells(cds_tcell, show_trajectory_graph = TRUE, color_cells_by = "partition")
cds_tcell <- learn_graph(cds_tcell, use_partition = FALSE)
cds_tcell <- order_cells(cds_tcell)

plot_cells(cds_tcell, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)


FeaturePlot(GSE_tcell, features = c("Stmn1"), split.by = "orig.ident", reduction = "cca")
FeaturePlot(GSE_tcell, features = c("Rpl23"), split.by = "orig.ident", reduction = "cca")
FeaturePlot(GSE_tcell, features = c("mt-Atp6"), split.by = "orig.ident", reduction = "cca")
FeaturePlot(GSE_tcell, features = c("Fau"), split.by = "orig.ident", reduction = "cca")

GSE_tcell3 <- subset(GSE_tcell, orig.ident == "D03")
GSE_tcell7 <- subset(GSE_tcell, orig.ident == "D07")
GSE_tcell14 <- subset(GSE_tcell, orig.ident == "D14")

FeaturePlot(GSE_tcell, features = c("Trbc2"), reduction = "cca", split.by = "orig.ident")

d_3Stmn1 <- FetchData(GSE_tcell3, vars = "Trbc2")
d_3Stmn1 <- subset(d_3Stmn1, Trbc2 > 0)
d_7Stmn1 <- FetchData(GSE_tcell7, vars = "Trbc2")
d_7Stmn1 <- subset(d_7Stmn1, Trbc2 > 0)
d_14Stmn1 <- FetchData(GSE_tcell14, vars = "Trbc2")
d_14Stmn1 <- subset(d_14Stmn1, Trbc2 > 0)
nrow(d_3Stmn1)
nrow(d_7Stmn1)
nrow(d_14Stmn1)



#Granulocytes
GSE_N <- subset(GSE, cell_type == "Neutrophils" | cell_type == "Basophils")
DimPlot(GSE_N, reduction = "cca", split.by = "orig.ident")
DimPlot(GSE_N, reduction = "cca", group.by = "orig.ident")
view(GSE_N)
view(GSE_g.markers)

GSE_Gs <- subset(GSE_N, Cd177 == 0)
view(GSE_Gs)
DimPlot(GSE_Gs, reduction = "cca", split.by = "orig.ident")


Granulocytes <- colnames(GSE_N)
DimPlot(GSE_N, reduction = "cca", split.by = "orig.ident")

ncol(subset(GSE_Gs, orig.ident == "D03"))
ncol(subset(GSE_Gs, orig.ident == "D07"))
ncol(subset(GSE_Gs, orig.ident == "D14"))

GSE_N3 <- subset(GSE_Gs, orig.ident == "D03")
GSE_N7 <- subset(GSE_Gs, orig.ident == "D07")
GSE_N14 <- subset(GSE_Gs, orig.ident == "D14")

FeaturePlot(GSE_Gs, features = "Arg2", split.by = "orig.ident", reduction = "cca")

d3_S100a9 <- FetchData(GSE_N3, vars = "Arg2")
d3_S100a9 <- subset(d3_S100a9, Arg2 > 0)
d7_S100a9 <- FetchData(GSE_N7, vars = "Arg2")
d7_S100a9 <- subset(d7_S100a9, Arg2 > 0)
d14_S100a9 <- FetchData(GSE_N14, vars = "Arg2")
d14_S100a9 <- subset(d14_S100a9, Arg2 > 0)
nrow(d3_S100a9)
nrow(d7_S100a9)
nrow(d14_S100a9)







GSE_g.markers[["day"]][which(GSE.markers$cluster == "16")] = "D03"

view(GSE_g.markers)
view(GSE_g)
DimPlot(GSE_g, reduction = "cca", group.by = "orig.ident")
GSE_g.markers <- FindAllMarkers(GSE_g, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)

GSE_g <- JoinLayers(GSE_g)
cds_g <- SeuratWrappers::as.cell_data_set(GSE_g)
cds_g <- cluster_cells(cds_g)

#plot_cells(cds_g, show_trajectory_graph = TRUE, color_cells_by = "partition")
cds_g <- learn_graph(cds_g, use_partition = FALSE)
cds_g <- order_cells(cds_g)

plot_cells(cds_g, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)

#Monocytes
view(GSE_m)
GSE_m <- subset(GSE, cell_type == "Monocytes")
ncol(subset(GSE_m, orig.ident == "D03"))
ncol(subset(GSE_m, orig.ident == "D07"))
ncol(subset(GSE_m, orig.ident == "D14"))
GSE_m <- subset(GSE, cell_type == "Monocytes")
GSE_m <- subset(GSE, cell_type == "Monocytes")
DimPlot(GSE_m, reduction = "cca", split.by = "orig.ident")
DimPlot(GSE_m, reduction = "cca", group.by = "orig.ident")

Monocytes <- colnames(GSE_m)
DimPlot(GSE_m, reduction = "cca", split.by = "orig.ident")

GSE_m3 <- subset(GSE_m, orig.ident == "D03")
GSE_m7 <- subset(GSE_m, orig.ident == "D07")
GSE_m14 <- subset(GSE_m, orig.ident == "D14")
FeaturePlot(GSE_m, features = "S100a9", split.by = "orig.ident", reduction = "cca")

d3_S100a9 <- FetchData(GSE_m3, vars = "S100a9")
d3_S100a9 <- subset(d3_S100a9, S100a9 > 0)
d7_S100a9 <- FetchData(GSE_m7, vars = "S100a9")
d7_S100a9 <- subset(d7_S100a9, S100a9 > 0)
d14_S100a9 <- FetchData(GSE_m14, vars = "S100a9")
d14_S100a9 <- subset(d14_S100a9, S100a9 > 0)
nrow(d3_S100a9)
nrow(d7_S100a9)
nrow(d14_S100a9)
view(GSE_m)

FeaturePlot(GSE_m, features = "Arg2", split.by = "orig.ident", reduction = "cca")
d3_Arg2 <- FetchData(GSE_m3, vars = "Arg2")
d3_Arg2 <- subset(d3_Arg2, Arg2 > 0)
d7_Arg2 <- FetchData(GSE_m7, vars = "Arg2")
d7_Arg2 <- subset(d7_Arg2, Arg2 > 0)
d14_Arg2 <- FetchData(GSE_m14, vars = "Arg2")
d14_Arg2 <- subset(d14_Arg2, Arg2 > 0)
nrow(d3_Arg2)
nrow(d7_Arg2)
nrow(d14_Arg2)

#Macrophages
view(GSE_M)
GSE_M <- subset(GSE, cell_type == "Macrophages")
ncol(subset(GSE_m, orig.ident == "D03"))
ncol(subset(GSE_m, orig.ident == "D07"))
ncol(subset(GSE_m, orig.ident == "D14"))



DimPlot(GSE_M, reduction = "cca", split.by = "orig.ident")
DimPlot(GSE_M, reduction = "cca", group.by = "orig.ident")

Macrophages <- colnames(GSE_M)
DimPlot(GSE_M, reduction = "cca", split.by = "orig.ident")

?FindNeighbors()
GSE_M <- FindNeighbors(GSE_M, dims = 1:10, reduction = "cca")
GSE_M <- FindClusters(GSE_M, resolution = 0.6, cluster.name = "cca_clusters",
                    algorithm = 4)
view(GSE_M)

GSE_M <- RunUMAP(GSE_M, dims = 1:34, reduction.name = "cca")

GSE_M3 <- subset(GSE_M, orig.ident == "D03")
ncol(GSE_M3)
GSE_M7 <- subset(GSE_M, orig.ident == "D07")
ncol(GSE_M7)
GSE_M14 <- subset(GSE_M, orig.ident == "D14")
ncol(GSE_M14)

FeaturePlot(GSE_M, features = "Tnfaip2", split.by = "orig.ident", reduction = "cca")
FeaturePlot(GSE_M, features = "Il10", split.by = "orig.ident", reduction = "cca")

d3_S100a9 <- FetchData(GSE_M3, vars = "Il10")
d3_S100a9 <- subset(d3_S100a9, Il10 > 0)
d7_S100a9 <- FetchData(GSE_M7, vars = "Il10")
d7_S100a9 <- subset(d7_S100a9, Il10 > 0)
d14_S100a9 <- FetchData(GSE_M14, vars = "Il10")
d14_S100a9 <- subset(d14_S100a9, Il10 > 0)
nrow(d3_S100a9)
nrow(d7_S100a9)
nrow(d14_S100a9)






ncol(subset(GSE_g, orig.ident == "D03"))
ncol(subset(GSE_g, orig.ident == "D07"))
ncol(subset(GSE_g, orig.ident == "D14"))

view(GSE_g.markers)

DimPlot(GSE_g, reduction = "cca", group.by = "orig.ident")
GSE_g.markers <- FindAllMarkers(GSE_g, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)

GSE_g <- JoinLayers(GSE_g)
cds_g <- SeuratWrappers::as.cell_data_set(GSE_g)
cds_g <- cluster_cells(cds_g)

plot_cells(cds_g, show_trajectory_graph = TRUE, color_cells_by = "partition")
cds_g <- learn_graph(cds_tcell, use_partition = FALSE)
cds_g <- order_cells(cds_tcell)

plot_cells(cds_g, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)


colnames(singleR_results)



GSE_d <- subset(GSE, cell_type == "Dendritic cells")
ncol(subset(GSE_d, orig.ident == "D03"))
ncol(subset(GSE_d, orig.ident == "D07"))
ncol(subset(GSE_d, orig.ident == "D14"))

GSE_M <- subset(GSE, cell_type == "Macrophages")
ncol(subset(GSE_M, orig.ident == "D03"))
ncol(subset(GSE_M, orig.ident == "D07"))
ncol(subset(GSE_M, orig.ident == "D14"))
DimPlot(GSE_M, split.by = "orig.ident", reduction = "cca")

GSE_MA <- subset(GSE, cell_type == "Macrophages activated")
ncol(subset(GSE_MA, orig.ident == "D03"))
ncol(subset(GSE_MA, orig.ident == "D07"))
ncol(subset(GSE_MA, orig.ident == "D14"))




