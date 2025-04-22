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
GSE <- FindClusters(GSE, resolution = 1.5, cluster.name = "cca_clusters",
                    algorithm = 4)

GSE <- RunUMAP(GSE, dims = 1:34, reduction.name = "cca")
DimPlot(GSE, reduction = "cca", group.by = "orig.ident")
DimPlot(GSE, reduction = "cca", group.by = "cca_clusters")

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
MouseRNAseqData()
BiocManager::install("celldex")

#reference_dataset <- LogNormalize(reference_dataset)
singleR_results <- SingleR(test = GSE@assays$RNA$counts, ref = reference_dataset,
                           labels = reference_dataset$label.fine)
view(singleR_results)
table(singleR_results$labels)
rownames(singleR_results) %in% colnames(GSE)  #double checking that everything is aligned
GSE[["cell_type"]] <- singleR_results$pruned.labels
view(GSE)

DimPlot(GSE, reduction = "umap", label = TRUE, group.by = "cell_type")




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
DimPlot(GSE_tcell, reduction = "umap", split.by = "orig.ident")
DimPlot(GSE_tcell, reduction = "umap", group.by = "orig.ident")
DimPlot(GSE_tcell, reduction = "umap", group.by = "orig.ident")

T_cells <- colnames(GSE_tcell)
DimPlot(GSE, reduction = "umap", split.by = "orig.ident", cells.highlight = T_cells)

ncol(subset(GSE_tcell, orig.ident == "D03"))
ncol(subset(GSE_tcell, orig.ident == "D07"))
ncol(subset(GSE_tcell, orig.ident == "D14"))

JoinLayers(GSE_tcell)
GSE_tcell.markers <- FindAllMarkers(GSE_tcell, test.use = "roc")
view(GSE_tcell.markers)


DimPlot(GSE_tcell, reduction = "umap", group.by = "orig.ident")
GSE_tcell.markers <- FindAllMarkers(GSE_tcell, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)

GSE_tcell11 <- subset(GSE_tcell.markers, GSE_tcell.markers$cluster == "11")
GSE_tcell15 <- subset(GSE_tcell.markers, GSE_tcell.markers$cluster == "15")
GSE_tcell18 <- subset(GSE_tcell.markers, GSE_tcell.markers$cluster == "18")
GSE_tcell_power.markers <- c(GSE_tcell11$gene[which.max(GSE_tcell11$power)],
                             GSE_tcell15$gene[which.max(GSE_tcell15$power)],
                             GSE_tcell18$gene[which.max(GSE_tcell18$power)])
GSE_tcell_power.markers[1]
GSE_tcell_power.markers[2]
GSE_tcell_power.markers[3]

new.cluster.ids <- c(GSE_tcell_power.markers[1],
                     GSE_tcell_power.markers[2],
                     GSE_tcell_power.markers[3])
names(new.cluster.ids) <- levels(GSE_tcell) #levels(GSE) uses markers from GSE seurat object
#that were discovered using ROC test. highest power marker is displayed by default
GSE_tcell <- RenameIdents(GSE_tcell,new.cluster.ids)
DimPlot(GSE_tcell, reduction = "umap", group.by = "orig.ident", label = FALSE)

GSE_tcell <- JoinLayers(GSE_tcell)
cds_tcell <- SeuratWrappers::as.cell_data_set(GSE_tcell)
cds_tcell <- cluster_cells(cds_tcell)

plot_cells(cds_tcell, show_trajectory_graph = TRUE, color_cells_by = "partition")
cds_tcell <- learn_graph(cds_tcell, use_partition = FALSE)
cds_tcell <- order_cells(cds_tcell)

plot_cells(cds_tcell, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)

FeaturePlot(GSE_tcell, features = c("Tpt1"))
FeaturePlot(GSE_tcell, features = c("Stmn1"), split.by = "orig.ident")
FeaturePlot(GSE_tcell, features = c("Rpl23"))
FeaturePlot(GSE_tcell, features = c("mt-Atp6"))
FeaturePlot(GSE_tcell, features = c("Fau"))
view(GSE_tcell)
view(GSE_tcell.markers)



#Granulocytes
GSE_g <- subset(GSE, cell_type == "Granulocytes")
DimPlot(GSE_g, reduction = "umap", split.by = "orig.ident")
DimPlot(GSE_g, reduction = "umap", group.by = "orig.ident")

view(GSE_g.markers)

Granulocytes <- colnames(GSE_g)
DimPlot(GSE, reduction = "umap", split.by = "orig.ident", cells.highlight = Granulocytes)

ncol(subset(GSE_g, orig.ident == "D03"))
ncol(subset(GSE_g, orig.ident == "D07"))
ncol(subset(GSE_g, orig.ident == "D14"))

GSE_g.markers[["day"]][which(GSE.markers$cluster == "16")] = "D03"

view(GSE_g.markers)
view(GSE_g)
DimPlot(GSE_g, reduction = "umap", group.by = "orig.ident")
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
ncol(subset(GSE_m, orig.ident == "D03"))
ncol(subset(GSE_m, orig.ident == "D07"))
ncol(subset(GSE_m, orig.ident == "D14"))
GSE_m <- subset(GSE, cell_type == "Monocytes")
GSE_m <- subset(GSE, cell_type == "Monocytes")
DimPlot(GSE_m, reduction = "umap", split.by = "orig.ident")
DimPlot(GSE_m, reduction = "umap", group.by = "orig.ident")

Monocytes <- colnames(GSE_m)
DimPlot(GSE_m, reduction = "umap", group.by = "orig.ident")

GSE_m3 <- subset(GSE_m, orig.ident == "D03")
GSE_m7 <- subset(GSE_m, orig.ident == "D07")
GSE_m14 <- subset(GSE_m, orig.ident == "D14")
FeaturePlot(GSE_m, features = "S100a9", split.by = "orig.ident")

d3_S100a9 <- FetchData(GSE_m3, vars = "S100a9")
d3_S100a9 <- subset(d3_S100a9, "S100a9" > 0)
nrow(d3_S100a9)
ncol(GSE_m3)

d7_S100a9 <- FetchData(GSE_m7, vars = "S100a9")
d7_S100a9 <- subset(d7_S100a9, "S100a9" > 0)
nrow(d7_S100a9)
ncol(GSE_m7)

d14_S100a9 <- FetchData(GSE_m14, vars = "S100a9")
d14_S100a9 <- subset(d14_S100a9, "S100a9" > 0)
nrow(d14_S100a9)
ncol(GSE_m14)


FeaturePlot(GSE_m, features = "Arg2", split.by = "orig.ident")





ncol(subset(GSE_g, orig.ident == "D03"))
ncol(subset(GSE_g, orig.ident == "D07"))
ncol(subset(GSE_g, orig.ident == "D14"))

view(GSE_g.markers)

DimPlot(GSE_g, reduction = "umap", group.by = "orig.ident")
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
DimPlot(GSE_m, split.by = "orig.ident")

GSE_MA <- subset(GSE, cell_type == "Macrophages activated")
ncol(subset(GSE_MA, orig.ident == "D03"))
ncol(subset(GSE_MA, orig.ident == "D07"))
ncol(subset(GSE_MA, orig.ident == "D14"))




