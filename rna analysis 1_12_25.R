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

#SEURAT CLUSTERING PART
D03 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D03_rna")
D07 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D07_rna")
D14 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D14_rna")

D03 <- CreateSeuratObject(counts = D03, project = "D03")
D07 <- CreateSeuratObject(counts = D07, project = "D07")
D14 <- CreateSeuratObject(counts = D14, project = "D14")

D03 <- NormalizeData(D03)
D07 <- NormalizeData(D07)
D14 <- NormalizeData(D14)

GSE <- merge(x = D03, y = c(D07, D14))

GSE[["percent.mt"]] <- PercentageFeatureSet(GSE, pattern = "^MT-")
FeatureScatter(GSE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
GSE <- subset(GSE, subset = nCount_RNA < 15000 & nFeature_RNA > 300 & nFeature_RNA < 3700 & percent.mt < 5)

GSE <- FindVariableFeatures(GSE, nfeatures = 2000)
GSE.varplot <- VariableFeaturePlot(GSE)
top10 <- head(VariableFeatures(GSE), n = 10)
LabelPoints(GSE.varplot, points = top10, repel = TRUE, xnudge = 0, ynudge = 0)

all.genes <- rownames(GSE)
GSE <- ScaleData(GSE, features = all.genes)

GSE <- RunPCA(GSE)
ElbowPlot(GSE, reduction = "pca", ndims = 30)

GSE <- FindNeighbors(GSE, dims = 1:24)
GSE <- FindClusters(GSE, resolution = 1)

GSE <- RunUMAP(GSE, dims = 1:24, reduction = "pca")
DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "seurat_clusters", split.by = "orig.ident")
DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "seurat_clusters")

DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "orig.ident")

GSE <- JoinLayers(GSE)
GSE.markers <- FindAllMarkers(GSE, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)
GSE.markers_power <- group_by(GSE.markers, cluster) %>% arrange(power) %>% slice(which.max(power))
view(GSE.markers_power)

GSE.markers_test <- FindAllMarkers(GSE, test.use = "roc")
GSE.markers_test_power <- group_by(GSE.markers_test, cluster) %>% arrange(power) %>% slice(which.max(power))
levels(GSE.markers_test_power$cluster) #bruh clusters 15-18 don't have markers??? actually gonna crash out
view(GSE.markers_test_power)
#new.cluster.ids <- c(filter(GSE.markers_power, cluster == 0)$gene,
#                     filter(GSE.markers_power, cluster == 1)$gene,
#                     filter(GSE.markers_power, cluster == 2)$gene,
#                    filter(GSE.markers_power, cluster == 3)$gene,
#                    filter(GSE.markers_power, cluster == 4)$gene,
#                    filter(GSE.markers_power, cluster == 5)$gene,
#                    filter(GSE.markers_power, cluster == 6)$gene,
#                    filter(GSE.markers_power, cluster == 7)$gene,
#                    filter(GSE.markers_power, cluster == 8)$gene,
#                    filter(GSE.markers_power, cluster == 9)$gene,
#                    filter(GSE.markers_power, cluster == 10)$gene,
#                    filter(GSE.markers_power, cluster == 11)$gene,
#                    filter(GSE.markers_power, cluster == 12)$gene,
#                    filter(GSE.markers_power, cluster == 13)$gene,
#                    filter(GSE.markers_power, cluster == 14)$gene,
#                    "15",
#                    "16",
#                    "17",
#                    "18")
names(new.cluster.ids) <- levels(GSE) #levels(GSE) uses markers from GSE seurat object
#that were discovered using ROC test. highest power marker is displayed by default
GSE <- RenameIdents(GSE,new.cluster.ids)
DimPlot(GSE, reduction = "umap", label = TRUE, repel = FALSE)
DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, 
        group.by = "seurat_clusters", split.by = "orig.ident")
FeaturePlot(object = GSE, features = "AW112010")

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
#reference_dataset <- LogNormalize(reference_dataset)
singleR_results <- SingleR(test = GSE@assays$RNA$counts, ref = reference_dataset,
                           labels = reference_dataset$label.fine)
view(singleR_results)
table(singleR_results$labels)
rownames(singleR_results) %in% colnames(GSE)  #double checking that everything is aligned
GSE[["cell_type"]] <- singleR_results$pruned.labels
view(GSE)

DimPlot(GSE, reduction = "umap", label = TRUE, group.by = "cell_type")





#Granulocytes
GSE_g <- subset(GSE, cell_type == "Granulocytes")
DimPlot(GSE_g, reduction = "umap", split.by = "orig.ident")
DimPlot(GSE_g, reduction = "umap", group.by = "orig.ident")

GSE_g.markers <- FindAllMarkers(GSE_g, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)
GSE_g.markers_power <- group_by(GSE_g.markers, cluster) %>% arrange(power) %>% slice(which.max(power))
view(GSE_g.markers_power)

GSE_g_power.markers <- c()

GSE_g.markers



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
DimPlot(GSE_tcell, reduction = "umap", group.by = "orig.ident", label = TRUE)

GSE_tcell <- JoinLayers(GSE_tcell)
cds_tcell <- SeuratWrappers::as.cell_data_set(GSE_tcell)
cds_tcell <- cluster_cells(cds_tcell)

plot_cells(cds_tcell, show_trajectory_graph = TRUE, color_cells_by = "partition")
cds_tcell <- learn_graph(cds_tcell, use_partition = FALSE)
cds_tcell <- order_cells(cds_tcell)

plot_cells(cds_tcell, color_cells_by = "pseudotime", label_branch_points = FALSE, label_leaves = FALSE)


