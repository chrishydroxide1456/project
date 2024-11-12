library(Seurat)
library(tidyverse)
remotes::install_github("satijalab/seurat-wrappers")
library(SeuratWrappers)
devtools::install_github("satijalab/seurat-data")

install.packages("BiocManager")
BiocManager::install("scRepertoire")
library(ggplot2)
library(scRepertoire)

install.packages("devtools", type = "win.binary")
library(devtools)
devtools::install_github("immunogenomics/presto", dependencies = FALSE) #this doesn't work
library(presto)

library(BiocManager)
BiocManager::install("rhdf5")
# bruh: subject1_blood <- h5read("..Downloads/GSE241739_Subject_5_Blood_molecule_info.h5")


#SEURAT CLUSTERING PART
D03 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D03_rna")
D07 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D07_rna")
D14 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D14_rna")

D03 <- CreateSeuratObject(counts = D03, project = "D03")
D07 <- CreateSeuratObject(counts = D07, project = "D07")
D14 <- CreateSeuratObject(counts = D14, project = "D14")

GSE <- merge(x = D03, y = c(D07, D14))

GSE[["percent.mt"]] <- PercentageFeatureSet(GSE, pattern = "^MT-")
FeatureScatter(GSE, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
GSE <- subset(GSE, subset = nFeature_RNA < 2400 & nFeature_RNA > 300 & percent.mt < 5)

GSE <- NormalizeData(GSE)

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

GSE <- RunUMAP(GSE, dims = 1:24)
DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "seurat_clusters", split.by = "orig.ident")
DimPlot(GSE, reduction = "umap", label = TRUE, repel = TRUE, group.by = "seurat_clusters")


GSE <- JoinLayers(GSE)
GSE.markers <- FindAllMarkers(GSE, test.use = "roc") %>% group_by(cluster) %>% filter(avg_log2FC > 1)
GSE.markers_power <- group_by(GSE.markers, cluster) %>% arrange(power) %>% slice(which.max(power))
view(GSE.markers_power)

GSE.markers_test <- FindAllMarkers(GSE, test.use = "roc")
GSE.markers_test_power <- group_by(GSE.markers_test, cluster) %>% arrange(power) %>% slice(which.max(power))
levels(GSE.markers_test_power$cluster) #clusters 15-18 don't have markers???
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

library(data.table)

view(GSE)
GSE_data <- GetAssayData(object = GSE, assay = "RNA", layer = "scale.data")
fwrite(GSE_data, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_data.csv")
write.csv(GSE_data, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_data.csv")
view(GSE.markers)
GSE_markers <- GetAssayData(object = GSE, assay = "RNA", layer = "scale.data")
fwrite(GSE.markers_power, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_markers.csv")
write.csv(GSE, "C:/Users/Chris Oh/Desktop/seuratdata/GSE242039_RAW/Results/GSE_markers")




#LOOKING AT T CELLS
library(rhdf5)
BiocManager::install("SingleR")
library(SingleR)
BiocManager::install("celldex")
library(celldex)

D03 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D03_rna")
D07 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D07_rna")
D14 <- Read10X(data.dir = "../Desktop/GSE242039 All/GSE242039_RNA/rna_data/D14_rna")

D03 <- CreateSeuratObject(counts = D03, project = "D03")
D07 <- CreateSeuratObject(counts = D07, project = "D07")
D14 <- CreateSeuratObject(counts = D14, project = "D14")

GSE <- merge(x = D03, y = c(D07, D14))
GSE <- NormalizeData(GSE)
GSE <- JoinLayers(GSE) #temporary solution
#GSE <- IntegrateLayers(GSE, method = CCAIntegration) - having problems w/ this

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

GSE_tcell <- subset(GSE, cell_type == "T cells")
FeatureScatter(GSE_tcell, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + geom_smooth(method = "lm")
GSE_tcell <- subset(GSE_tcell, subset = nCount_RNA < 15000 & nFeature_RNA > 1000) #correlation = 0.94
GSE_tcell <- NormalizeData(GSE_tcell)
view(GSE_tcell)

colnames(GSE_tcell)





