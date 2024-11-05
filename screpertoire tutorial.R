library(scRepertoire)
library(Seurat)
library(tidyverse)
library(data.table)


#DATA LOADING/MANIPULATION

#load in annotated T cell contig data like this
S_10X <- read.csv(file = "../Desktop/TCR data/filtered_contig_annotations.csv")

#if barcodes have prefixes, remove prefixes like this
barcodes <- read.csv(file = "whatever data")
barcodes <- stripBarcode(barcodes)


#for the tutorial, we're gonna use example data built into scRepertoire
data("contig_list")
View(contig_list)

#combineTCR: pairs sequences of alpha and beta chain from the same cell
    #samples: sample ID, equivalent of 'project' in Seurat. labels data from each element in the original list
    #removeMulti: remove cells with more than 2 chains?
    #filterMulti: only take 2 most expressed chains from cells with  > 2 chains?
    #removeNA: remove all rows with any 'NA' value?
    #ID: optional additional sample labeling, adds suffix to "samples"
combined_TCR <- combineTCR(input.data = contig_list, samples = c("P17", "P17", "P18",
                                                                 "P18", "P19", "P19",
                                                                 "P20", "P20"),
                           ID = c("B", "L", "B", "L", "B", "L", "B", "L"),
                           removeMulti = FALSE,
                           filterMulti = FALSE,
                           removeNA = FALSE)

View(combined_TCR)
#use "combined_TCR$samplename" to look up data
P17_B <- combined_TCR$P17_B
View(P17_B)

#exportClones: export TCR object
    #format: paired (default value) = export with paired chain seqs
view(exportClones(combined_TCR, format = "paired", write.file = FALSE))
    #write.file: TRUE = save as file, FALSE = show dataframe (i guess to view)
view(exportClones(combined_TCR, write.file = FALSE))

#add variables like this
addvar <- addVariable(combined_TCR, variable.name = "type", variables = c("skibidi", "sigma", "skibidi", "sigma", "skibidi", "sigma", "skibidi", "sigma"))
View(addvar$P17_L)
combined_TCR <- addVariable(combined_TCR, variable.name = "aura", variables = c("skibidi", "sigma", "skibidi", "sigma", "skibidi", "sigma", "skibidi", "sigma"))
view(combined_TCR$P17_B)

#subset contigs with "subsetClones"
only19 <- subsetClones(combined_TCR, name = "sample", variables = c("P19"))
View(only19)
only19B <- subsetClones(only19, name = "ID", variables = c("B"))


#VISUALIZATIONS

#arguments
#cloneCall: how you categorize cells.
    #"gene" = v(d)jc gene seq b4 v(d)j recombination, "nt" = nucleotide seq (TCR seq after v(d)j recombination),
    #"strict" = both (i think), "aa" = amino acid seq
    #note: clonalLength only accepts "nt" and "aa" as arguments for cloneCall
#chain: "both" = alpha and beta, "TRA" = alpha, "TRB" = beta
    #note: vizGenes doesn't accept "both"
#scale: show relative distributions
#exportTable: TRUE = return data frame, FALSE = return graph

#clonalAbundance: produce a # of clones vs abundance graph
    #run with exportTable = TRUE to get a data frame
    #run with exportTable = FALSE to get a graph
cA_aa <- clonalAbundance(combined_TCR, cloneCall = "aa", chain = "both",
                         exportTable = FALSE)
view(cA_aa)
cA_nt <- clonalAbundance(combined_TCR, cloneCall = "nt", chain = "both",
                        exportTable = TRUE)
view(cA_nt)
cA_gene <- clonalAbundance(combined_TCR, cloneCall= "gene", chain = "both",
                           exportTable = TRUE)
view(cA_gene)
#run clonalAbundance without exportTable = FALSE to get a graph
clonalAbundance(combined_TCR, cloneCall = "aa", chain = "both", exportTable = FALSE)

#clonalQuant: clonotype quantification
    #cloneCall = "strict" is the most specific
    #scale = FALSE returns # of unique clonotypes in the same
clonalQuant(combined_TCR, cloneCall = "strict", chain = "both")
clonalQuant(combined_TCR, cloneCall = "strict", chain = "both", exportTable = T)
    #scale = TRUE returns % of clonotypes in the sample that are unique
clonalQuant(combined_TCR, cloneCall = "strict", chain = "both", scale = T)
clonalQuant(combined_TCR, cloneCall = "strict", chain = "both", scale = T, exportTable = T)
    #group.by: factors of this column from the data frame are column names in graph. default value is "samples"
clonalQuant(combined_TCR, cloneCall = "strict", chain = "both", group.by = "aura", scale = T)

#checkContig: check NA values
#WHERE TF IS THIS??? WHY DID THEY RENAME IT???

#clonalLength: graph of # of CDR3 AAs vs length. tells you how many CDR3 seqs
# (AA seq or NT seq, depends on which one you choose) have a certain length
clonalLength(combined_TCR, cloneCall = "aa", chain = "both")
clonalLength(combined_TCR, cloneCall = "nt", chain = "both", scale = T)
clonalLength(combined_TCR, cloneCall = "gene", chain = "both") #this doesn't work

#clonalCompare: compare clonotypes. determine which clonotypes are unique to a
#certain sample or cell type, compare proportions of a clonotype within different
#samples or cell types, determine how similar or different cell types/samples are
  #top.clones: use top (x) abundant clonotypes
top10_P17 <- clonalAbundance(combined_TCR, cloneCall = "strict", chain = "both", exportTable = TRUE)
top10a_P17 <- arrange(top10_P17, Abundance)
top10a_P17 <- top10a_P17[rev(1:nrow(top10a_P17)), ] %>% head(, n = 10)
view(top10a_P17)
#why don't these two return the same thing??? i trust the 1st one tho
clonalCompare(combined_TCR, cloneCall = "strict", chain = "both",
              samples = c("P17_L", "P17_B"), top.clones = 10, graph = "alluvial")
clonalCompare(combined_TCR, cloneCall = "strict", chain = "both", samples = c("P17_L", "P17_B"), graph = "alluvial", clones = top10a_P17$CTstrict)

clonalCompare(combined_TCR, cloneCall = "aa", chain = "both", samples = c("P17_B", "P17_L"), graph = "alluvial", top.clones = 10)
#no shared clonotypes bc different samples
clonalCompare(combined_TCR, cloneCall = "strict", chain = "both", samples = c("P17_B", "P18_B"), graph = "alluvial", top.clones = 10)

#vizGenes: visualize gene usage
    #gene = "V", "D" (only for beta chains), "J", or "C"
    #chain = TRB(V, D, J) or TRA(V, J)
    #plot: barplot or heatmap
    #order: order x axis by variance or gene?
    #scale: TRUE = get %, FALSE = get the #
#***don't put an argument for y axis if you only want to visualize 1 gene segment
#default y axis for barplot is gene usage
vizGenes(combined_TCR, x.axis = "TRBJ", plot = "barplot", order = "variance", scale = TRUE, exportTable = FALSE)
#default y axis for heatmap is sample
vizGenes(combined_TCR, x.axis = "TRBV", plot = "heatmap", order = "variance", scale = TRUE, exportTable = FALSE)

#find which VJ combinations are common in the beta chains of P17_B and P17_L
vizGenes(combined_TCR[[1]], x.axis = "TRBV", y.axis = "TRBJ", plot = "heatmap", order = "gene")
vizGenes(combined_TCR[[2]], x.axis = "TRBV", y.axis = "TRBJ", plot = "heatmap", order = "gene")

#ex. only visualizing V gene segment usage in alpha chains from B samples
vizGenes(combined_TCR[c(1, 3, 5, 7)], x.axis = "TRAV", y.axis = "TRAJ", plot = "heatmap", order = "variance", scale = TRUE)
#ex. only visualizing V gene segment usage in alpha chains from L samples
vizGenes(combined_TCR[c(2, 4, 6, 8)], x.axis = "TRAV", plot = "heatmap", order = "variance", scale = TRUE)

#clonalProportion: visualize clonal proportion. "clonal indices" = abundance
    #clonalSplit: determines the range of abundance that each color covers
clonalProportion(combined_TCR, cloneCall = "strict")
clonalProportion(combined_TCR, cloneCall = "strict", clonalSplit = c(1, 500, 1000, 5000, 10000, 30000, 100000))

#clonalHomeostasis: calculates the space that is occupied by certain size ranges
#of cells
clonalHomeostasis(combined_TCR, cloneCall = "strict",
                  cloneSize = c(Rare = 1e-04, Small = 0.001, Medium = 0.01,
                   Large = 0.1, Hyperexpanded = 1)) #this is the default setting for cloneSize

#clonalOverlap: measure similarity between different samples
    #method: "overlap" = overlap coefficient, "jaccard", "morisita",
    #"raw" = exact #, "cosine" = cosine
    #***"clone" is determined by cloneCall (categorization criteria)
clonalOverlap(combined_TCR, cloneCall = "strict", method = "overlap") #categorize cells by V(D)JC seq and CDR3 nucleotide (after V(D)JC recombination)
clonalOverlap(combined_TCR, cloneCall = "strict", method = "raw") #exact number of overlapping clones
clonalOverlap(combined_TCR, cloneCall = "gene", method = "overlap")
clonalOverlap(combined_TCR, cloneCall = "nt", method = "overlap")
clonalOverlap(combined_TCR, cloneCall = "aa", method = "overlap") #categorize cells by amino acid sequence (building blocks of proteins and TCR)

#clonalSizeDistribution: cluster samples based on clone size distribution
clonalSizeDistribution(combined_TCR, cloneCall = "strict", chain = "both", method = "ward.D2")


#DIVERSITY ANALYSIS
#1.) Shannon: estimate baseline diversity
#2.) inverse Simpson: estimate baseline diversity
#3.) Chao: estimate richness of samples
#4.) ACE: abundance-based coverage estimator. estimate richness of samples
#5.) inverse Pielou: measure clonotype evenness

#clonalDiversity: get diversity index values for multiple systems.
    #bootstrapping: identify T cell clones that change significantly in response to therapy
#???? why doesnt the x axis get split up by sample???? wtf
clonalDiversity(combined_TCR, cloneCall = "gene", group.by = "sample", x.axis = "ID", n.boots = 10)

#clonalScatter: compare abundance of clonotypes from different samples#clonasample()lScatter: compare abundance of clonotypes from different samples
clonalScatter(combined_TCR, cloneCall = "strict", chain = "both", x.axis = "P17_B",
              y.axis = "P17_L", dot.size = "total", graph = "proportion")


#INTEGRATING WITH SC OBJECTS
seuratO <- readRDS("../Desktop/TCR data/screp_FullExample.rds")
view(seuratO)

DimPlot(seuratO, label = TRUE, group.by = "seurat_clusters", reduction = "umap")
DimPlot(seuratO, label = TRUE, group.by = "integrated_snn_res.0.5", reduction = "umap")
DimPlot(seuratO, label = TRUE, group.by = "Type", reduction = "umap")

#get the # of cells in each cluster
table(Idents(seuratO))

#FeaturePlot: visualize expression of a gene in the clusters
FeaturePlot(seuratO, features = c("CD4", "CD8B"), label = TRUE)

#combineExpression: combine seurat scRNA data and scTCR data
combined_TCR$P17_B[["patient"]] <- paste(combined_TCR$P17_B$sample, combined_TCR$P17_B$ID, sep = "_")
combined_TCR$P17_L[["patient"]] <- paste(combined_TCR$P17_L$sample, combined_TCR$P17_L$ID, sep = "_")
combined_TCR$P18_B[["patient"]] <- paste(combined_TCR$P18_B$sample, combined_TCR$P18_B$ID, sep = "_")
combined_TCR$P18_L[["patient"]] <- paste(combined_TCR$P18_L$sample, combined_TCR$P18_L$ID, sep = "_")
combined_TCR$P19_B[["patient"]] <- paste(combined_TCR$P19_B$sample, combined_TCR$P19_B$ID, sep = "_")
combined_TCR$P19_L[["patient"]] <- paste(combined_TCR$P19_L$sample, combined_TCR$P19_L$ID, sep = "_")
combined_TCR$P20_B[["patient"]] <- paste(combined_TCR$P20_B$sample, combined_TCR$P20_B$ID, sep = "_")
combined_TCR$P20_L[["patient"]] <- paste(combined_TCR$P20_L$sample, combined_TCR$P20_L$ID, sep = "_")
seurat_combined <- combineExpression(input.data = combined_TCR, sc.data = seuratO,
                                     cloneCall = "gene", group.by = "Patient", split.by = "Type", proportion = FALSE,
                                     cloneSize = c(Single=1, Small=5, Medium=20, Large=100, Hyperexpanded=500))

#split.by = "Type": split a sample into 2 cell types with
#group.by = "clonoType" (measures clonal proliferation. rare, small, large, etc.): different colors to identify different clonal proliferation
DimPlot(seurat_combined, group.by = "cloneType", split.by = "Type")

#clonalOverlay: visualize clonal expansion (proliferation of a clonotype)
    #
clonalOverlay(seurat_combined, reduction = "umap", freq.cutpoint = 30, bins = 12, facet = "Type" + guides(color = "none"))

#clonalNetwork: track T-cells through proliferation and differentiation
#ex. generate a clonal network map for cluster 2
clonalNetwork(seurat_combined, reduction = "umap", identity = "ident",
              filter.ident = "C2", cloneCall = "aa")

#highlightClones: find clonotypes with the highest frequency of a certain aa/nt/gene sequence
#ex. highlight for a certain amino acid sequence
seurat_combined <- highlightClones(seurat_combined, cloneCall = "aa",
                                   sequence = "CAVNGGSQGNLIF_CSAEREDTDTQYF")

#occupiedscRepertoire(): show proportions of clonotypes among single cell clusters

#StartracDiversity(): visualize clonal expansion, cross-tissue migration, and
                      #state transition of clonotypes.
                      #use group.by to determine x axis

#clonotypeBias(): 

#clonalBias(): show how individual clones are skewed to a certain cluster
    #x axis = clone size, y axis = clone bias

#alluvialClones(): generate an alluvial graph to show correlations between multiple
                  #y axes (ex. y.axes = c("Patient", "ident", "Type"))
                  #color = "": select individual genes from V(D)JC segment



#VISUALIZE INTERCONNECTIONS OF CLUSTERS
install.packages("circlize")
library(circlize)
install.packages("scales")
library(scales)

#getCirclize(): generate a data frame to be used in a chord diagram
                #group.by = "idents": show connections between clusters

#assign normal colors to each cluster like this
grid.cols <- scales::hue_pal()(length(unique(seurat_combined@active.ident)))
names(grid.cols) <- levels(seurat_combined@active.ident)

#graph the chord diagram
circlize::chordDiagram(circles, self.link = 1, grid.col = grid.cols)
subset_SC <- subset(seurat_combined, Type == "T")
circles <- getCirclize(subset_SC, group.by = "ident")

chordDiagram(circles, self.link = 1, grid.col = grid.cols, directional = 1,
             direction.type = "arrows", linke.arr.type = "big.arrow")