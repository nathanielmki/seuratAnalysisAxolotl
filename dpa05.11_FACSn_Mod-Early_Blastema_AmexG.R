#Load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(RColorBrewer)

setwd("/compbio/analysis/PrayagMurawala/")
# annotation file
anno = read.delim("/compbio/analysis/PrayagMurawala/AnnotationFile_AmexG_v6_chr_unscaffolded_CherryGFP_v1.1.csv",sep = "\t",header = F)
rownames(anno) = anno[,1]

unannotated.transcripts = rownames(anno)[anno$V2 %in% anno$V2[grep("PEG10",anno$V2)] ]
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("L1TD1",anno$V2)] ], unannotated.transcripts))
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("RTL",anno$V2)] ], unannotated.transcripts))
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("GIN1",anno$V2)] ], unannotated.transcripts))
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("L1\\-RT",anno$V2)] ], unannotated.transcripts))
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("^N\\/A",anno$V2)] ], unannotated.transcripts))
unannotated.transcripts = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("^AMEX",anno$V2)] ], unannotated.transcripts))

anno = anno[!rownames(anno) %in% unannotated.transcripts,]


rp.genes_AmexG = rownames(anno)[anno$V2 %in% anno$V2[grep("RPL",anno$V2)]]
rp.genes_AmexG = unique(c(rownames(anno)[anno$V2 %in% anno$V2[grep("RPS",anno$V2)] ], rp.genes_AmexG))

anno = anno[!rownames(anno) %in% rp.genes_AmexG,]

library(tidyverse)

dpa05_FACSn_GER043 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER043/GER043_11dpa_limbBL_Solo.out/Gene/raw/")
dpa05_FACSn_GER043 = dpa05_FACSn_GER043[rownames(anno), ,drop = FALSE]

dpa05_FACSn_GER045 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER045/GER045_11dpa_limbBL_Solo.out/Gene/raw/")
dpa05_FACSn_GER045 = dpa05_FACSn_GER045[rownames(anno), ,drop = FALSE]

dpa11_FACSn_GER047 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER047/GER047_11dpa_limbBL_Solo.out/Gene/raw/")
dpa11_FACSn_GER047 = dpa11_FACSn_GER047[rownames(anno), ,drop = FALSE]

library(dplyr)
library(Seurat)
#library(harmony)

#load in Seurat

dpa05_FACSn_GER043 <- CreateSeuratObject(counts = dpa05_FACSn_GER043, project = "dpa05_FACSn_GER043", min.cells = 3, min.features = 200)
dpa05_FACSn_GER043
dpa05_FACSn_GER045 <- CreateSeuratObject(counts = dpa05_FACSn_GER045, project = "dpa05_FACSn_GER045", min.cells = 3, min.features = 200)
dpa05_FACSn_GER045
dpa11_FACSn_GER047 <- CreateSeuratObject(counts = dpa11_FACSn_GER047, project = "dpa11_FACSn_GER047", min.cells = 3, min.features = 200)
dpa11_FACSn_GER047


#merge ojects
#merge(x = NULL, y = NULL, add.cell.ids = NULL, merge.data = TRUE, project = "SeuratProject", ...)

BL_dpa05.11 = merge(dpa05_FACSn_GER043, y = c(dpa05_FACSn_GER045, dpa11_FACSn_GER047), add.cell.ids = c("GER043", "GER045", "GER047"), project = "dpa05.11_Limb_BL")

mito.genes <- c("ND2","ND1","ND3","ND4","ND4L","ND5","ND6")
mito.contig <- intersect(c(rownames(anno)[anno$V2 %in% mito.genes ] , rownames(anno)[anno$V2 %in% anno$V2[grep("^COX",anno$V2)]] ) , rownames(BL_dpa05.11))


BL_dpa05.11[["percent.mt"]] <- PercentageFeatureSet(BL_dpa05.11, features = mito.contig)


#parallize scaling
library(future)
plan("multiprocess", workers = 24)
plan()

#set maximum to 50GB for each worker
options(future.globals.maxSize = 50000 * 1024^2)

#normalize data

BL_dpa05.11 = NormalizeData(BL_dpa05.11)

plot1 <- FeatureScatter(BL_dpa05.11, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_dpa05.11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


pdf("./BL_dpa05.11_QC_PreFilter.pdf",width=10,height=5)
VlnPlot(BL_dpa05.11, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

BL_dpa05.11 <- subset(BL_dpa05.11, subset = nCount_RNA < 20000  & nCount_RNA > 500  & percent.mt < 10 )

plot1 <- FeatureScatter(BL_dpa05.11, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_dpa05.11, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("./BL_dpa05.11_QC_PostFilter.pdf",width=10,height=5)
VlnPlot(BL_dpa05.11, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#get cell cycle genes

g2m.genes <- cc.genes$g2m.genes
g2m.contig <- intersect(rownames(anno)[anno$V2 %in% g2m.genes ] , rownames(BL_dpa05.11) )

s.genes <- cc.genes$s.genes
s.contig <- intersect(rownames(anno)[anno$V2 %in% s.genes ] , rownames(BL_dpa05.11) )

#cell cycle scoring
BL_dpa05.11 <- CellCycleScoring(BL_dpa05.11, s.features = s.contig, g2m.features = g2m.contig, set.ident = TRUE)


#scale
#ScaleData(  object,  features = NULL,  assay = NULL,  vars.to.regress = NULL,  split.by = NULL,  model.use = "linear",  use.umi = FALSE,  do.scale = TRUE,  do.center = TRUE,  scale.max = 10,  block.size = 1000,  min.cells.to.block = 3000,  verbose = TRUE)

all.genes <- rownames(BL_dpa05.11)

BL_dpa05.11 = ScaleData(  BL_dpa05.11, features = all.genes ,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","S.Score", "G2M.Score","eGFP","mCherry"))

#run PCA on all genes
#RunPCA(object, assay = NULL, features = NULL,  npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,  ndims.print = 1:5, nfeatures.print = 30, reduction.name = "pca",  reduction.key = "PC_", seed.use = 42)

BL_dpa05.11 = RunPCA(BL_dpa05.11,features = all.genes,npcs = 100)

#plot PCA heatmaps
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

library(RColorBrewer)

pdf("BL_dpa05.11_PCAheatmaps.pdf",width=30,height=20)
DimHeatmap(BL_dpa05.11, dims = 1:16, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_dpa05.11, dims = 17:32, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_dpa05.11, dims = 33:48, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa05.11, dims = 49:64, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa05.11, dims = 65:80, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa05.11, dims = 81:96, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
dev.off()

#elbow plot
pdf("BL_dpa05.11_PCelbow.pdf",width=20,height=10)
ElbowPlot(BL_dpa05.11, ndims = 100)
dev.off()

#Cluster cells
BL_dpa05.11<- FindNeighbors(BL_dpa05.11, dims = 1:30)
BL_dpa05.11 <- FindClusters(BL_dpa05.11, resolution = 0.5)

#run UMAP
#min.dist deafult is 0.3. Can go down to 0.001. Play around to modify clustering
BL_dpa05.11 = RunUMAP(BL_dpa05.11,dims = 1:30)

#run tSNE
BL_dpa05.11 = RunTSNE(BL_dpa05.11,dims = 1:30)

#Subset Clusters
#Subset on 5,17

Cluster.5 <- subset(x = BL_dpa05.11_SeuratObj, subset = cluster == "5")

#run Harmony
#BL_dpa05.11 <- RunHarmony(BL_dpa05.11_SeuratObj, "dataset")
#BL_dpa05.11 <- RunUMAP(BL_dpa05.11_SeuratObj, reduction = "harmony", dims = 1:30)


#plot result
#DimPlot(object, dims = c(1, 2), cells = NULL, cols = NULL,  pt.size = NULL, reduction = NULL, group.by = NULL,  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,  label.size = 4, repel = FALSE, cells.highlight = NULL,  cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",  combine = TRUE)

pdf("BL_dpa05.11_Embedding_PC30_res0.5.pdf",width=10,height=10)
DimPlot(object = BL_dpa05.11, reduction = 'umap', pt.size = 2)
DimPlot(object = BL_dpa05.11, reduction = 'umap', pt.size = 2,group.by = "orig.ident")
DimPlot(object = BL_dpa05.11, reduction = 'tsne', pt.size = 2)
dev.off()

#some QC features

pdf("BL_dpa05.11_feature.pdf",width=8,height=10)
FeaturePlot(BL_dpa05.11, reduction = 'umap', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
FeaturePlot(BL_dpa05.11, reduction = 'tsne', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
dev.off()


plan("multiprocess", workers = 6)
library(tidyverse)

markers <- FindAllMarkers(BL_dpa05.11, only.pos = TRUE,  logfc.threshold = 0.3)
markers.anno = markers
markers.anno = merge(markers.anno,anno, by.x="gene" , by.y="V1")
markers.anno$ID = markers.anno$gene

#markers.anno = markers.anno[order(as.numeric(markers.anno$cluster)),]
markers.anno = markers.anno  %>% arrange(cluster , desc(avg_logFC))

write.csv(markers.anno,"BL_dpa05.11_allMarker.csv")

saveRDS(BL_dpa05.11, file = "BL_dpa05.11_SeuratObj.RDS")


#plot some canonical markers to roughly identfy cell types


gene_ids = c("AMEX60DD027986","AMEX60DD056342","AMEX60DD045921","AMEX60DD035908","AMEX60DD022398","AMEX60DD018450","AMEX60DD020580","AMEX60DD024035","AMEX60DD006619","AMEX60DD013910","AMEX60DD042097","AMEX60DD043936","AMEX60DD025155","AMEX60DD052070","AMEX60DD031414","AMEX60DD025537","AMEX60DD032898", "eGFP", "mCherry")

gene_names = c("MYLPF","DES","S100P","EPCAM","COL1A2","PRRX1","MYH11","GP9","VWF","PLVAP","GZMA","PAX5","CTSW","C1QB","LGALS3BP","ALAS2","RHAG", "eGFP", "mCherry")

gg_Fig <- FeaturePlot(BL_dpa05.11, pt.size = 0.5, features = gene_ids,order = T, cols = brewer.pal(9,"Greys")[3:9], repel=T)
gg_Fig <- lapply( 1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=gene_names[x]) })

pdf("BL_dpa05.11_UMAP_feature_ClustMarker.pdf",width=14,height=10)
CombinePlots( gg_Fig )
dev.off()

# MYLPF	AMEX60DD027986		Muscle
# DES	AMEX60DD056342		Muscle
# S100P	AMEX60DD045921		Epidermis
# EPCAM	AMEX60DD035908		Epidermis
# COL1A2	AMEX60DD022398		CT
# PRRX1	AMEX60DD018450		CT
# PHOSPHO1	AMEX60DD009968		Bone
# MYH11	AMEX60DD020580		Pericyte
# GP9	AMEX60DD024035		Pericyte
# VWF	AMEX60DD006619		Endothel
# PLVAP	AMEX60DD013910		Endothel
# GZMA	AMEX60DD042097		Immune
# PAX5	AMEX60DD043936		Immune
# CTSW	AMEX60DD025155    Immune
# C1QB	AMEX60DD052070		Macrophage
# LGALS3BP	AMEX60DD031414		Macrophage
# ALAS2	AMEX60DD025537		Red blood cell
# RHAG	AMEX60DD032898		Red blood cell
