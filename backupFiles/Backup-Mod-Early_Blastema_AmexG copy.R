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

BL_5dpa_A <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER004/GER004_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_A = BL_5dpa_A[rownames(anno), ,drop = FALSE]

BL_5dpa_B <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER005/GER005_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_B = BL_5dpa_B[rownames(anno), ,drop = FALSE]

BL_5dpa_C <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER006/GER006_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_C = BL_5dpa_C[rownames(anno), ,drop = FALSE]

BL_5dpa_D <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER007/GER007_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_D = BL_5dpa_D[rownames(anno), ,drop = FALSE]

BL_5dpa_E <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER014/GER014_11dpa_limbBL_L007_Solo.out/Gene/raw/")
BL_5dpa_E = BL_5dpa_E[rownames(anno), ,drop = FALSE]

BL_5dpa_F <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER015/GER015_11dpa_limbBL_L008_Solo.out/Gene/raw/")
BL_5dpa_F = BL_5dpa_F[rownames(anno), ,drop = FALSE]

BL_5dpa_G <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER017/GER017_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_G = BL_5dpa_G[rownames(anno), ,drop = FALSE]

BL_5dpa_H <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER018/GER018_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_H = BL_5dpa_H[rownames(anno), ,drop = FALSE]

BL_5dpa_I1 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER019/GER019_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_I1 = BL_5dpa_I1[rownames(anno), ,drop = FALSE]

BL_5dpa_I2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER019/GER019_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_I2 = BL_5dpa_I2[rownames(anno), ,drop = FALSE]

BL_5dpa_J1 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER020/GER020_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_J1 = BL_5dpa_J1[rownames(anno), ,drop = FALSE]

BL_5dpa_J2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER020/GER020_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_J2 = BL_5dpa_J2[rownames(anno), ,drop = FALSE]

BL_5dpa_K <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER024/GER024_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_K = BL_5dpa_K[rownames(anno), ,drop = FALSE]

BL_5dpa_L1 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER025/GER025_11dpa_limbBL_L001_Solo.out/Gene/raw/")
BL_5dpa_L1 = BL_5dpa_L1[rownames(anno), ,drop = FALSE]

BL_5dpa_L2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER025/GER025_11dpa_limbBL_L002_Solo.out/Gene/raw/")
BL_5dpa_L2 = BL_5dpa_L2[rownames(anno), ,drop = FALSE]

BL_5dpa_M <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER043/GER043_11dpa_limbBL_Solo.out/Gene/raw/")
BL_5dpa_M = BL_5dpa_M[rownames(anno), ,drop = FALSE]

BL_5dpa_N <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER044/GER044_11dpa_limbBL_Solo.out/Gene/raw/")
BL_5dpa_N = BL_5dpa_N[rownames(anno), ,drop = FALSE]

BL_5dpa_O <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER045/GER045_11dpa_limbBL_Solo.out/Gene/raw/")
BL_5dpa_O = BL_5dpa_O[rownames(anno), ,drop = FALSE]

BL_5dpa_P <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER046/GER046_11dpa_limbBL_Solo.out/Gene/raw/")
BL_5dpa_P = BL_5dpa_P[rownames(anno), ,drop = FALSE]

BL_5dpa_Q <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER047/GER047_11dpa_limbBL_Solo.out/Gene/raw/")
BL_5dpa_Q = BL_5dpa_Q[rownames(anno), ,drop = FALSE]

library(dplyr)
library(Seurat)

#load in Seurat

BL_5dpa_A <- CreateSeuratObject(counts = BL_5dpa_A, project = "BL_5dpa_A", min.cells = 3, min.features = 200)
BL_5dpa_A
BL_5dpa_B <- CreateSeuratObject(counts = BL_5dpa_B, project = "BL_5dpa_B", min.cells = 3, min.features = 200)
BL_5dpa_B
BL_5dpa_C <- CreateSeuratObject(counts = BL_5dpa_C, project = "BL_5dpa_C", min.cells = 3, min.features = 200)
BL_5dpa_C
BL_5dpa_D <- CreateSeuratObject(counts = BL_5dpa_D, project = "BL_5dpa_D", min.cells = 3, min.features = 200)
BL_5dpa_D
BL_5dpa_E <- CreateSeuratObject(counts = BL_5dpa_E, project = "BL_5dpa_E", min.cells = 3, min.features = 200)
BL_5dpa_E
BL_5dpa_F <- CreateSeuratObject(counts = BL_5dpa_F, project = "BL_5dpa_F", min.cells = 3, min.features = 200)
BL_5dpa_F
BL_5dpa_G <- CreateSeuratObject(counts = BL_5dpa_G, project = "BL_5dpa_G", min.cells = 3, min.features = 200)
BL_5dpa_G
BL_5dpa_H <- CreateSeuratObject(counts = BL_5dpa_H, project = "BL_5dpa_H", min.cells = 3, min.features = 200)
BL_5dpa_H
BL_5dpa_I1 <- CreateSeuratObject(counts = BL_5dpa_I1, project = "BL_5dpa_I1", min.cells = 3, min.features = 200)
BL_5dpa_I1
BL_5dpa_I2 <- CreateSeuratObject(counts = BL_5dpa_I2, project = "BL_5dpa_I2", min.cells = 3, min.features = 200)
BL_5dpa_I2
BL_5dpa_J1 <- CreateSeuratObject(counts = BL_5dpa_J1, project = "BL_5dpa_J1", min.cells = 3, min.features = 200)
BL_5dpa_J1
BL_5dpa_J2 <- CreateSeuratObject(counts = BL_5dpa_J2, project = "BL_5dpa_J2", min.cells = 3, min.features = 200)
BL_5dpa_J2
BL_5dpa_K <- CreateSeuratObject(counts = BL_5dpa_K, project = "BL_5dpa_K", min.cells = 3, min.features = 200)
BL_5dpa_K
BL_5dpa_L1 <- CreateSeuratObject(counts = BL_5dpa_L1, project = "BL_5dpa_L1", min.cells = 3, min.features = 200)
BL_5dpa_L1
BL_5dpa_L2 <- CreateSeuratObject(counts = BL_5dpa_L2, project = "BL_5dpa_L2", min.cells = 3, min.features = 200)
BL_5dpa_L2
BL_5dpa_M <- CreateSeuratObject(counts = BL_5dpa_M, project = "BL_5dpa_M", min.cells = 3, min.features = 200)
BL_5dpa_M
BL_5dpa_N <- CreateSeuratObject(counts = BL_5dpa_N, project = "BL_5dpa_N", min.cells = 3, min.features = 200)
BL_5dpa_N
BL_5dpa_O <- CreateSeuratObject(counts = BL_5dpa_O, project = "BL_5dpa_O", min.cells = 3, min.features = 200)
BL_5dpa_O
BL_5dpa_P <- CreateSeuratObject(counts = BL_5dpa_P, project = "BL_5dpa_P", min.cells = 3, min.features = 200)
BL_5dpa_P
BL_5dpa_Q <- CreateSeuratObject(counts = BL_5dpa_Q, project = "BL_5dpa_Q", min.cells = 3, min.features = 200)
BL_5dpa_Q


#merge ojects
#merge(x = NULL, y = NULL, add.cell.ids = NULL, merge.data = TRUE, project = "SeuratProject", ...)

BL_5dpa = merge(BL_5dpa_A, y = c(BL_5dpa_B, BL_5dpa_C, BL_5dpa_D, BL_5dpa_E, BL_5dpa_F, BL_5dpa_G, BL_5dpa_H, BL_5dpa_I1, BL_5dpa_I2, BL_5dpa_J1, BL_5dpa_J2, BL_5dpa_K, BL_5dpa_L1, BL_5dpa_L2, BL_5dpa_M, BL_5dpa_N, BL_5dpa_O, BL_5dpa_P, BL_5dpa_Q), add.cell.ids = c("A", "B", "C", "D", "E", "F", "G", "H", "I1", "I2", "J1", "J2", "K", "L1", "L2", "M", "N", "O", "P", "Q"), project = "Limb_BL_5dpa")
#BL_5dpa = merge(BL_5dpa_M, y = c(BL_5dpa_N, BL_5dpa_O, BL_5dpa_P, BL_5dpa_Q), add.cell.ids = c("M", "N", "O", "P", "Q"), project = "Limb_BL_5dpa")

mito.genes <- c("ND2","ND1","ND3","ND4","ND4L","ND5","ND6")
mito.contig <- intersect(c(rownames(anno)[anno$V2 %in% mito.genes ] , rownames(anno)[anno$V2 %in% anno$V2[grep("^COX",anno$V2)]] ) , rownames(BL_5dpa))


BL_5dpa[["percent.mt"]] <- PercentageFeatureSet(BL_5dpa, features = mito.contig)


#parallize scaling
library(future)
plan("multiprocess", workers = 50)
plan()

#set maximum to 50GB for each worker
options(future.globals.maxSize = 50000 * 1024^2)

#normalize data

BL_5dpa = NormalizeData(BL_5dpa)

plot1 <- FeatureScatter(BL_5dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_5dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


pdf("./BL_5dpa_QC_PreFilter.pdf",width=10,height=5)
VlnPlot(BL_5dpa, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

BL_5dpa <- subset(BL_5dpa, subset = nCount_RNA < 20000  & nCount_RNA > 500  & percent.mt < 10 )

plot1 <- FeatureScatter(BL_5dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_5dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("./BL_5dpa_QC_PostFilter.pdf",width=10,height=5)
VlnPlot(BL_5dpa, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#get cell cycle genes

g2m.genes <- cc.genes$g2m.genes
g2m.contig <- intersect(rownames(anno)[anno$V2 %in% g2m.genes ] , rownames(BL_5dpa) )

s.genes <- cc.genes$s.genes
s.contig <- intersect(rownames(anno)[anno$V2 %in% s.genes ] , rownames(BL_5dpa) )

#cell cycle scoring
BL_5dpa <- CellCycleScoring(BL_5dpa, s.features = s.contig, g2m.features = g2m.contig, set.ident = TRUE)


#scale
#ScaleData(  object,  features = NULL,  assay = NULL,  vars.to.regress = NULL,  split.by = NULL,  model.use = "linear",  use.umi = FALSE,  do.scale = TRUE,  do.center = TRUE,  scale.max = 10,  block.size = 1000,  min.cells.to.block = 3000,  verbose = TRUE)

all.genes <- rownames(BL_5dpa)

BL_5dpa = ScaleData(  BL_5dpa, features = all.genes ,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","S.Score", "G2M.Score","eGFP","mCherry"))

#run PCA on all genes
#RunPCA(object, assay = NULL, features = NULL,  npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,  ndims.print = 1:5, nfeatures.print = 30, reduction.name = "pca",  reduction.key = "PC_", seed.use = 42)

BL_5dpa = RunPCA(BL_5dpa,features = all.genes,npcs = 100)

#plot PCA heatmaps
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

library(RColorBrewer)

pdf("BL_5dpa_PCAheatmaps.pdf",width=30,height=20)
DimHeatmap(BL_5dpa, dims = 1:16, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_5dpa, dims = 17:32, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_5dpa, dims = 33:48, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_5dpa, dims = 49:64, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_5dpa, dims = 65:80, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_5dpa, dims = 81:96, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
dev.off()

#elbow plot
pdf("BL_5dpa_PCelbow.pdf",width=20,height=10)
ElbowPlot(BL_5dpa, ndims = 100)
dev.off()

#Cluster cells
BL_5dpa<- FindNeighbors(BL_5dpa, dims = 1:30)
BL_5dpa <- FindClusters(BL_5dpa, resolution = 0.5)

#run UMAP
#min.dist deafult is 0.3. Can go down to 0.001. Play around to modify clustering
BL_5dpa = RunUMAP(BL_5dpa,dims = 1:30)

#run tSNE
BL_5dpa = RunTSNE(BL_5dpa,dims = 1:30)


#plot result
#DimPlot(object, dims = c(1, 2), cells = NULL, cols = NULL,  pt.size = NULL, reduction = NULL, group.by = NULL,  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,  label.size = 4, repel = FALSE, cells.highlight = NULL,  cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",  combine = TRUE)

pdf("BL_5dpa_Embedding_PC30_res0.5.pdf",width=10,height=10)
DimPlot(object = BL_5dpa, reduction = 'umap', pt.size = 2)
DimPlot(object = BL_5dpa, reduction = 'umap', pt.size = 2,group.by = "orig.ident")
DimPlot(object = BL_5dpa, reduction = 'tsne', pt.size = 2)
dev.off()

#some QC features

pdf("BL_5dpa_feature.pdf",width=8,height=10)
FeaturePlot(BL_5dpa, reduction = 'umap', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
FeaturePlot(BL_5dpa, reduction = 'tsne', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
dev.off()


plan("multiprocess", workers = 6)
library(tidyverse)

markers <- FindAllMarkers(BL_5dpa, only.pos = TRUE,  logfc.threshold = 0.3)
markers.anno = markers
markers.anno = merge(markers.anno,anno, by.x="gene" , by.y="V1")
markers.anno$ID = markers.anno$gene

#markers.anno = markers.anno[order(as.numeric(markers.anno$cluster)),]
markers.anno = markers.anno  %>% arrange(cluster , desc(avg_logFC))

write.csv(markers.anno,"BL_5dpa_allMarker.csv")

saveRDS(BL_5dpa, file = "BL_5dpa_SeuratObj.RDS")


#plot some canonical markers to roughly identfy cell types


gene_ids = c("AMEX60DD027986","AMEX60DD056342","AMEX60DD045921","AMEX60DD035908","AMEX60DD022398","AMEX60DD018450","AMEX60DD020580","AMEX60DD024035","AMEX60DD006619","AMEX60DD013910","AMEX60DD042097","AMEX60DD043936","AMEX60DD025155","AMEX60DD052070","AMEX60DD031414","AMEX60DD025537","AMEX60DD032898")

gene_names = c("MYLPF","DES","S100P","EPCAM","COL1A2","PRRX1","MYH11","GP9","VWF","PLVAP","GZMA","PAX5","CTSW","C1QB","LGALS3BP","ALAS2","RHAG")

gg_Fig <- FeaturePlot(BL_5dpa, pt.size = 0.5, features = gene_ids,order = T, cols = brewer.pal(9,"Greys")[3:9], repel=T)
gg_Fig <- lapply( 1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=gene_names[x]) })

pdf("BL_5dpa_UMAP_feature_ClustMarker.pdf",width=14,height=10)
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
