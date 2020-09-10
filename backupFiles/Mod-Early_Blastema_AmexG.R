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

dpa18_FACSy_GER004 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER004/GER004_11dpa_limbBL_L001_Solo.out/Gene/raw/")
dpa18_FACSy_GER004 = dpa18_FACSy_GER004[rownames(anno), ,drop = FALSE]

dpa18_FACSy_GER005 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER005/GER005_11dpa_limbBL_L002_Solo.out/Gene/raw/")
dpa18_FACSy_GER005 = dpa18_FACSy_GER005[rownames(anno), ,drop = FALSE]

dpa00_FACSn_GER006 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER006/GER006_11dpa_limbBL_L001_Solo.out/Gene/raw/")
dpa00_FACSn_GER006 = dpa00_FACSn_GER006[rownames(anno), ,drop = FALSE]

dpa00_FACSn_GER007 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER007/GER007_11dpa_limbBL_L002_Solo.out/Gene/raw/")
dpa00_FACSn_GER007 = dpa00_FACSn_GER007[rownames(anno), ,drop = FALSE]

dpa00_FACSy_GER014 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER014/GER014_11dpa_limbBL_L007_Solo.out/Gene/raw/")
dpa00_FACSy_GER014 = dpa00_FACSy_GER014[rownames(anno), ,drop = FALSE]

dpa00_FACSy_GER015 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER015/GER015_11dpa_limbBL_L008_Solo.out/Gene/raw/")
dpa00_FACSy_GER015 = dpa00_FACSy_GER015[rownames(anno), ,drop = FALSE]

dpa18_FACSy_GER017 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER017/GER017_11dpa_limbBL_L001_Solo.out/Gene/raw/")
dpa18_FACSy_GER017 = dpa18_FACSy_GER017[rownames(anno), ,drop = FALSE]

dpa18_FACSy_GER018 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER018/GER018_11dpa_limbBL_L002_Solo.out/Gene/raw/")
dpa18_FACSy_GER018 = dpa18_FACSy_GER018[rownames(anno), ,drop = FALSE]

dpa25_FACSy_GER019 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER019/GER019_11dpa_limbBL_L012_Solo.out/Gene/raw/")
dpa25_FACSy_GER019 = dpa25_FACSy_GER019[rownames(anno), ,drop = FALSE]

#dpa25_FACSy_GER019_l2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER019/GER019_11dpa_limbBL_L002_Solo.out/Gene/raw/")
#dpa25_FACSy_GER019_l2 = dpa25_FACSy_GER019_l2[rownames(anno), ,drop = FALSE]

dpa25_FACSy_GER020 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER020/GER020_11dpa_limbBL_L012_Solo.out/Gene/raw/")
dpa25_FACSy_GER020 = dpa25_FACSy_GER020[rownames(anno), ,drop = FALSE]

#dpa25_FACSy_GER020_l2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER020/GER020_11dpa_limbBL_L002_Solo.out/Gene/raw/")
#dpa25_FACSy_GER020_l2 = dpa25_FACSy_GER020_l2[rownames(anno), ,drop = FALSE]

dpa38_FACSy_GER024 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER024/GER024_11dpa_limbBL_L001_Solo.out/Gene/raw/")
dpa38_FACSy_GER024 = dpa38_FACSy_GER024[rownames(anno), ,drop = FALSE]

dpa38_FACSy_GER025 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER025/GER025_11dpa_limbBL_L012_Solo.out/Gene/raw/")
dpa38_FACSy_GER025 = dpa38_FACSy_GER025[rownames(anno), ,drop = FALSE]

#dpa38_FACSy_GER025_l2 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER025/GER025_11dpa_limbBL_L002_Solo.out/Gene/raw/")
#dpa38_FACSy_GER025_l2 = dpa38_FACSy_GER025_l2[rownames(anno), ,drop = FALSE]

dpa05_FACSn_GER043 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER043/GER043_11dpa_limbBL_Solo.out/Gene/raw/")
dpa05_FACSn_GER043 = dpa05_FACSn_GER043[rownames(anno), ,drop = FALSE]

dpa05_FACSy_GER044 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER044/GER044_11dpa_limbBL_Solo.out/Gene/raw/")
dpa05_FACSy_GER044 = dpa05_FACSy_GER044[rownames(anno), ,drop = FALSE]

dpa05_FACSn_GER045 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER045/GER045_11dpa_limbBL_Solo.out/Gene/raw/")
dpa05_FACSn_GER045 = dpa05_FACSn_GER045[rownames(anno), ,drop = FALSE]

dpa11_FACSy_GER046 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER046/GER046_11dpa_limbBL_Solo.out/Gene/raw/")
dpa11_FACSy_GER046 = dpa11_FACSy_GER046[rownames(anno), ,drop = FALSE]

dpa11_FACSn_GER047 <- Read10X(data.dir = "/compbio/analysis/PrayagMurawala/GER047/GER047_11dpa_limbBL_Solo.out/Gene/raw/")
dpa11_FACSn_GER047 = dpa11_FACSn_GER047[rownames(anno), ,drop = FALSE]

library(dplyr)
library(Seurat)

#load in Seurat

dpa18_FACSy_GER004 <- CreateSeuratObject(counts = dpa18_FACSy_GER004, project = "dpa18_FACSy_GER004", min.cells = 3, min.features = 200)
dpa18_FACSy_GER004
dpa18_FACSy_GER005 <- CreateSeuratObject(counts = dpa18_FACSy_GER005, project = "dpa18_FACSy_GER005", min.cells = 3, min.features = 200)
dpa18_FACSy_GER005
dpa00_FACSn_GER006 <- CreateSeuratObject(counts = dpa00_FACSn_GER006, project = "dpa00_FACSn_GER006", min.cells = 3, min.features = 200)
dpa00_FACSn_GER006
dpa00_FACSn_GER007 <- CreateSeuratObject(counts = dpa00_FACSn_GER007, project = "dpa00_FACSn_GER007", min.cells = 3, min.features = 200)
dpa00_FACSn_GER007
dpa00_FACSy_GER014 <- CreateSeuratObject(counts = dpa00_FACSy_GER014, project = "dpa00_FACSy_GER014", min.cells = 3, min.features = 200)
dpa00_FACSy_GER014
dpa00_FACSy_GER015 <- CreateSeuratObject(counts = dpa00_FACSy_GER015, project = "dpa00_FACSy_GER015", min.cells = 3, min.features = 200)
dpa00_FACSy_GER015
dpa18_FACSy_GER017 <- CreateSeuratObject(counts = dpa18_FACSy_GER017, project = "dpa18_FACSy_GER017", min.cells = 3, min.features = 200)
dpa18_FACSy_GER017
dpa18_FACSy_GER018 <- CreateSeuratObject(counts = dpa18_FACSy_GER018, project = "dpa18_FACSy_GER018", min.cells = 3, min.features = 200)
dpa18_FACSy_GER018
dpa25_FACSy_GER019 <- CreateSeuratObject(counts = dpa25_FACSy_GER019, project = "dpa25_FACSy_GER019", min.cells = 3, min.features = 200)
dpa25_FACSy_GER019
#dpa25_FACSy_GER019_l2 <- CreateSeuratObject(counts = dpa25_FACSy_GER019_l2, project = "dpa25_FACSy_GER019_l2", min.cells = 3, min.features = 200)
#dpa25_FACSy_GER019_l2
dpa25_FACSy_GER020 <- CreateSeuratObject(counts = dpa25_FACSy_GER020, project = "dpa25_FACSy_GER020", min.cells = 3, min.features = 200)
dpa25_FACSy_GER020
#dpa25_FACSy_GER020_l2 <- CreateSeuratObject(counts = dpa25_FACSy_GER020_l2, project = "dpa25_FACSy_GER020_l2", min.cells = 3, min.features = 200)
#da25_FACSy_GER020_l2
dpa38_FACSy_GER024 <- CreateSeuratObject(counts = dpa38_FACSy_GER024, project = "dpa38_FACSy_GER024", min.cells = 3, min.features = 200)
dpa38_FACSy_GER024
dpa38_FACSy_GER025 <- CreateSeuratObject(counts = dpa38_FACSy_GER025, project = "dpa38_FACSy_GER025", min.cells = 3, min.features = 200)
dpa38_FACSy_GER025
#dpa38_FACSy_GER025_l2 <- CreateSeuratObject(counts = dpa38_FACSy_GER025_l2, project = "dpa38_FACSy_GER025_l2", min.cells = 3, min.features = 200)
#dpa38_FACSy_GER025_l2
dpa05_FACSn_GER043 <- CreateSeuratObject(counts = dpa05_FACSn_GER043, project = "dpa05_FACSn_GER043", min.cells = 3, min.features = 200)
dpa05_FACSn_GER043
dpa05_FACSy_GER044 <- CreateSeuratObject(counts = dpa05_FACSy_GER044, project = "dpa05_FACSy_GER044", min.cells = 3, min.features = 200)
dpa05_FACSy_GER044
dpa05_FACSn_GER045 <- CreateSeuratObject(counts = dpa05_FACSn_GER045, project = "dpa05_FACSn_GER045", min.cells = 3, min.features = 200)
dpa05_FACSn_GER045
dpa11_FACSy_GER046 <- CreateSeuratObject(counts = dpa11_FACSy_GER046, project = "dpa11_FACSy_GER046", min.cells = 3, min.features = 200)
dpa11_FACSy_GER046
dpa11_FACSn_GER047 <- CreateSeuratObject(counts = dpa11_FACSn_GER047, project = "dpa11_FACSn_GER047", min.cells = 3, min.features = 200)
dpa11_FACSn_GER047


#merge ojects
#merge(x = NULL, y = NULL, add.cell.ids = NULL, merge.data = TRUE, project = "SeuratProject", ...)

BL_dpa = merge(dpa18_FACSy_GER004, y = c(dpa18_FACSy_GER005, dpa00_FACSn_GER006, dpa00_FACSn_GER007, dpa00_FACSy_GER014, dpa00_FACSy_GER015, dpa18_FACSy_GER017, dpa18_FACSy_GER018, dpa25_FACSy_GER019, dpa25_FACSy_GER020, dpa38_FACSy_GER024, dpa38_FACSy_GER025, dpa05_FACSn_GER043, dpa05_FACSy_GER044, dpa05_FACSn_GER045, dpa11_FACSy_GER046, dpa11_FACSn_GER047), add.cell.ids = c("GER004", "GER005", "GER006", "GER007", "GER014", "GER015", "GER017", "GER018", "GER019", "GER020", "GER024", "GER025", "GER043", "GER044", "GER045", "GER046", "GER047"), project = "Limb_BL_dpa")

mito.genes <- c("ND2","ND1","ND3","ND4","ND4L","ND5","ND6")
mito.contig <- intersect(c(rownames(anno)[anno$V2 %in% mito.genes ] , rownames(anno)[anno$V2 %in% anno$V2[grep("^COX",anno$V2)]] ) , rownames(BL_dpa))


BL_dpa[["percent.mt"]] <- PercentageFeatureSet(BL_dpa, features = mito.contig)


#parallize scaling
library(future)
plan("multiprocess", workers = 50)
plan()

#set maximum to 50GB for each worker
options(future.globals.maxSize = 50000 * 1024^2)

#normalize data

BL_dpa = NormalizeData(BL_dpa)

plot1 <- FeatureScatter(BL_dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")


pdf("./BL_dpa_QC_PreFilter.pdf",width=10,height=5)
VlnPlot(BL_dpa, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()

BL_dpa <- subset(BL_dpa, subset = nCount_RNA < 20000  & nCount_RNA > 500  & percent.mt < 10 )

plot1 <- FeatureScatter(BL_dpa, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(BL_dpa, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

pdf("./BL_dpa_QC_PostFilter.pdf",width=10,height=5)
VlnPlot(BL_dpa, features = c("nCount_RNA","nFeature_RNA","percent.mt"),pt.size = -1)
CombinePlots(plots = list(plot1, plot2))
dev.off()


#get cell cycle genes

g2m.genes <- cc.genes$g2m.genes
g2m.contig <- intersect(rownames(anno)[anno$V2 %in% g2m.genes ] , rownames(BL_dpa) )

s.genes <- cc.genes$s.genes
s.contig <- intersect(rownames(anno)[anno$V2 %in% s.genes ] , rownames(BL_dpa) )

#cell cycle scoring
BL_dpa <- CellCycleScoring(BL_dpa, s.features = s.contig, g2m.features = g2m.contig, set.ident = TRUE)


#scale
#ScaleData(  object,  features = NULL,  assay = NULL,  vars.to.regress = NULL,  split.by = NULL,  model.use = "linear",  use.umi = FALSE,  do.scale = TRUE,  do.center = TRUE,  scale.max = 10,  block.size = 1000,  min.cells.to.block = 3000,  verbose = TRUE)

all.genes <- rownames(BL_dpa)

BL_dpa = ScaleData(  BL_dpa, features = all.genes ,vars.to.regress = c("nCount_RNA","nFeature_RNA","percent.mt","S.Score", "G2M.Score","eGFP","mCherry"))

#run PCA on all genes
#RunPCA(object, assay = NULL, features = NULL,  npcs = 50, rev.pca = FALSE, weight.by.var = TRUE, verbose = TRUE,  ndims.print = 1:5, nfeatures.print = 30, reduction.name = "pca",  reduction.key = "PC_", seed.use = 42)

BL_dpa = RunPCA(BL_dpa,features = all.genes,npcs = 100)

#plot PCA heatmaps
#DimHeatmap(pbmc, dims = 1, cells = 500, balanced = TRUE)

library(RColorBrewer)

pdf("BL_dpa_PCAheatmaps.pdf",width=30,height=20)
DimHeatmap(BL_dpa, dims = 1:16, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_dpa, dims = 17:32, cells = 1000, balanced = TRUE, ncol = 4,  fast = F) 
DimHeatmap(BL_dpa, dims = 33:48, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa, dims = 49:64, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa, dims = 65:80, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
DimHeatmap(BL_dpa, dims = 81:96, cells = 1000, balanced = TRUE,  ncol = 4, fast = F) 
dev.off()

#elbow plot
pdf("BL_dpa_PCelbow.pdf",width=20,height=10)
ElbowPlot(BL_dpa, ndims = 100)
dev.off()

#Cluster cells
BL_dpa<- FindNeighbors(BL_dpa, dims = 1:30)
BL_dpa <- FindClusters(BL_dpa, resolution = 0.5)

#run UMAP
#min.dist deafult is 0.3. Can go down to 0.001. Play around to modify clustering
BL_dpa = RunUMAP(BL_dpa,dims = 1:30)

#run tSNE
BL_dpa = RunTSNE(BL_dpa,dims = 1:30)


#plot result
#DimPlot(object, dims = c(1, 2), cells = NULL, cols = NULL,  pt.size = NULL, reduction = NULL, group.by = NULL,  split.by = NULL, shape.by = NULL, order = NULL, label = FALSE,  label.size = 4, repel = FALSE, cells.highlight = NULL,  cols.highlight = "red", sizes.highlight = 1, na.value = "grey50",  combine = TRUE)

pdf("BL_dpa_Embedding_PC30_res0.5.pdf",width=10,height=10)
DimPlot(object = BL_dpa, reduction = 'umap', pt.size = 2)
DimPlot(object = BL_dpa, reduction = 'umap', pt.size = 2,group.by = "orig.ident")
DimPlot(object = BL_dpa, reduction = 'tsne', pt.size = 2)
dev.off()

#some QC features

pdf("BL_dpa_feature.pdf",width=8,height=10)
FeaturePlot(BL_dpa, reduction = 'umap', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
FeaturePlot(BL_dpa, reduction = 'tsne', pt.size = 0.5, features = c("nFeature_RNA","nCount_RNA","percent.mt","mCherry","eGFP"),order = T, cols = c(brewer.pal(9,"Greys")[9:2],brewer.pal(9,"Reds")[2:9]))
dev.off()


plan("multiprocess", workers = 6)
library(tidyverse)

markers <- FindAllMarkers(BL_dpa, only.pos = TRUE,  logfc.threshold = 0.3)
markers.anno = markers
markers.anno = merge(markers.anno,anno, by.x="gene" , by.y="V1")
markers.anno$ID = markers.anno$gene

#markers.anno = markers.anno[order(as.numeric(markers.anno$cluster)),]
markers.anno = markers.anno  %>% arrange(cluster , desc(avg_logFC))

write.csv(markers.anno,"BL_dpa_allMarker.csv")

saveRDS(BL_dpa, file = "BL_dpa_SeuratObj.RDS")


#plot some canonical markers to roughly identfy cell types


gene_ids = c("AMEX60DD027986","AMEX60DD056342","AMEX60DD045921","AMEX60DD035908","AMEX60DD022398","AMEX60DD018450","AMEX60DD020580","AMEX60DD024035","AMEX60DD006619","AMEX60DD013910","AMEX60DD042097","AMEX60DD043936","AMEX60DD025155","AMEX60DD052070","AMEX60DD031414","AMEX60DD025537","AMEX60DD032898")

gene_names = c("MYLPF","DES","S100P","EPCAM","COL1A2","PRRX1","MYH11","GP9","VWF","PLVAP","GZMA","PAX5","CTSW","C1QB","LGALS3BP","ALAS2","RHAG")

gg_Fig <- FeaturePlot(BL_dpa, pt.size = 0.5, features = gene_ids,order = T, cols = brewer.pal(9,"Greys")[3:9], repel=T)
gg_Fig <- lapply( 1:length(gene_ids), function(x) { gg_Fig[[x]] + labs(title=gene_names[x]) })

pdf("BL_dpa_UMAP_feature_ClustMarker.pdf",width=14,height=10)
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
