library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(umap)
library(parallel)
library(ggplot2)
#########################################################################
#                                                                       #
#                        Single cell RNA-seq                            #
#                                                                       #
#########################################################################

#1.Identification of protoplasting-sensitive genes by bulk RNA-seq####
#Reading data and calculating the transcript per million(TPM) of each gene
#"Gene.length" is the file about the length of represented transcipt of each gene
genelength <- read.table("Gene.length",sep="\t",header = T)[,c(1,6)]
head(genelength)
str(genelength)
genelength$Geneid <- as.character(genelength$Geneid)
#The data of post-protoplasting in replicate1
rep1_AP <- read.table("HTSeq_out_rep1_AP.txt")
rep1_AP <- rep1_AP[grep("AT",rep1_AP$V1),]
rep1_AP$V1 <- as.character(rep1_AP$V1)
colnames(rep1_AP)[2] <- "AP.rep1.count"
mer <- merge(rep1_AP,genelength,by.x="V1",by.y="Geneid")
head(mer)
colnames(mer)[1] <- "gene"
mer$rpk_AP <- mer[,2]/mer[,3]
mer$AP.rep1 <- mer$rpk_AP*1e6/sum(mer[,4])
bulk_rice <- mer[,c("gene","AP.rep1")]
#The data of pre-protoplasting in replicate1
rep1_BP <- read.table("HTSeq_out_rep1_BP.txt")
rep1_BP <- rep1_BP[grep("AT",rep1_BP$V1),]
rep1_BP$V1 <- as.character(rep1_BP$V1)
colnames(rep1_BP)[2] <- "BP.rep1.count"
mer <- merge(rep1_BP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_BP <- mer[,2]/mer[,3]
mer$BP.rep1 <- mer$rpk_BP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","BP.rep1")],by="gene")
#The data of post-protoplasting in replicate2
rep2_AP <- read.table("HTSeq_out_rep2_AP.txt")
rep2_AP <- rep2_AP[grep("AT",rep2_AP$V1),]
rep2_AP$V1 <- as.character(rep2_AP$V1)
colnames(rep2_AP)[2] <- "AP.rep2.count"
mer <- merge(rep2_AP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_AP <- mer[,2]/mer[,3]
mer$AP.rep2 <- mer$rpk_AP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","AP.rep2")],by="gene")
#The data of pre-protoplasting in replicate2
rep2_BP <- read.table("HTSeq_out_rep2_BP.txt")
rep2_BP <- rep2_BP[grep("AT",rep2_BP$V1),]
rep2_BP$V1 <- as.character(rep2_BP$V1)
colnames(rep2_BP)[2] <- "BP.rep2.count"
mer <- merge(rep2_BP,genelength,by.x="V1",by.y="Geneid")
colnames(mer)[1] <- "gene"
mer$rpk_BP <- mer[,2]/mer[,3]
mer$BP.rep2 <- mer$rpk_BP*1e6/sum(mer[,4])
bulk_rice <- merge(bulk_rice,mer[,c("gene","BP.rep2")],by="gene")

bulk_rice$mean_AP <- rowMeans(bulk_rice[,c("AP.rep1","AP.rep2")])
bulk_rice$mean_BP <- rowMeans(bulk_rice[,c("BP.rep1","BP.rep2")])
#Identification of protoplasting-sensitive genes
bulk_rice$fc_rep1 <- bulk_rice$AP.rep1/bulk_rice$BP.rep1
head(bulk_rice)
bulk_rice$fc_rep2 <- bulk_rice$AP.rep2/bulk_rice$BP.rep2
head(bulk_rice)
ggplot(bulk_rice,aes(log2(fc_rep1),log2(fc_rep2)))+
  geom_point(size=.5)+my_theme
cutoff <- 3
bulk_rice$group<- "N"
bulk_rice$group[(log2(bulk_rice$fc_rep1)>cutoff&log2(bulk_rice$fc_rep2)>cutoff)|
                  (log2(bulk_rice$fc_rep1)< -cutoff&log2(bulk_rice$fc_rep2)< -cutoff)]<- "Y"

#2.Seurat####
#Reading data
all_integ <- as.data.frame(fread("all.csv"))
rownames(all) <- all$V1
all <- all[,-1]
#Excluding protoplasting responding genes
all_integ <- all_integ[!(rownames(all_integ) %in% bulk_rice[bulk_rice$group=="Y",]$gene),]
#Creating the Seurat object
sample <- separate(data = data.frame(cell=colnames(all_integ)), col = "cell", 
                   into = c("barcode", "sample"), sep = "-")$sample
mydata <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
mydata@meta.data$sample <- sample
future::plan("multiprocess", workers = 10)
options(future.globals.maxSize = 80 * 1024^3)
#SCTransform
rice.list <- SplitObject(mydata, split.by = "sample")
for (i in 1:length(rice.list)) {
  rice.list[[i]] <- SCTransform(rice.list[[i]], verbose = T)
}
#Integration
rice.features <- SelectIntegrationFeatures(object.list = rice.list,nfeatures = 3000)
rice.list <- PrepSCTIntegration(object.list = rice.list, anchor.features = rice.features, 
                                verbose = T)
reference_dataset <- which(names(rice.list) == "Ctrl1")
rice.anchors <- FindIntegrationAnchors(object.list = rice.list, normalization.method = "SCT", 
                                       anchor.features = rice.features, reference = reference_dataset)
rice.integrated <- IntegrateData(anchorset = rice.anchors, normalization.method = "SCT")
#Dimensionality reduction
rice.integrated <- RunPCA(object = rice.integrated, verbose = FALSE, npcs = 100)
rice.integrated <- RunUMAP(object = rice.integrated, dims = 1:10, 
                           min.dist = 0.05, n.neighbors = 5, seed.use = 100)
data <- as.data.frame(rice.integrated@reductions[["pca"]]@cell.embeddings)
umap <- umap::umap(data[,1:100],
                   n_neighbors=10,metric="pearson2",
                   min_dist=0.01,random_state=39)
plt <- as.data.frame(umap$layout)
plt$sample <- sample
plt$cell <- rownames(plt)
data <- as.matrix(plt[colnames(rice.integrated),c("V1","V2")])
rownames(data) <- plt[colnames(rice.integrated),]$cell
colnames(data) <- c("UMAP_1","UMAP_2")
rice.integrated@reductions$umap@cell.embeddings <- data
#Clustering
DefaultAssay(rice.integrated) <- "integrated"
rice.integrated <- FindNeighbors(rice.integrated, reduction = "pca", dims = 1:100)
rice.integrated <- FindClusters(rice.integrated, resolution = 0.75, n.start = 10)
plt$cluster <- rice.integrated@meta.data$seurat_clusters
plt$cell <- rownames(plt)
plt <- merge(plt,data.frame(c2=c(0:28),stringsAsFactors = F,
                            cluster=c(24,13,3,23,7,21,17,1,10,9,8,5,6,14,27,20,19,22,4,12,2,0,15,26,11,18,25,16,28)),by="cluster")
plt$cluster <- plt$c2
ggplot(plt,aes(x=V1,y=V2,color=cluster))+
  geom_point(alpha=0.5,size=.1)
#log-normalization
mydata_allgene <- CreateSeuratObject(counts = all_integ, project = "mydata_scRNAseq")
DefaultAssay(mydata_allgene) <- "RNA"
mydata_allgene <- NormalizeData(mydata_allgene, normalization.method = "LogNormalize",scale.factor = 1e6)
mydata_allgene$sample <- plt[colnames(mydata_allgene),]$sample
mydata_allgene <- FindVariableFeatures(mydata_allgene, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(mydata_allgene)
mydata_allgene <- ScaleData(mydata_allgene, features = all.genes)
mydata_allgene <- RunPCA(object = mydata_allgene, verbose = FALSE, npcs = 10)
mydata_allgene <- RunUMAP(object = mydata_allgene, dims = 1:10, 
                          min.dist = 0.05, n.neighbors = 5, seed.use = 100)
data <- as.matrix(plt[colnames(mydata_allgene),c("V1","V2")])
rownames(data) <- plt[colnames(mydata_allgene),]$cell
colnames(data) <- c("UMAP_1","UMAP_2")
mydata_allgene@reductions$umap@cell.embeddings <- data
mydata_allgene@meta.data$seurat_clusters <- plt[colnames(mydata_allgene),]$cluster