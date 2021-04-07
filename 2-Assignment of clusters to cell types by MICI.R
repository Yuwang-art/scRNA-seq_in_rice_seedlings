library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(parallel)
library(patchwork)
#########################################################
#                                                       #
#                          MICI                         #
#                                                       #
#########################################################
#Leaf####
marker <- as.data.frame(fread('Table S5',header = T))
marker$`Gene ID` <- gsub("_","-",marker$`Gene ID`)
rownames(marker) <- 1:dim(marker)[1]
colnames(marker)[1] <- "gene"
marker1 <-  marker[marker$Location %in% c("Both","Leaf"),]
dim(marker1)
marker2 <- marker1 %>% separate_rows(Leaf,sep = ",")#The column "Leaf" is the tissue of leaf
dim(marker2)
sam <- c("S_A1","S_A2")
expr_data <- as.data.frame(mydata_allgene@assays$RNA@data[
  intersect(rownames(mydata_allgene@assays$RNA@data),unique(marker1$gene)),plt[plt$sample %in% sam,]$cell])
#1.Weight estimation for each marker genes####
weight_data <- data.frame(stringsAsFactors = F)
for (i in 0:28){
  for (j in 1:length(marker1$gene)) {
    weight_data[j,i+1] <- log2(rowMeans(expm1(expr_data[marker1$gene[j],plt[plt$cluster==i&plt$sample %in% sam,]$cell]))+1)
  }
}
rownames(weight_data) <- marker1$gene
weight_data$weight <- apply(weight_data,1,var)
weight_data$gene <- rownames(weight_data)
weight_data <- merge(weight_data,marker1[,c("gene","Tissue","Leaf")],by="gene")
#2.Estimation of MICI for each cell and for each cell type####
expr_scale_data <- as.data.frame(t(scale(t(expr_data))))
expr_scale_data <- merge(marker2,expr_scale_data,by.x='gene',by.y="row.names")
expr_scale_data <- merge(weight_data[,c("gene","weight")],expr_scale_data,by="gene")
expr_scale_data[,plt[plt$sample %in% sam,]$cell] <- expr_scale_data[,plt[plt$sample %in% sam,]$cell]*expr_scale_data$weight
fun <- function(x) {
  tmp <- expr_scale_data[expr_scale_data$Leaf==unique(expr_scale_data$Leaf)[x],plt[plt$sample %in% sam,]$cell]
  tmp1 <- as.data.frame(t(as.data.frame(colSums(tmp))))
  rownames(tmp1) <- unique(dat$Leaf)[x]
  return(tmp1)
}
MICI_out <- do.call('rbind',parallel::mclapply(1:length(unique(expr_scale_data$Leaf)),
                                               function(x){fun(x)},mc.cores = length(unique(expr_scale_data$Leaf))))
MICI_out$cell_type <- as.character(rownames(MICI_out))
MICI_result <-apply(MICI_out[,-dim(MICI_out)[2]], 2, function(x){MICI_out[,dim(MICI_out)[2]][which(x==max(x))]}) %>% as.data.frame()
colnames(MICI_result) <- c("identified_ct")
MICI_result <- merge(plt,MICI_result,by.x="cell",by.y="row.names")
MICI_result$identified_ct <- as.character(MICI_result$identified_ct)
tmp <- as.data.frame(t(MICI_out[,-dim(MICI_out)[2]]))
unct <- rownames(tmp[rowMax(as.matrix(tmp))<2,])
MICI_result[MICI_result$cell %in% unct,]$identified_ct <- "Unknown"
#3.Assignment of clusters to cell type####
dat1 <- as.data.frame(table(MICI_result$identified_ct,MICI_result$cluster))
colnames(dat1) <- c("cell_type","cluster","number_of_cell")
dat2 <- as.data.frame(table(MICI_result$cluster))
colnames(dat2) <- c("cluster","number_of_cluster")
dat3 <- merge(dat1,dat2,by="cluster")
dat3$ratio <- signif(dat3$number_of_cell/dat3$number_of_cluster,3)
dat4 <- dat3 %>% group_by(cluster) %>% top_n(1,ratio) %>% as.data.frame()
dat5 <- reshape2::dcast(dat3,cluster~cell_type,value.var = "ratio")
cluster_assign <- merge(dat4[,c("cluster","cell_type")],dat5,by="cluster")
#Root####
marker <- as.data.frame(fread('Table S5',header = T))
marker$`Gene ID` <- gsub("_","-",marker$`Gene ID`)
rownames(marker) <- 1:dim(marker)[1]
colnames(marker)[1] <- "gene"
marker1 <-  marker[marker$Location %in% c("Both","Root"),]
dim(marker1)
marker2 <- marker1 %>% separate_rows(Leaf,sep = ",")#The column "Root" is the tissue of root
dim(marker2)
sam <- c("S_R1","S_R2")
clu <- c(0:22,25,26:28)  
expr_data <- as.data.frame(mydata_allgene@assays$RNA@data[
  intersect(rownames(mydata_allgene@assays$RNA@data),unique(marker1$gene)),plt[plt$sample %in% sam&plt$cluster %in% clu,]$cell])
#1.Weight estimation for each marker genes####
weight_data <- data.frame(stringsAsFactors = F)
for (i in 1:length(clu)){
  for (j in 1:length(marker1$gene)) {
    tmp[j,i] <- log2(rowMeans(expm1(dat[marker1$gene[j],plt[plt$cluster==clu[i]&plt$sample %in% sam,]$cell]))+1)
  }
}
rownames(weight_data) <- marker1$gene
weight_data$weight <- apply(weight_data,1,var)
weight_data$gene <- rownames(weight_data)
weight_data <- merge(weight_data,marker1[,c("gene","Tissue","Leaf")],by="gene")
#2.Estimation of MICI for each cell and for each cell type####
expr_scale_data <- as.data.frame(t(scale(t(expr_data))))
expr_scale_data <- merge(marker2,expr_scale_data,by.x='gene',by.y="row.names")
expr_scale_data <- merge(weight_data[,c("gene","weight")],expr_scale_data,by="gene")
expr_scale_data[,plt[plt$sample %in% sam,]$cell] <- expr_scale_data[,plt[plt$sample %in% sam,]$cell]*expr_scale_data$weight
fun <- function(x) {
  tmp <- expr_scale_data[expr_scale_data$Leaf==unique(expr_scale_data$Leaf)[x],plt[plt$sample %in% sam,]$cell]
  tmp1 <- as.data.frame(t(as.data.frame(colSums(tmp))))
  rownames(tmp1) <- unique(dat$Leaf)[x]
  return(tmp1)
}
MICI_out <- do.call('rbind',parallel::mclapply(1:length(unique(expr_scale_data$Leaf)),
                                               function(x){fun(x)},mc.cores = length(unique(expr_scale_data$Leaf))))
MICI_out$cell_type <- as.character(rownames(MICI_out))
MICI_result <-apply(MICI_out[,-dim(MICI_out)[2]], 2, function(x){MICI_out[,dim(MICI_out)[2]][which(x==max(x))]}) %>% as.data.frame()
colnames(MICI_result) <- c("identified_ct")
MICI_result <- merge(plt,MICI_result,by.x="cell",by.y="row.names")
MICI_result$identified_ct <- as.character(MICI_result$identified_ct)
tmp <- as.data.frame(t(MICI_out[,-dim(MICI_out)[2]]))
unct <- rownames(tmp[rowMax(as.matrix(tmp))<2,])
MICI_result[MICI_result$cell %in% unct,]$identified_ct <- "Unknown"
#3.Assignment of clusters to cell type####
dat1 <- as.data.frame(table(MICI_result$identified_ct,MICI_result$cluster))
colnames(dat1) <- c("cell_type","cluster","number_of_cell")
dat2 <- as.data.frame(table(MICI_result$cluster))
colnames(dat2) <- c("cluster","number_of_cluster")
dat3 <- merge(dat1,dat2,by="cluster")
dat3$ratio <- signif(dat3$number_of_cell/dat3$number_of_cluster,3)
dat4 <- dat3 %>% group_by(cluster) %>% top_n(1,ratio) %>% as.data.frame()
dat5 <- reshape2::dcast(dat3,cluster~cell_type,value.var = "ratio")
cluster_assign <- merge(dat4[,c("cluster","cell_type")],dat5,by="cluster")