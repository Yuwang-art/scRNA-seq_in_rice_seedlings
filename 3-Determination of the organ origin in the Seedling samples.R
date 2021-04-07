library(data.table)
library(tidyverse)
library(Seurat)
library(reshape2)
library(Matrix)
library(dplyr)
library(parallel)
###########################################################################
#                                                                         #
#        Determination of the organ origin in the Seedling samples        #
#                                                                         #
###########################################################################
#1.Selection of leaf-biased genes and root-biased genes####
#batch1
Idents(mydata_allgene) <- "sample"
ae.ro1 <- subset(mydata_allgene,idents = c("S_A1","S_R1"))
ae.ro1@meta.data$Sample <- substr(ae.ro1@meta.data$sample,1,3)
ae.ro1$cluster.stress <- paste(ae.ro1$seurat_clusters, ae.ro1$Sample, sep = "_")
Idents(ae.ro1) <- "cluster.stress"
unique(ae.ro1$cluster.stress)
DefaultAssay(ae.ro1) <- "RNA"
ae.ro1$cluster.stress <- as.character(ae.ro1$cluster.stress)
fun <- function(x){
  response.ae.ro1 <- FindMarkers(ae.ro1, slot="data",test.use = "wilcox",
                                 ident.1 = as.character(paste0(x,"_S_A")), 
                                 ident.2 = as.character(paste0(x,"_S_R")), 
                                 verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro1$cluster <- x
  response.ae.ro1$gene <- rownames(response.ae.ro1)
  rownames(response.ae.ro1) <- NULL
  return(response.ae.ro1)
}
clu <- c(0:22,26:28)  
deg.ae.ro1 <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
deg.ae.ro1$batch <- "batch1"
#batch2
Idents(mydata_allgene) <- "sample"
ae.ro2 <- subset(mydata_allgene,idents = c("S_A2","S_R2"))
ae.ro2@meta.data$Sample <- substr(ae.ro2@meta.data$sample,1,3)
ae.ro2$cluster.stress <- paste(ae.ro2$seurat_clusters, ae.ro2$Sample, sep = "_")
Idents(ae.ro2) <- "cluster.stress"
DefaultAssay(ae.ro2) <- "RNA"
ae.ro2$cluster.stress <- as.character(ae.ro2$cluster.stress)
fun <- function(x){
  response.ae.ro2 <- FindMarkers(ae.ro2, slot="data",test.use = "wilcox",
                                 ident.1 = as.character(paste0(x,"_S_A")), 
                                 ident.2 = as.character(paste0(x,"_S_R")), 
                                 verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro2$cluster <- x
  response.ae.ro2$gene <- rownames(response.ae.ro2)
  rownames(response.ae.ro2) <- NULL
  return(response.ae.ro2)
}
clu <- c(0:22,25:28)  
deg.ae.ro2 <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
table(deg.ae.ro2$cluster)
deg.ae.ro2$batch <- "batch2"
#batch1 without cluster25
deg.ae.ro1 <- rbind(deg.ae.ro1,deg.ae.ro2[deg.ae.ro2$cluster==25,])
deg.ae.ro1[deg.ae.ro1$cluster==25,]$batch <- "batch1"
deg.ae.ro <- rbind(deg.ae.ro1,deg.ae.ro2)
colnames(deg.ae.ro)[3:4] <- c("pct.aerial","pct.root")
deg.ae.ro$direction <- "ae"
deg.ae.ro[deg.ae.ro$avg_logFC <0,]$direction <- "ro"
deg.ae.ro$gene_cluster <- paste0(deg.ae.ro$gene,"-",deg.ae.ro$cluster)
deg.ae.ro$log2fc_avg <- log2(exp(deg.ae.ro$avg_logFC))
prop <- 0.1
deg.ae.ro[deg.ae.ro$avg_logFC>0&deg.ae.ro$pct.root<0.1,]$pct.root <- prop
deg.ae.ro[deg.ae.ro$avg_logFC<0&deg.ae.ro$pct.aerial<0.1,]$pct.aerial <- prop
deg.ae.ro$log2fc_pct <- log2(deg.ae.ro$pct.aerial/deg.ae.ro$pct.root)
deg.ae.ro$score <- deg.ae.ro$log2fc_avg+deg.ae.ro$log2fc_pct
head(deg.ae.ro)
#2.Selection of genes that unchaged in stress condition####
clu <- c(0:22,25:28)  
deg.ae.ro$Ct <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$Ct <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="Ct",]$cell]))+1)
  }
}
deg.ae.ro$HS <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$HS <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="HS",]$cell]))+1)
  }
}
deg.ae.ro$ID <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$ID <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="ID",]$cell]))+1)
  }
}
deg.ae.ro$LN <- 0
for (i in clu){
  for(j in c("batch1","batch2")){
    deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$LN <- log2(rowMeans(expm1(mydata_allgene@assays$RNA@data[deg.ae.ro[deg.ae.ro$cluster==i&deg.ae.ro$batch==j,]$gene,plt[plt$cluster==i&plt$Sample=="LN",]$cell]))+1)
  }
}

deg.twobatch <- reshape2::dcast(deg.ae.ro[deg.ae.ro$p_val<0.01,],gene_cluster+cluster+gene+Ct+HS+ID+LN+direction~batch,value.var = "score")
deg.twobatch <- na.omit(deg.twobatch)
deg.twobatch$mean_score <- rowMeans(deg.twobatch[,c("batch1","batch2")])
deg.twobatch$is_stressresponse <-"Y"
deg.twobatch[abs(deg.twobatch$HS-deg.twobatch$Ct)<1&
               abs(deg.twobatch$ID-deg.twobatch$Ct)<1&
               abs(deg.twobatch$LN-deg.twobatch$Ct)<1, ]$is_stressresponse <-"N"
table(deg.twobatch$is_stressresponse)

num <- 3
marker.ae <- deg.twobatch[deg.twobatch$direction=="ae"&deg.twobatch$is_stressresponse=="N",] %>% 
  group_by(cluster) %>% top_n(num,mean_score) %>% as.data.frame()
length(marker.ae$gene)
length(unique(marker.ae$gene))
marker.ro <- deg.twobatch[deg.twobatch$direction=="ro"&deg.twobatch$is_stressresponse=="N",] %>% 
  group_by(cluster) %>% top_n(num,abs(mean_score)) %>% as.data.frame()
length(marker.ro$gene)
length(unique(marker.ro$gene))
deg.twobatch$is_deg <- "N"
deg.twobatch[deg.twobatch$direction=="ae"&deg.twobatch$gene_cluster %in% marker.ae$gene_cluster,]$is_deg <- "Y"
deg.twobatch[deg.twobatch$direction=="ro"&deg.twobatch$gene_cluster %in% marker.ro$gene_cluster,]$is_deg <- "Y"
table(deg.twobatch$is_deg,deg.twobatch$cluster,deg.twobatch$direction)

deg.twobatch <- deg.twobatch[order(deg.twobatch$is_deg),]
#3.wilcox test####
Idents(mydata_allgene) <- "sample"
ae.ro <- subset(mydata_allgene,idents = c("S_A1","S_A2","S_R1","S_R2"))
ae.ro
ae.ro@meta.data$Sample <- substr(ae.ro@meta.data$sample,1,3)
ae.ro$cluster.stress <- paste(ae.ro$seurat_clusters, ae.ro$Sample, sep = "_")
Idents(ae.ro) <- "cluster.stress"
unique(ae.ro$cluster.stress)
DimPlot(ae.ro,reduction = 'umap',split.by = 'sample',label = T,group.by = 'seurat_clusters')
DefaultAssay(ae.ro) <- "RNA"
ae.ro$cluster.stress <- as.character(ae.ro$cluster.stress)
fun <- function(x){
  response.ae.ro <- FindMarkers(ae.ro, slot="data",test.use = "wilcox",
                                ident.1 = as.character(paste0(x,"_S_A")), 
                                ident.2 = as.character(paste0(x,"_S_R")), 
                                verbose = T, logfc.threshold = 0, pseudocount.use=1)
  response.ae.ro$cluster <- x
  response.ae.ro$gene <- rownames(response.ae.ro)
  rownames(response.ae.ro) <- NULL
  return(response.ae.ro)
}
clu <- c(0:25,27:28)
deg.ae.ro.mix2rep <- do.call('rbind',parallel::mclapply(clu,function(x){fun(x)},mc.cores = 28))
table(deg.ae.ro.mix2rep$cluster)
head(deg.ae.ro.mix2rep)
dim(deg.ae.ro.mix2rep)

deg.ae.ro.mix2rep$log2fc_avg <- log2(exp(deg.ae.ro.mix2rep$avg_logFC))
deg.ae.ro.mix2rep$direction <- "ae"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$avg_logFC <0,]$direction <- "ro"
deg.ae.ro.mix2rep$gene_cluster <- paste0(deg.ae.ro.mix2rep$gene,"-",deg.ae.ro.mix2rep$cluster)
deg.ae.ro.mix2rep$is_deg <- "N"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$direction=="ae"&deg.ae.ro.mix2rep$gene_cluster %in% marker.ae$gene_cluster,]$is_deg <- "Y"
deg.ae.ro.mix2rep[deg.ae.ro.mix2rep$direction=="ro"&deg.ae.ro.mix2rep$gene_cluster %in% marker.ro$gene_cluster,]$is_deg <- "Y"
deg.ae.ro.mix2rep <- deg.ae.ro.mix2rep[order(deg.ae.ro.mix2rep$is_deg),]

#4.Estimation the average scaled expression level of the three leaf-biased (root-biased) genes as the index of leaf (root) for the cell####
plt$sample2 <- plt$sample
plt[plt$sample %in% c("S_A1","S_A2","S_R1","S_R2"),]$sample2 <- "Small"

plt.ae.ro <- data.frame(stringsAsFactors = F)
deg <- deg.ae.ro.mix2rep
score <- function(i) {
  tmp <- data.frame(stringsAsFactors = F)
  for (j in sample2){
    plt.tmp <- plt[plt$cluster==i&plt$sample2==j,]
    gene.ae <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC>0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ae,plt.tmp$cell]))
    if (length(gene.ae)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"]) 
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$aerial <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$aerial <- scale(dat)
    } else {
      plt.tmp$aerial <- -5
    }
    }
    
    gene.ro <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC<0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ro,plt.tmp$cell]))
    if (length(gene.ro)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"]) 
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$root <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$root <- scale(dat)
    } else {
      plt.tmp$root <- -5
    }
    }
    plt.tmp$aerial <- plt.tmp$aerial-min(plt.tmp$aerial)
    plt.tmp$root <- plt.tmp$root-min(plt.tmp$root)
    
    tmp <- rbind(tmp,plt.tmp)
  }
  return(tmp)
}
clu1 <- c(0:13,15:22,25,27:28)
plt.ae.ro <- do.call('rbind',parallel::mclapply(clu1,function(i){score(i)},mc.cores = 28))
plt.ae.ro[plt.ae.ro$sample %in% c("S_A1","S_A2"),]$Sample <- "Aerial"
plt.ae.ro[plt.ae.ro$sample %in% c("S_R1","S_R2"),]$Sample <- "Root"
plt.score.1 <- plt.ae.ro
plt.score.1$tissue <- "N"
plt.score.1[plt.score.1$root< plt.score.1$aerial,]$tissue <- "Aerial"
plt.score.1[plt.score.1$root> plt.score.1$aerial,]$tissue <- "Root"

plt.ae.ro <- data.frame(stringsAsFactors = F)
deg <- deg.ae.ro.mix2rep
score <- function(i) {
  tmp <- data.frame(stringsAsFactors = F)
  sample3 <- c("Ctrl1","Ctrl2","Ctrl3","Ctrl4","ID1","ID2","ID3","ID4","LN1","LN2","LN3","LN4","Small")
  for (j in sample3){
    plt.tmp <- plt[plt$cluster==i&plt$sample2==j,]
    gene.ae <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC>0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ae,plt.tmp$cell]))
    if (length(gene.ae)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"]) 
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$aerial <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$aerial <- scale(dat)
    } else {
      plt.tmp$aerial <- -5
    }
    }
    
    gene.ro <- deg[deg$cluster==i&deg$is_deg=="Y"&deg$avg_logFC<0,]$gene
    dat <- expm1(as.data.frame(mydata_allgene@assays$RNA@data[gene.ro,plt.tmp$cell]))
    if (length(gene.ro)>1){
      if (dim(dat[rowSums(dat)>0,])[1]>0){
        dat <- as.data.frame(t(scale(t(dat))))
        dat[dat=="NaN"] <- min(dat[dat!="NaN"]) 
      } else {
        dat[dat==0] <- -5
      }
      plt.tmp$root <- as.numeric(colMeans(dat[,plt.tmp$cell]))
    } else { if (length(dat[dat>0])>0){
      plt.tmp$root <- scale(dat)
    } else {
      plt.tmp$root <- -5
    }
    }
    plt.tmp$aerial <- plt.tmp$aerial-min(plt.tmp$aerial)
    plt.tmp$root <- plt.tmp$root-min(plt.tmp$root)
    tmp <- rbind(tmp,plt.tmp)
  }
  return(tmp)
}
clu2 <- c(14,26)
plt.ae.ro <- do.call('rbind',parallel::mclapply(clu2,function(i){score(i)},mc.cores = 8))
plt.ae.ro[plt.ae.ro$sample %in% c("S_A1","S_A2"),]$Sample <- "Aerial"
plt.ae.ro[plt.ae.ro$sample %in% c("S_R1","S_R2"),]$Sample <- "Root"
plt.score.2 <- plt.ae.ro
plt.score.2$tissue <- "N"
plt.score.2[plt.score.2$root< plt.score.2$aerial,]$tissue <- "Aerial"
plt.score.2[plt.score.2$root> plt.score.2$aerial,]$tissue <- "Root"

plt.score <- rbind(plt.score.1,plt.score.2)
plt.score$tissue <- "N"
plt.score[plt.score$root< plt.score$aerial,]$tissue <- "Aerial"
plt.score[plt.score$root> plt.score$aerial,]$tissue <- "Root"


