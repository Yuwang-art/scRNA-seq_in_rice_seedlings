library(data.table)
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)
library(tidyverse)
library(monocle)
library(parallel)
library(ggpubr)
#################################################################################
#                                                                               #
#          Reconstruction of developmental trajectory by monocle                #
#                                                                               #
#################################################################################
#Data input
data_input <- as.matrix(rice.integrated@assays[["integrated"]]@scale.data)
cell_input <- plt.score[plt.score$cluster %in% c(0,1,17:22)&plt.score$tissue=="Aerial"&
                plt.score$Sample!="Aerial"&plt.score$Sample!="Root",]
cell_input$random <- sample(1:dim(cell_input)[1],dim(cell_input)[1])
cell_input <- cell_input[order(cell_input$random),]
data_input <- data_input[,colnames(data_input) %in% cell_input$cell[1:10000]]
gene_annotation <- data.frame(row.names = rownames(data_input),gene_short_name =rownames(data_input))
gene_annotation$gene_short_name <- gsub("LOC-","",gene_annotation$gene_short_name)
raw_select <- as.matrix(data_input)
sample_sheet <- data.frame(row.names = colnames(raw_select),Library= rep("mes",dim(raw_select)[2]))
fd <- new("AnnotatedDataFrame", data = gene_annotation)
pd <- new("AnnotatedDataFrame", data = sample_sheet)
mono <- newCellDataSet(as(raw_select, "sparseMatrix"),phenoData = pd,featureData = fd)
mono <- estimateSizeFactors(mono)
mono <- detectGenes(mono, min_expr = 0.1)
#Dimensional reduction
ordering_genes <- rownames(df)
mono <- setOrderingFilter(mono, ordering_genes)
mono <- reduceDimension(mono, reduction_method = "DDRTree",norm_method="none",pseudo_expr=0)
#Calculation of pseudotime
mono <- orderCells(mono)
plot_cell_trajectory(mono, color_by="Pseudotime")
plot_cell_trajectory(mono, color_by="State")
mono <- orderCells(mono,root_state = 2)
plot_cell_trajectory(mono, color_by="Pseudotime")
