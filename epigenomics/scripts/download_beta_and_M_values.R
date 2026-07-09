#!/bin/bash

# Methylation array: data download, Beta-values and M-values
# Dataset:GSE48325

library(GEOquery)
library(limma)
library(ggplot2)
library(Biobase)
library(minfi)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)

# 1. Download data from NCBI

gset <- getGEO("GSE48325", GSEMatrix =TRUE, getGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL13534", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

keep<-c("GSM1174837","GSM1174838","GSM1174839","GSM1174840","GSM1174841","GSM1174842","GSM1174843","GSM1174844","GSM1174845","GSM1174847","GSM1174848","GSM1174856","GSM1174858","GSM1174865","GSM1174867","GSM1174869","GSM1174872","GSM1174873","GSM1174874","GSM1174875","GSM1174877","GSM1174880","GSM1174882","GSM1174889","GSM1174891","GSM1174893","GSM1174894","GSM1174901","GSM1174903","GSM1174910","GSM1174912","GSM1174913","GSM1174916")
gset <- gset_1[, sampleNames(gset_1)%in%keep]
matrix <- read.table("GSE48325/GSE48325_signal_intensities.txt", header=T, sep="\t", check.names=F)

#select only column NASH/CONTROL
keep_id <- c(1,2,3,4,5,6,7,8,9,12,13,21,23,30,32,34,37,38,39,40,42,45,47,54,56,58,59,66,68,75,77,78,81)
pattern<-paste0("^(", paste(keep_id, collapse="|"), ")\\.") 
cols_keep<-grep(
	pattern, 
	colnames(matrix), 
	value=T
)
matrix_filtered<-matrix[,c("ID-REF", cols_keep)]

# 2. Quality filter: p-value < 0.05
pval_cols<-grep(
	"Detection", 
	colnames(matrix_filtered), 
	value=T
)
matrix_pvalue<-matrix_filtered[, pval_cols] #create matrix only with Pval, then filter

filter<- apply(matrix_pvalue, 1, function(x) all(x <0.05)) #logic to apply the filter
matrix_qual <- matrix_filtered[filter, ] #filter
rownames(matrix_qual)<-matrix_qual$"ID-REF"

#3. beta and M-values
#create beta-values and assign the correct name to the column
A_cols <- grep("Signal_A", colnames(matrix_qual), value=T)
B_cols <- grep("Signal_B", colnames(matrix_qual), value=T)
beta<- matrix_qual[, B_cols] / (matrix_qual[,A_cols]+matrix_qual[, B_cols]+100)
colnames(beta)<-c("GSM1174837","GSM1174843","GSM1174838","GSM1174844","GSM1174839","GSM1174845","GSM1174840","GSM1174841","GSM1174842","GSM1174847","GSM1174848","GSM1174858","GSM1174869","GSM1174882","GSM1174893","GSM1174894","GSM1174901","GSM1174916","GSM1174856","GSM1174865","GSM1174867","GSM1174877","GSM1174872","GSM1174873","GSM1174874","GSM1174880","GSM1174875","GSM1174891","GSM1174889","GSM1174910","GSM1174903","GSM1174912","GSM1174913")

#avoid 0 and 1 value because in the M values formula if B=0 M=-infinite and if B=1 M=+infinite
beta[beta == 0] <- 1e-6
beta[beta == 1] <- 1 - 1e-6

#create M values
Mval <- log2(beta / (1 - beta))
