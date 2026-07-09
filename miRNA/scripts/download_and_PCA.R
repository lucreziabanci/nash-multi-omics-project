#!/bin/bash

# miRNA analysis: download data and PCA
# Dataset:GSE59492

library(GEOquery)
library(Biobase)
library(limma)
library(umap)
library(ggplot2)

# 1. Download data from NCBI
# load series and platform data from GEO
gset <- getGEO(
	"GSE59492", 
	GSEMatrix =TRUE, 
	getGPL=FALSE3
)
if (length(gset) > 1) idx <- grep("GPL16384", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]

keep <- pData(gset)$"immunotolerance group:ch1" %in% c("NASH-CH", "CTRL")
expr <- exprs(gset)[, keep]
group <- factor(pData(gset)$"immunotolerance group:ch1"[keep])
levels(group)

# 2. PCA
v <- apply(expr, 1, var, na.rm=T)
expr_mat <- expr[v > 0, ]

pca <- prcomp(t(expr_mat),scale=TRUE)

percentVar <- pca$sdev^2 / sum(pca$sdev^2)

xlab = paste0("PC1: ", round(percentVar[1]*100, 1), "% variance")
ylab = paste0("PC2: ", round(percentVar[2]*100, 1), "% variance")

pca_df <- data.frame(
	PC1 = pca$x[,1],
	PC2 = pca$x[,2],
	Condition = group
)
pca_df<-pca_df[complete.cases(pca_df), ]

x_min <- min(pca_df[,1])
x_max <- max(pca_df[,1])
y_min <- min(pca_df[,2])
y_max <- max(pca_df[,2])
max <- abs(round(max(x_max,y_max) + 10))
min <- abs(round(min(x_min,y_min) - 10))
range=c(-min,max)

png("plots/pca.png", width = 500, height = 500) #save the plot
plot.new() #open
par(bg = 'white')
condition_colors<-c(
	"NASH-CH"="blue",
	"CTRL"="red"
)

cols<-condition_colors[group]
plot(
	pca_df[,1], pca_df[,2],
        xlim = range,
        ylim = range,
        xlab = xlab,
        ylab = ylab,
	main = "GSE59492",
        pch = 19,
	col=cols
)

legend(
	"bottomleft",
	legend=names(condition_colors),
	col=condition_colors,
	pch=19,
	title="Condition"
)
dev.off()
