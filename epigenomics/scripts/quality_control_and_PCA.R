#!/bin/bash

# Methylation array quality control and PCA
# Dataset:GSE48325

# 1. Quality control
## Mean detection p-value
groups <- pData(gset)$group  
pal <- c("Control"="#2c7bb6", "NASH"="#d7191c")
col<-pal[as.factor(group)]
par(mfrow=c(1,2))

png(file="plots/Mean_detection_p-value.png", 2500, 2500, res = 500)
barplot(
	colMeans(matrix_pvalue), 
	col=col, 
	las=2, 
	cex.names=0.6, 
	cex.axis=0.7, 
	ylim=c(0, 0.002), 
	ylab=""
)
title(
	ylab="Mean detection p-values", 
	line=3.2
)
abline(
	h=0.05,
	col="red"
)
legend(
	"topleft", 
	legend=names(pal), 
	fill=pal, 
	bg="white"
)
dev.off()

# Density
#beta values
normalBeta <- beta[, group=="NORMAL"]
nashBeta <- beta[, group=="NASH"]

png(file="plots/Beta_density.png", 2500, 2500, res = 500)
plot(
	density(
		rowMeans(normalBeta)
	), 
	main="Density of Beta Values", 
	col="#0067E6",
	lwd=2.5
)
lines(
	density(
		rowMeans(nashBeta)
	), 
	main="Density of Beta Values", 
	col="#E50068", 
	lwd=2.5
)
legend(
	'topright', 
	legend=c("Control", "NASH"), 
	fill=c("#0067E6","#E50068"), 
	cex=1
)
dev.off()

#M values
normalM<-Mval[, group=="NORMAL"]
nashM<-Mval[, group=="NASH"]

png(file="plots/M_density.png", 2500, 2500, res = 500)
plot(
	density(
		rowMeans(normalM)
	), 
	main="Density of M Values", 
	col="#0067E6",
	lwd=2.5
)
lines(
	density(
		rowMeans(nashM)
	), 
	main="Density of M Values", 
	col="#E50068", 
	lwd=2.5
)
legend(
	'topright', 
	legend=c("Control", "NASH"), 
	fill=c("#0067E6","#E50068"), 
	cex=1
)
dev.off()

# 2. PCA
cols<-ifelse(group == "NORMAL", "#0067E6", "#E50068")
Mval_var <- apply(Mval, 1, var)
Mval_pca <- Mval[Mval_var > 0, ]

pca <- prcomp(t(Mval_pca))

percentVar <- (pca$sdev^2 / sum(pca$sdev^2)) * 100

xlab <- paste0("PC1: ", round(percentVar[1], 1), "% variance")
ylab <- paste0("PC2: ", round(percentVar[2], 1), "% variance")

png(file="plots/PCA.png", 2500, 2500, res = 500)
plot(
	pca$x[,1], 
	pca$x[,2], 
	col=cols, 
	pch=19, 
	xlab=xlab, 
	ylab=ylab,  
	main="GSE48325"
)
legend(
	"bottomright", 
	legend=c("Control", "NASH"), 
	col=c("#0067E6", "#E50068"), 
	pch=19, 
	cex=1.0
)
dev.off()
