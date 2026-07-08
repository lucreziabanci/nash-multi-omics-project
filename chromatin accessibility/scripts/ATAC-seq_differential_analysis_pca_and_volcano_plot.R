#!/bin/bash

# ATAC-seq differential analysis, pca and volcano plot
# Dataset:PRJNA725028

# 1. Differential analysis
# move from directory try to directory peaks
for i in try/SRR*/*.broadPeak; do 
	mv $i peak; 
done

## import in R
library(GenomicRanges)
library(rtracklayer)
library(ChIPpeakAnno) # peak annotation
library(EnsDb.Hsapiens.v75)
library(csaw)
library(DESeq2)
library(ggplot2)
path_to_peaks <- "peak"

n1_peak <- import(paste0(path_to_peaks, "SRR14321", 495, "_peaks.broadPeak"))
n2_peak <- import(paste0(path_to_peaks, "SRR14321", 497, "_peaks.broadPeak"))
n3_peak <- import(paste0(path_to_peaks, "SRR14321", 498, "_peaks.broadPeak"))
n4_peak <- import(paste0(path_to_peaks, "SRR14321", 499, "_peaks.broadPeak"))
n5_peak <- import(paste0(path_to_peaks, "SRR14321", 500, "_peaks.broadPeak"))
n6_peak <- import(paste0(path_to_peaks, "SRR14321", 501, "_peaks.broadPeak"))
n7_peak <- import(paste0(path_to_peaks, "SRR14321", 502, "_peaks.broadPeak"))
n8_peak <- import(paste0(path_to_peaks, "SRR14321", 503, "_peaks.broadPeak"))
n9_peak <- import(paste0(path_to_peaks, "SRR14321", 504, "_peaks.broadPeak"))
n10_peak <- import(paste0(path_to_peaks, "SRR14321", 505, "_peaks.broadPeak"))
n11_peak <- import(paste0(path_to_peaks, "SRR14321", 506, "_peaks.broadPeak"))
n12_peak <- import(paste0(path_to_peaks, "SRR14321", 508, "_peaks.broadPeak"))
n13_peak <- import(paste0(path_to_peaks, "SRR14321", 509, "_peaks.broadPeak"))
n14_peak <- import(paste0(path_to_peaks, "SRR14321", 510, "_peaks.broadPeak"))
n15_peak <- import(paste0(path_to_peaks, "SRR14321", 511, "_peaks.broadPeak"))
n16_peak <- import(paste0(path_to_peaks, "SRR14321", 512, "_peaks.broadPeak"))
n17_peak <- import(paste0(path_to_peaks, "SRR14321", 513, "_peaks.broadPeak"))
n18_peak <- import(paste0(path_to_peaks, "SRR14321", 514, "_peaks.broadPeak"))
n19_peak <- import(paste0(path_to_peaks, "SRR14321", 536, "_peaks.broadPeak"))
n20_peak <- import(paste0(path_to_peaks, "SRR14321", 537, "_peaks.broadPeak"))
n21_peak <- import(paste0(path_to_peaks, "SRR14321", 538, "_peaks.broadPeak"))

c1_peak <- import(paste0(path_to_peaks, "SRR14321", 507, "_peaks.broadPeak"))
c2_peak <- import(paste0(path_to_peaks, "SRR14321", 518, "_peaks.broadPeak"))
c3_peak <- import(paste0(path_to_peaks, "SRR14321", 529, "_peaks.broadPeak"))
c4_peak <- import(paste0(path_to_peaks, "SRR14321", 530, "_peaks.broadPeak"))

## intersect
nu1_peak <- union(n1_peak, n2_peak)
nu2_peak <- union (nu1_peak, n3_peak)
nu3_peak <- union (nu2_peak, n4_peak)
nu4_peak <- union (nu3_peak, n5_peak)
nu5_peak <- union (nu4_peak, n6_peak)
nu6_peak <- union (nu5_peak, n7_peak)
nu7_peak <- union (nu6_peak, n8_peak)
nu8_peak <- union (nu7_peak, n9_peak)
nu9_peak <- union (nu8_peak, n10_peak)
nu10_peak <- union (nu9_peak, n11_peak)
nu11_peak <- union (nu10_peak, n12_peak)
nu12_peak <- union (nu11_peak, n13_peak)
nu13_peak <- union (nu12_peak, n13_peak)
nu14_peak <- union (nu13_peak, n14_peak)
nu15_peak <- union (nu14_peak, n15_peak)
nu16_peak <- union (nu15_peak, n16_peak)
nu17_peak <- union (nu16_peak, n17_peak)
nu18_peak <- union (nu17_peak, n18_peak)
nu19_peak <- union (nu18_peak, n19_peak)
nu20_peak <- union (nu19_peak, n20_peak)
n_peak <- union (nu20_peak, n21_peak)

cu1_peak <- union(c1_peak, c2_peak)
cu2_peak <- union(cu1_peak, c3_peak)
c_peak <- union (cu2_peak, c4_peak)

intersect_peaks <- intersect(n_peak, c_peak)
intersect_peaks

## Annotation
annoData <- toGRanges(EnsDb.Hsapiens.v75)
library(GenomeInfoDb)
annoData <- keepSeqlevels(annoData, 
                          c(paste0("chr", c(1:22, "X", "Y")), "chrMT"), 
                          pruning.mode = "coarse")
seqlevels(annoData) <- gsub("chrMT", "chrM", seqlevels(annoData))

annotated_peaks <- annoPeaks(
	intersect_peaks, 
	annoData, 
	bindingRegion=c(-1000, 1000)
)

## Preparation for counting peaks
path_to_bam <- "mapped"
bamfiles <- paste0(path_to_bam, list.files(path=path_to_bam, pattern="*_shift.bam$"))

# Define read parameters
chr <- paste0("chr", c(1:22, "X", "Y")) # extract reads only from standard chromosomes

param <- readParam(max.frag=250, pe="both", restrict=chr)
bpparam <- BiocParallel::MulticoreParam(workers=16)

#order the controls and nash
controls <- grep("SRR14321507|SRR14321518|SRR14321529|SRR14321530", bamfiles)
nash <- setdiff(1:25, controls)
bamfiles_ordered <- bamfiles[c(nash, controls)]
basename(bamfiles_ordered)

peak_counts <- regionCounts(
	bamfiles_ordered, 
	annotated_peaks, 
	param=param, 
	BPPARAM=bpparam
)

## Count matrix for differential analysis
countdata <- assays(peak_counts)$counts
dim(countdata)

sel <- match(
	rownames(countdata), 
	annotated_peaks$peak
)
rownames(countdata) <- annotated_peaks$gene_name[sel]

colnames(countdata) <- c(paste("NASH", 1:21, sep=""), paste("control", 1:4 , sep=""))
countdata <- assays(peak_counts)$counts

condition <- c(rep("NASH", 21), rep("control", 4))
samples <- data.frame(
	condition, 
	row.names = colnames(countdata)
)
samples$condition <- factor(samples$condition)

# Differential analysis
# Create a DESeqDataSet object
design <- model.matrix(~0+samples$condition)

dds <- DESeqDataSetFromMatrix(countData=countdata,
                              colData=samples, design=design)

# This transforms data to regularize the variance and allow PCA
vst_d <- vst(dds) # variance stabilizing transformation
vst_mat<-assay(vst_d)
v <- apply(vst_mat, 1, var)
vst_mat <- vst_mat[v > 0, ]
pca <- prcomp(
	t(vst_mat), 
	scale= TRUE
)
percentVar <- pca$sdev^2 / sum(pca$sdev^2)
pca_df <- data.frame(
	PC1 = pca$x[,1],
	PC2 = pca$x[,2]
)
pca_df$condition <- samples[rownames(pca_df), "condition"]

xlab = paste0("PC1: ", round(percentVar[1]*100, 1), "% variance")
ylab = paste0("PC2: ", round(percentVar[2]*100, 1), "% variance")

x_min <- min(pca_df[,1])
x_max <- max(pca_df[,1])
y_min <- min(pca_df[,2])
y_max <- max(pca_df[,2])
max <- abs(round(max(x_max,y_max) + 10))
min <- abs(round(min(x_min,y_min) - 10))
range=c(-min, max)

condition_colors<-c(
	"NASH"="blue",
	"control"="red"
)
cols<-condition_colors[as.character(samples$condition)]

plot.new() #open
par(bg = 'white') #white background

#create the plot
png(file="plots/PCA.png", 2500, 2500, res = 500)
plot(
	pca_df[,1], pca_df[,2],
        xlim = range,
        ylim = range,
        xlab = xlab,
        ylab = ylab,
	main = "PRJNA725028",
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

# select rows (genes) with sufficient information
keep <- rowSums(counts(dds)) >= 10
dds <- dds[keep, ]

# run the analysis
dds_d <- DESeq(dds)
resultsNames(dds_d)

res <- data.frame(
	results(
		dds_d, 
		contrast=list("samples.conditionNASH", "samples.conditioncontrol")), 
		lfcThreshold=1, 
		alpha=0.05
	)
)
summary(res)
res <- res[complete.cases(res), ]
res$peak <- rownames(res)
peak2gene <- data.frame(
	peak = annotated_peaks$peak, 
	gene = annotated_peaks$gene_name
)
res_anno <- merge(
	res, 
	peak2gene, 
	by="peak", 
	all.x=TRUE
)
res_gene <- res_anno[c(10, 3, 7)]
res_gene<-res_gene[!is.na(res_gene$gene), ]

# now we can extract differentially expressed genes
res_sel <- res_gene[res_gene$padj<0.05 & abs(res_gene$log2FoldChange)>=0.5, ]
res_sel <- res_sel[!is.na(res_sel$gene), ] #remove NA
write.table(res_sel, file="results/res_sel", sep="\t", quote=F)
dim(res_sel)

write.table(res_gene, file="results/gene_FC_padj.txt", sep="\t", col.names=F, quote=F)
write.table(res_gene, file="results/ATAC_genelist.txt", sep="\t", quote=F)

# 2. Volcano plot
analysis_id<-"res_gene"
p_value_thr<-0.05
log2_fc_thr<- 0.5
df <- res_gene
filename <- "PRJNA725028_volcano_plot"
df_clean<-df[complete.cases(df), ]
df_clean<-df_clean[df_clean$padj > 0,]
df_clean$diffexpressed <- "NO"
df_clean$diffexpressed[df_clean$log2FoldChange > log2_fc_thr & df_clean$padj < p_value_thr] <- "UP"
df_clean$diffexpressed[df_clean$log2FoldChange < -log2_fc_thr & df_clean$padj < p_value_thr] <- "DOWN"
df_clean<-rbind(df_clean[df_clean$diffexpressed == "NO",], df_clean[df_clean$diffexpressed == "UP",], df_clean[df_clean$diffexpressed == "DOWN",])
df_clean$diffexpressed <- factor(
  df_clean$diffexpressed,
  levels = c("DOWN", "NO", "UP")
)
color_values <- c(
  "DOWN" = "#00AFBB",
  "NO"   = "grey",
  "UP"   = "#bb0c00"
)
df_clean$neglogp <- -log10(df_clean$padj)
CAP_Y <- 50                 # you can set 30–100 based on dataset
df_clean$neglogp_cap <- pmin(df_clean$neglogp, CAP_Y)

xlimit <- ceiling(max(c(
	abs(max(df_clean$log2FoldChange)),
	abs(min(df_clean$log2FoldChange))
)) / 2.5) * 2.5
ylimit <- -log10(min(df_clean$padj))

p<-ggplot(
	data = df_clean, 
	aes(
		x = log2FoldChange,
		y = -log10(padj),
		# y = neglogp_cap,
		col = diffexpressed
		# label = .data[["delabel"]]
	)
) +

geom_vline(xintercept = c(-log2_fc_thr, log2_fc_thr), col = "gray", linetype = 'dashed') + #border up and down
geom_hline(yintercept = -log10(p_value_thr), col = "gray", linetype = 'dashed')+ #significatività 

geom_point(size = 0.5) + #desing a gene like a 0.5 size dot

scale_color_manual( #color the dots
	values = color_values
) +
coord_cartesian(ylim = c(0, ylimit), xlim = c(-xlimit, xlimit)) + #limits of the axes
labs( #label
	color = 'Expression', #legend_title
	x = expression("log"[2]*"FC"),
	y = expression("-log"[10]*"p-value adj.")
) +

ggtitle("PRJNA725028") #assign the title to the plot

print(p)

filepath<-file.path("results",filename)

ggsave(
	filename = paste0(filepath,".png"),
	plot = p,
	dpi = 300,
	width = 1800, height = 2400, units = "px"
)

ggsave(
	filename = paste0(filepath,".tiff"),
	plot = p,
	dpi = 300,
	width = 1800, height = 2400, units = "px",
	compression = "lzw"
)

ggsave(
	filename = paste0(filepath,".pdf"),
	plot = p,
	width = 1800, height = 2400, units = "px"
)

up_genes<-df_clean$gene[df_clean$diffexpressed=="UP"]
down_genes<-df_clean$gene[df_clean$diffexpressed=="DOWN"]

length(up_genes)
length(down_genes)
