#!/bin/bash

# miRNA analysis: differential analysis and volcano plot
# Dataset:GSE59492

library(limma)
library(ggplot2)

# 1. Differential analysis

design <- model.matrix(~ group)
fit <- lmFit(
	expr, 
	design
)
fit <- eBayes(fit)
res <- topTable(
	fit, 
	coef="groupNASH-CH", 
	number=Inf, 
	adjust.method = "BH"
)
# annotation 
mirtar <- read.delim(
	"mirtarbase_sel.tsv", 
	header=T
)
mirna_sel <- mirtar[mirtar$miRNA %in% rownames(res)]

res$miRNA_clean <- gsub("_st", "", rownames(res))
res$miRNA_clean <- gsub("-star", "-3p", res$miRNA_clean) #-star = -3p
res_gene <- merge(
	res, 
	mirtar, 
	by.x ="miRNA_clean", 
	by.y ="miRNA"
)
res_gene <- res_gene[, c("Target.Gene", "logFC", "adj.P.Val", "t", "miRNA_clean", "miRTarBase.ID")]
res_gene <- res_gene[grepl("^hsa-", res_gene$miRNA_clean),]
res_gene_filt<-res_gene[res_gene$adj.P.Val < 0.05 & abs(res_gene$logFC) > 0.5, ] #599

write.table(res_gene_filt, file="mirna_gene_logFC", sep="\t", quote=F, col.names=T, row.names=F)

# 2. Volcano plot
analysis_id<-"res_gene"
p_value_thr<-0.05
log2_fc_thr<- 0.5
df <- res_gene
filename <- "GSE59492_volcano_plot"
df_clean<-df[complete.cases(df), ]
df_clean<-df_clean[df_clean$adj.P.Val > 0,]
df_clean$diffexpressed <- "NO"
df_clean$diffexpressed[df_clean$logFC > log2_fc_thr & df_clean$adj.P.Val < p_value_thr] <- "UP"
df_clean$diffexpressed[df_clean$logFC < -log2_fc_thr & df_clean$adj.P.Val < p_value_thr] <- "DOWN"
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
df_clean$neglogp <- -log10(df_clean$adj.P.Val)
CAP_Y <- 50                 # you can set 30–100 based on dataset
df_clean$neglogp_cap <- pmin(df_clean$neglogp, CAP_Y)

xlimit <- ceiling(max(c(
	abs(max(df_clean$logFC)),
	abs(min(df_clean$logFC))
)) / 2.5) * 2.5
ylimit <- -log10(min(df_clean$adj.P.Val))

p<-ggplot(
	data = df_clean, 
	aes(
		x = logFC,
		y = -log10(adj.P.Val),
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

ggtitle("GSE59492") #assign the title to the plot

print(p)

filepath<-file.path(filename)

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

up_genes<-df_clean$Target.Gene[df_clean$diffexpressed=="UP"]
down_genes<-df_clean$Target.Gene[df_clean$diffexpressed=="DOWN"]

length(up_genes)
length(down_genes)
#these are up- and down- regulated genes not the miRNA!
