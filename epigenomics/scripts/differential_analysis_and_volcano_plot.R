#!/bin/bash

# Methylation array differential analysis and volcano plot
# Dataset:GSE48325

library(limma)
library(ggplot2)


# 1. Differential analysis
group <- pData(gset)$group
names(group)<-pData(gset)$geo_accession
group<-group[colnames(Mval)]

group<-factor(group, levels=c("NORMAL", "NASH"))
design<-model.matrix(~ group)
rownames(design) <- colnames(Mval)
all(colnames(Mval)==rownames(design)) #must be TRUE

# annotation
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
colnames(anno450k)

fit <-lmFit(
	Mval,
	design
)
fit<-eBayes(fit)
summary(decideTests(fit)) #Down, NotSig and Up
res<-topTable(
	fit, 
	number=Inf, 
	coef=2, 
	adjust="BH"
)
resf<-res[res$adj.P.Val < 0.05 & abs(res$logFC) > 0.5, ]
resf$CpG<-rownames(resf)
anno450kSub<-merge(
	resf, 
	anno450k, 
	by.x="CpG", 
	by.y="Name"
)

anno_clean <- anno450kSub[anno450kSub$UCSC_RefGene_Name != "", ]
anno_clean$Gene <- sapply(strsplit(anno_clean$UCSC_RefGene_Name, ";"), `[`, 1)

res_gene <- anno_clean[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val", "B", "CpG")]

write.table(res_gene, file="/methylation/reads/GSE48325/res_gene.rnk", sep="\t", quote=FALSE, row.names=FALSE)

res_gene_rank <- res_gene[order(res_gene$adj.P.Val),]
write.table(res_gene_rank, file="methylation/reads/GSE48325/res_gene_rank.rnk", sep="\t", quote=FALSE, row.names=FALSE)

anno_df <- as.data.frame(anno_clean)
anno_filtered <- anno_df[abs(anno_df$logFC) > 0.5 & anno_df$Gene !="" & anno_df$Relation_to_Island %in% c("OpenSea", "Island"),]

##separate differentially methylated regions between Island and OpenSea
anno_filtered<-anno_filtered[, c("Gene", "logFC", "Relation_to_Island")]
anno_filt_opensea <- anno_filtered[anno_filtered$Relation_to_Island == "OpenSea",] 
anno_filt_island <- anno_filtered[anno_filtered$Relation_to_Island == "Island", ] 

write.table(anno_filtered, file="methylation_gene_logFC.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(anno_filt_opensea, file="methylation_opensea_gene_logFC.txt", sep="\t", quote=F, col.names=T, row.names=F)
write.table(anno_filt_island, file="methylation_island_gene_logFC.txt", sep="\t", quote=F, col.names=T, row.names=F)


# 2. Volcano plot
analysis_id<-"res"
p_value_thr<-0.05
log2_fc_thr<- 0.5
df <- res
filename <- "GSE48325_volcano_plot"
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

ggtitle("GSE48325") #assign the title to the plot

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



df_clean$CpG <- rownames(df_clean)

df_clean <- merge(
  df_clean,
  anno_clean[, c("CpG", "Gene", "Relation_to_Island")],
  by = "CpG",
  all.x = TRUE
)

df_clean_filt <- na.omit(df_clean)

up_genes<-df_clean_filt$Gene[df_clean_filt$diffexpressed=="UP"]
down_genes<-df_clean_filt$Gene[df_clean_filt$diffexpressed=="DOWN"]

