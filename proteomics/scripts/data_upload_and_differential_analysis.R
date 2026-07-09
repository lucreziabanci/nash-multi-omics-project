#!/bin/bash

# Proteomics data from MaxQuant (v2.7.5.0) uploading and differential analysis
# Dataset:PXD026717

library(limma)
library(ggplot2)

# 1. Upload proteinGroups.txt file in R
proteingroups_table <- read.delim(
	"/PXD026717/combined/txt/proteinGroups.txt", 
	header=TRUE, 
	sep="\t", 
	quote=""
)

# remove "+" from Only.identified.by.site, Reverse and Potential.contaminat
data<-proteingroups_table[proteingroups_table$Only.identified.by.site !="+" & proteingroups_table$Reverse !="+" & proteingroups_table$Potential.contaminant !="+",] 

# remove Healthy_2 e Healthy_4
data_filtered <-data[, c(1, 38, 40, 42, 43, 44, 45)] 
data_filtered[ , -1] <- log2(data_filtered[ , -1])
data_filtered <- data_filtered[apply(data_filtered[ , -1], 1, var) != 0, ]
colnames(data_filtered)<-c("Protein_ID", "Healthy_1", "Healthy_3", "NASH_1", "NASH_2", "NASH_3", "NASH_4")

# 2. Differential analysis
groups<- factor(gsub("_", "", gsub("\\d", "", colnames(data_filtered)[-1])))
design<-model.matrix(~0+groups)
colnames(design) <- levels(groups)
contrast.matrix <- makeContrasts(
	NASH - Healthy,
	levels=design
)
fit <- lmFit(
	data_filtered[ , -1], 
	design
)
fit <- contrasts.fit(
	fit, 
	contrast.matrix
)
fit <- eBayes(fit) 
res<-topTable(
	fit,
	n=Inf
)

res_clean <- res[!is.na(res$logFC), ]
res_clean$Protein_ID <- data_filtered$Protein_ID[as.numeric(rownames(res_clean))]
res_clean$Uniprot <- sapply(strsplit(res_clean$Protein_ID, ";"), `[`, 1)
res_sig <- res_clean[res_clean$adj.P.Val < 0.05 & abs(res_clean$logFC) > 0.5, ]
res_sig <- res_sig[!is.na(res_sig$logFC) & !is.na(res_sig$t), ]
sig_proteins <- unique(res_sig$Uniprot) # extract significant proteins


# 3. Volcano plot
analysis_id<-"res_clean"
p_value_thr<-0.05
log2_fc_thr<- 0.5
df <- res_clean
filename <- "PXD026717_volcano_plot"
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

ggtitle("PXD026717") #assign the title to the plot

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

up_genes<-df_clean$Uniprot[df_clean$diffexpressed=="UP"]
down_genes<-df_clean$Uniprot[df_clean$diffexpressed=="DOWN"]

length(up_genes)

length(down_genes)
