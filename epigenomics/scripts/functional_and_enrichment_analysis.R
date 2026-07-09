#!/bin/bash

# Methylation array functional analysis
# Dataset:GSE48325

library(clusterProfiler)
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(IlluminaHumanMethylation450kmanifest)
library(IlluminaHumanMethylation450kanno.ilmn12.hg19)
library(dplyr)
library(KEGGREST)
library(VennDiagram)

# 1. ORA
log2fc_list <- list()
results_list <- list()
enrich_list <- list()

heat_colors <- colorRampPalette(c("green", "black", "red"))(150)
gene_sets <- clusterProfiler::read.gmt("/ref/c5.go.bp.v2024.1.Hs.symbols.gmt")

results_list[["NASH_vs_NORMAL"]] <- res_gene_rank
log2fc_list[["NASH_vs_NORMAL"]] <- res_gene_rank$logFC
names(log2fc_list[["NASH_vs_NORMAL"]]) <- res_gene_rank$Gene

#split genes by direction
up_genes <- res_gene_rank$Gene[res_gene_rank$logFC > 0]
down_genes <- res_gene_rank$Gene[res_gene_rank$logFC < 0]

up_genes_island <- anno_filt_island$Gene[anno_filt_island$logFC > 0]
down_genes_island <- anno_filt_island$Gene[anno_filt_island$logFC < 0]

up_genes_opensea<- anno_filt_opensea$Gene[anno_filt_opensea$logFC > 0]
down_genes_opensea <- anno_filt_opensea$Gene[anno_filt_opensea$logFC < 0]

#run enrichment
#all
enriched_up <- clusterProfiler::enricher(
	gene = up_genes,
	TERM2GENE = gene_sets
) 

enriched_down <- clusterProfiler::enricher(
	gene = down_genes,
	TERM2GENE = gene_sets
)

enriched_all <- clusterProfiler::enricher(
	gene = res_gene_rank$Gene,
	TERM2GENE = gene_sets
)

#island
enriched_up_island <- clusterProfiler::enricher(
	gene = up_genes_island,
	TERM2GENE = gene_sets
)

enriched_down_island <- clusterProfiler::enricher(
	gene = down_genes_island,
	TERM2GENE = gene_sets
)

enriched_all_island <- clusterProfiler::enricher(
	gene = anno_filt_island$Gene,
	TERM2GENE = gene_sets
)

#opensea
enriched_up_opensea <- clusterProfiler::enricher(
	gene = up_genes_opensea,
	TERM2GENE = gene_sets
)

enriched_down_opensea <- clusterProfiler::enricher(
	gene = down_genes_opensea,
	TERM2GENE = gene_sets
)

enriched_all_opensea <- clusterProfiler::enricher(
	gene = anno_filt_opensea$Gene,
	TERM2GENE = gene_sets
)

#general pathways
if (!is.null(enriched_all) && nrow(enriched_all@result) > 0) { 
      all_df <- data.frame(
         ID = enriched_all@result$ID,
         Count = enriched_all@result$Count,
         GeneRatio = enriched_all@result$GeneRatio,
         pvalue = enriched_all@result$pvalue,
         padj = enriched_all@result$p.adjust
  )

  prefix <- paste0("EN_ALL_NASH_vs_NORMAL")
  write.csv(all_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

  enrich_list[["NASH_vs_NORMAL_ALL"]] <- all_df$ID[all_df$padj < 0.05]

} else {
  enrich_list[["NASH_vs_NORMAL_ALL"]] <- character(0)
}

#Up pathtway
if (!is.null(enriched_up) && nrow(enriched_up@result) > 0) {
      up_df <- data.frame(
        ID = enriched_up@result$ID,
        Count = enriched_up@result$Count,
        GeneRatio = enriched_up@result$GeneRatio,
        pvalue = enriched_up@result$pvalue,
        padj = enriched_up@result$p.adjust
      )

      prefix <- paste0("EN_UP_NASH_vs_NORMAL")
      write.csv(up_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_UP"]] <- up_df$ID[up_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_UP"]] <- character(0)
    }

#Down pathway
 if (!is.null(enriched_down) && nrow(enriched_down@result) > 0) {
      down_df <- data.frame(
        ID = enriched_down@result$ID,
        Count = enriched_down@result$Count,
        GeneRatio = enriched_down@result$GeneRatio,
        pvalue = enriched_down@result$pvalue,
        padj = enriched_down@result$p.adjust
      )

      prefix <- paste0("EN_DOWN_NASH_vs_NORMAL")
      write.csv(down_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- down_df$ID[down_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- character(0)
    }
  
#general pathways
if (!is.null(enriched_all_island) && nrow(enriched_all_island@result) > 0) { 
      all_df <- data.frame(
         ID = enriched_all_island@result$ID,
         Count = enriched_all_island@result$Count,
         GeneRatio = enriched_all_island@result$GeneRatio,
         pvalue = enriched_all_island@result$pvalue,
         padj = enriched_all_island@result$p.adjust
  )

  prefix <- paste0("EN_ALL_NASH_vs_NORMAL_island")
  write.csv(all_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

  enrich_list[["NASH_vs_NORMAL_ALL"]] <- all_df$ID[all_df$padj < 0.05]

} else {
  enrich_list[["NASH_vs_NORMAL_ALL"]] <- character(0)
}

#Up pathtway
if (!is.null(enriched_up_island) && nrow(enriched_up_island@result) > 0) {
      up_df <- data.frame(
        ID = enriched_up_island@result$ID,
        Count = enriched_up_island@result$Count,
        GeneRatio = enriched_up_island@result$GeneRatio,
        pvalue = enriched_up_island@result$pvalue,
        padj = enriched_up_island@result$p.adjust
      )

      prefix <- paste0("EN_UP_NASH_vs_NORMAL_island")
      write.csv(up_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_UP"]] <- up_df$ID[up_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_UP"]] <- character(0)
    }

#Down pathway
 if (!is.null(enriched_down_island) && nrow(enriched_down_island@result) > 0) {
      down_df <- data.frame(
        ID = enriched_down_island@result$ID,
        Count = enriched_down_island@result$Count,
        GeneRatio = enriched_down_island@result$GeneRatio,
        pvalue = enriched_down_island@result$pvalue,
        padj = enriched_down_island@result$p.adjust
      )

      prefix <- paste0("EN_DOWN_NASH_vs_NORMAL_island")
      write.csv(down_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- down_df$ID[down_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- character(0)
    }

#general pathways
if (!is.null(enriched_all_opensea) && nrow(enriched_all_opensea@result) > 0) { 
      all_df <- data.frame(
         ID = enriched_all_opensea@result$ID,
         Count = enriched_all_opensea@result$Count,
         GeneRatio = enriched_all_opensea@result$GeneRatio,
         pvalue = enriched_all_opensea@result$pvalue,
         padj = enriched_all_opensea@result$p.adjust
  )

  prefix <- paste0("EN_ALL_NASH_vs_NORMAL_opensea")
  write.csv(all_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

  enrich_list[["NASH_vs_NORMAL_ALL"]] <- all_df$ID[all_df$padj < 0.05]

} else {
  enrich_list[["NASH_vs_NORMAL_ALL"]] <- character(0)
}

#Up pathtway
if (!is.null(enriched_up_opensea) && nrow(enriched_up_opensea@result) > 0) {
      up_df <- data.frame(
        ID = enriched_up_opensea@result$ID,
        Count = enriched_up_opensea@result$Count,
        GeneRatio = enriched_up_opensea@result$GeneRatio,
        pvalue = enriched_up_opensea@result$pvalue,
        padj = enriched_up_opensea@result$p.adjust
      )

      prefix <- paste0("EN_UP_NASH_vs_NORMAL_opensea")
      write.csv(up_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_UP"]] <- up_df$ID[up_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_UP"]] <- character(0)
    }

#Down pathway
 if (!is.null(enriched_down_opensea) && nrow(enriched_down_opensea@result) > 0) {
      down_df <- data.frame(
        ID = enriched_down_opensea@result$ID,
        Count = enriched_down_opensea@result$Count,
        GeneRatio = enriched_down_opensea@result$GeneRatio,
        pvalue = enriched_down_opensea@result$pvalue,
        padj = enriched_down_opensea@result$p.adjust
      )

      prefix <- paste0("EN_DOWN_NASH_vs_NORMAL_opensea")
      write.csv(down_df, paste0("DE/", prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- down_df$ID[down_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- character(0)
    }


#Pathway presence/absence between contrasts
all_pathways <-c(enrich_list[["NASH_vs_NORMAL_UP"]], enrich_list[["NASH_vs_NORMAL_DOWN"]], enrich_list[["NASH_vs_NORMAL_ALL"]])

pathway_matrix <- data.frame(Pathway = all_pathways, UP = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_UP"]]), DOWN = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_DOWN"]]), ALL = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_ALL"]]))

head(pathway_matrix)

# dotplot
png("plots/dotplot_enrich_down_island.png", width=1200, height=1200, res=150)
dotplot(enriched_down_island)
dev.off()

# barplot
png("plots/barplot_all_island.png", width=1200, height=1200, res=150)
barplot(
	enriched_all_island, 
	title = "Island GO Biological Process Enriched"
)
dev.off()

write.csv(pathway_matrix, "results/EN_NASH_vs_NORMAL_pathway_presence_absence.csv", row.names = FALSE)

genes<- res_gene_rank$Gene
write.csv(genes, paste0("results/genes_NASH_vs_NORMAL"), row.names=FALSE)

# 2. GSEA
#set the organism
organism="org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#prepare input
res$CpG<-rownames(res)
anno450kSub<-merge(
	res, 
	anno450k, 
	by.x="CpG", 
	by.y="Name"
)
anno_clean <- anno450kSub[anno450kSub$UCSC_RefGene_Name != "", ]
anno_clean$Gene <- sapply(strsplit(anno_clean$UCSC_RefGene_Name, ";"), `[`, 1)
res_gene <- anno_clean[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val", "B", "CpG", "Relation_to_Island")]

res_gene_unique <- aggregate(t ~ Gene, data = res_gene, FUN = max)
original_gene_list <- res_gene_unique$t
names(original_gene_list) <- res_gene_unique$Gene
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing =TRUE)

# gene set enrichment
gse <- gseGO(
	geneList=gene_list, 
	ont="BP", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

# GSEA plot
png("plots/GSEA.png", width = 1200, height = 800)
gseaplot(
	gse , 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

# GSEA sign p val
res$sign_p <- res$adj.P.Val * sign(res$logFC)
res$CpG<-rownames(res)
anno450k <- getAnnotation(IlluminaHumanMethylation450kanno.ilmn12.hg19)
anno450kSub<-merge(res, anno450k, by.x="CpG", by.y="Name")
anno_clean <- anno450kSub[anno450kSub$UCSC_RefGene_Name != "", ]
anno_clean$Gene <- sapply(strsplit(anno_clean$UCSC_RefGene_Name, ";"), `[`, 1)
head(anno_clean$Gene)
res_gene_sign <- anno_clean[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val", "B", "sign_p", "CpG")]

res_gene_unique_2 <- aggregate(sign_p ~ Gene, data = res_gene_sign, FUN = max)
original_gene_list_2 <- res_gene_unique_2$sign_p
names(original_gene_list_2) <- res_gene_unique_2$Gene
gene_list_2 <- na.omit(original_gene_list_2)
gene_list_2 <- sort(gene_list_2, decreasing =TRUE)

# gene set enrichment
gse_2 <- gseGO(
	geneList=gene_list_2, 
	ont="ALL", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

#GSEA plot
png("plots/GSEA_signp.png", width = 1200, height = 800)
gseaplot(gse_2 , by="all", title=gse_2$Description[1], geneSetID=1)
dev.off()

# GSEA logFC
res_gene_unique_fc <- aggregate(logFC ~ Gene, data = res_gene, FUN = max)
original_gene_list_fc <- res_gene_unique_fc$logFC
names(original_gene_list_fc) <- res_gene_unique_fc$Gene
gene_list_fc <- na.omit(original_gene_list_fc)
gene_list_fc <- sort(gene_list_fc, decreasing =TRUE)

#gene set enrichment
gse_fc <- gseGO(
	geneList=gene_list_fc, 
	ont="ALL", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

#GSEA plot
png("plots/GSEA_fc.png", width = 1200, height = 800)
gseaplot(gse_fc , by="all", title=gse_fc$Description[1], geneSetID=1)
dev.off()

##GSEA random
set.seed(56)
gene_list_random <- sample(gene_list)
names(gene_list_random) <- names(gene_list)
gene_list_random <- sort(gene_list_random, decreasing=TRUE)

#gene set enrichment
gse_random <- gseGO(
	geneList=gene_list_random, 
	ont="ALL", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

#GSEA plot
png("plots/GSEA_random.png", width = 1200, height = 800)
gseaplot(gse_random , by="all", title=gse_random$Description[1], geneSetID=1)
dev.off()

#GSEA with only Island and OpenSea regions
#select Island and OpenSea
res_gene_Is <- res_gene[res_gene$Relation_to_Island == "Island", ]
res_gene_OpS <- res_gene[res_gene$Relation_to_Island == "OpenSea", ]

##check
dim(res_gene_Is)

dim(res_gene_OpS)

#gene list Island
res_gene_Is_unique <- aggregate(t ~ Gene, data = res_gene_Is, FUN = max)
gene_list_Is <- res_gene_Is_unique$t
names(gene_list_Is) <- res_gene_Is_unique$Gene
gene_list_Is <- na.omit(gene_list_Is)
gene_list_Is <- sort(gene_list_Is, decreasing =TRUE)

#gene set enrichment Island
gse_Is <- gseGO(
	geneList=gene_list_Is, 
	ont="BP", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

#GSEA plot Island
png("plots/GSEA_Is.png", width = 1200, height = 800)
gseaplot(
	gse_Is, 
	by="all", 
	title=gse_Is$Description[1], 
	geneSetID=1
)
dev.off()

#gene list OpenSea
res_gene_OpS_unique <- aggregate(t ~ Gene, data = res_gene_OpS, FUN = max)
gene_list_OpS <- res_gene_OpS_unique$t
names(gene_list_OpS) <- res_gene_OpS_unique$Gene
gene_list_OpS <- na.omit(gene_list_OpS)
gene_list_OpS <- sort(gene_list_OpS, decreasing =TRUE)

#gene set enrichment Opensea
gse_OpS <- gseGO(
	geneList=gene_list_OpS, 
	ont="BP", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

#GSEA plot OpenSea
png("plots/GSEA_OpS.png", width = 1200, height = 800)
gseaplot(
	gse_OpS, 
	by="all", 
	title=gse_OpS$Description[1], 
	geneSetID=1
)
dev.off()

write.csv(gse, file="results/GSEA_results_all.csv", row.names=F)
write.csv(gse_Is, file="results/GSEA_results_island_BP", row.names=F)
write.csv(gse_OpS, file="results/GSEA_results_open_sea_BP", row.names=F)


# common pathway between island and opensea
common_path <- intersect(
	gse_Is$ID, 
	gse_OpS$ID
)
#length(common_path)

path_island <- setdiff(gse_Is$ID, gse_OpS$ID)
#length(path_island)
path_opensea <- setdiff(gse_OpS$ID, gse_Is$ID)
#length(path_opensea)

#compare path
all_path <- unique(c(gse_Is$ID, gse_OpS$ID))
#length(all_path)

comparison_path <- data.frame(Pathway = all_path, Island = as.integer(all_path %in% results_is_BP$ID), OpenSea = as.integer(all_path %in% results_ops_BP$ID))

write.csv(comparison_path, "results/GSEA_BP_island_vs_opensea.csv", row.names=F)

# ridgeplot
png("plots/ridgeplot_all.png", width = 1200, height = 800)
ridgeplot(gse, showCategory=15) + 
labs(x = "Enrichment distribution") + 
theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14))
dev.off()

png("plots/ridgeplot_island.png", width = 1200, height = 800)
ridgeplot(gse_Is, showCategory=15) + 
labs(x = "Island enrichment distribution") + 
theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14))
dev.off()

png("plots/ridgeplot_opensea.png", width = 1200, height = 800)
ridgeplot(gse_OpS, showCategory=15) + 
labs(x = "Opensea enriachment distribution") + 
theme(axis.text.y = element_text(size = 10), axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14))
dev.off()

# dotplot
png("plots/dotplot_gsea.png", width=1200, height=1200, res=150)
dotplot(
	gse, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10,
	title = "Biological Process Enriched"
)
dev.off()

png("plots/dotplot_is_gsea.png", width=1200, height=1200, res=150)
dotplot(
	gse_Is, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10,
	title = "Island Biological Process Enriched"
)
dev.off()

png("plots/dotplot_ops_gsea.png", width=1200, height=1200, res=150)
dotplot(
	gse_OpS, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10,
	title = "OpenSea Biological Process Enriched"
)
dev.off()

# methylated and unmethylated regions

#order res based on chromosome and position
res_dmr <- anno_clean[, c("Gene", "logFC", "t", "P.Value", "adj.P.Val", "B", "CpG", "Relation_to_Island", "chr", "pos")]
res_dmr <- res_dmr[order(res_dmr$chr, res_dmr$pos), ]

#create cluster
breaks <- c(TRUE, diff(res_dmr$pos) > 500 | res_dmr$chr[-1] != res_dmr$chr[-length(res_dmr$chr)])

res_dmr$cluster <- cumsum(breaks)

#define regions
res_dmr <- as.data.frame(res_dmr)
regions <- res_dmr %>% group_by(chr, cluster) %>% summarise(start=min(pos), end=max(pos), mean_fc=mean(logFC), .groups="drop")
regions$length <- regions$end - regions$start

#separate hypo- and hyper- methylated regions
hypo <- regions[regions$mean_fc < 0, ]
hyper <- regions[regions$mean_fc > 0, ]

#boxplot
png("plots/boxplot.png", width = 1200, height = 800)
boxplot(
	hypo$length, 
	hyper$length, 
	horizzontal=F, 
	col=c("green", "coral"), 
	names=c("Hypo", "Hyper")
)
dev.off()

# 3. KEGG
#all
Human.GRCh38.p13.annot <- read.table("/ref/Human.GRCh38.p13.annot.tsv", header=T, sep="\t", fill= T, comment.char ="#")

res_gene_entrez <- merge(
	res_gene_unique, 
	Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], 
	by.x ="Gene", 
	by.y = "Symbol", 
	all.x =T
)

## filter NA values
res_gene_entrez_filt <- res_gene_entrez[!is.na(res_gene_entrez$GeneID), ]

# KEGG analysis
#https://www.kegg.jp/pathway/hsa04932
#pathway=hsa04932

path <- keggGet("hsa04932") #get the pathway
nafld_genes<- path[[1]]$GENE #extract the genes

nafld_entrez <- nafld_genes[seq(1, length(nafld_genes), 2)] #extract the entrez
nafld_entrez <- sub(";.*", "", nafld_entrez) #remove ;

#overlap to see if the gene I obtained are involved in the pathology
gene_list <- unique(na.omit(res_gene_entrez$GeneID)) #extract the gene obtained
overlap <- intersect(
	gene_list, 
	nafld_entrez
)
length(overlap)

png("plots/barplot_kegg.png", width = 1200, height = 800)
barplot(c(Overlap = length(overlap),Not_overlap = length(gene_list) - length(overlap)))
dev.off()

png("plots/venn_kegg.png", width = 800, height = 800)
venn.plot <- venn.diagram(
	x = list(
		"Epigenomic" = gene_list, 
		"NAFLD pathway" = nafld_entrez
	),
	filename = NULL, 
	scaled=FALSE, 
	fill = c("#FF9900", "#3276FF"), 
	alpha = 0.5, 
	cex = 2, 
	cat.cex = 2, 
	cat.pos = c(-20, 20), 
	cat.dist = c(0.05, 0.05), 
	main = "Overlap Epigenomic vs KEGG NAFLD", 
	main.cex=2
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

overlap_symbols <- Human.GRCh38.p13.annot[Human.GRCh38.p13.annot$GeneID %in% overlap, c("GeneID", "Symbol", "Description")]
overlap_symbols <- overlap_symbols[order(overlap_symbols$Symbol), ]
write.csv(overlap_symbols, "NASH_overlap_genes.csv", row.names = FALSE, quote=F)

#Island
res_gene_entrez_is <- merge(
	res_gene_Is_unique, 
	Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], 
	by.x ="Gene", 
	by.y = "Symbol", 
	all.x =T
)

## filter NA values
res_gene_entrez_filt_is <- res_gene_entrez_is[!is.na(res_gene_entrez_is$GeneID), ]

# KEGG analysis
#https://www.kegg.jp/pathway/hsa04932
#pathway=hsa04932

#overlap to see if the gene I obtained are involved in the pathology
gene_list_is <- unique(na.omit(res_gene_entrez_is$GeneID)) #extract the gene obtained
overlap_is <- intersect(
	gene_list_is, 
	nafld_entrez
)
length(overlap_is)

png("plots/barplot_is_kegg.png", width = 1200, height = 800)
barplot(c(Overlap = length(overlap_is),Not_overlap = length(gene_list_is) - length(overlap_is)))
dev.off()

png("plots/venn_is_kegg.png", width = 800, height = 800)
venn.plot <- venn.diagram(
	x = list(
		"Island Methylation" = gene_list_is, 
		"NAFLD pathway" = nafld_entrez
	),
	filename = NULL, 
	scaled=FALSE, 
	fill = c("#FF9900", "#3276FF"), 
	alpha = 0.5, 
	cex = 2, 
	cat.cex = 2, 
	cat.pos = c(-20, 20), 
	cat.dist = c(0.05, 0.05), 
	main = "Overlap Island Methylation vs KEGG NAFLD", 
	main.cex=2
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

overlap_symbols_is <- Human.GRCh38.p13.annot[Human.GRCh38.p13.annot$GeneID %in% overlap_is, c("GeneID", "Symbol", "Description")]
overlap_symbols_is <- overlap_symbols_is[order(overlap_symbols_is$Symbol), ]
write.csv(overlap_symbols_is, "NASH_overlap_genes_is.csv", row.names = FALSE, quote=F)

#OpenSea
res_gene_entrez_ops <- merge(
	res_gene_OpS_unique, 
	Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], 
	by.x ="Gene", 
	by.y = "Symbol", 
	all.x =T
)

## filter NA values
res_gene_entrez_filt_ops <- res_gene_entrez_ops[!is.na(res_gene_entrez_ops$GeneID), ]

# KEGG analysis
#https://www.kegg.jp/pathway/hsa04932
#pathway=hsa04932

#overlap to see if the gene I obtained are involved in the pathology
gene_list_ops <- unique(na.omit(res_gene_entrez_ops$GeneID)) #extract the gene obtained
overlap_ops <- intersect(gene_list_ops, nafld_entrez)
length(overlap_ops)

png("plots/barplot_ops_kegg.png", width = 1200, height = 800)
barplot(c(Overlap = length(overlap_ops),Not_overlap = length(gene_list_ops) - length(overlap_ops)))
dev.off()

png("plots/venn_ops_kegg.png", width = 800, height = 800)
venn.plot <- venn.diagram(
	x = list(
		"Open Sea Methylation" = gene_list_ops, 
		"NAFLD pathway" = nafld_entrez
	),
	filename = NULL, 
	scaled=FALSE, 
	fill = c("#FF9900", "#3276FF"), 
	alpha = 0.5, 
	cex = 2, 
	cat.cex = 2, 
	cat.pos = c(-10, 10), 
	cat.dist = c(0.05, 0.05), 
	main = "Overlap Open Sea Methylation vs KEGG NAFLD", 
	main.cex=2
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

overlap_symbols_ops <- Human.GRCh38.p13.annot[Human.GRCh38.p13.annot$GeneID %in% overlap_ops, c("GeneID", "Symbol", "Description")]
overlap_symbols_ops <- overlap_symbols_ops[order(overlap_symbols_ops$Symbol), ]
write.csv(overlap_symbols_ops, "NASH_overlap_genes_ops.csv", row.names = FALSE, quote=F)