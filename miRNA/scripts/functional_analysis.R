#!/bin/bash

# miRNA analysis: functional analysis
# Dataset:GSE59492

library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(org.Hs.eg.db)
library(KEGGREST)
library(VennDiagram)

# 1. ORA
log2fc_list <- list()
results_list <- list()
enrich_list <- list()

heat_colors <- colorRampPalette(c("green", "black", "red"))(150)
gene_sets <- clusterProfiler::read.gmt("/media/matteo/SSD5/lucrezia/ref/c5.go.bp.v2024.1.Hs.symbols.gmt")

results_list[["NASH_vs_NORMAL"]] <- res_gene_filt
log2fc_list[["NASH_vs_NORMAL"]] <- res_gene_filt$logFC
names(log2fc_list[["NASH_vs_NORMAL"]]) <- res_gene_filt$Target.Gene

# split genes by direction
up_genes <- res_gene_filt$Target.Gene[res_gene_filt$logFC > 0]
down_genes <- res_gene_filt$Target.Gene[res_gene_filt$logFC < 0]

# check
length(up_genes) #mirna up -> target gene down
length(down_genes) #mirna down -> target gene up

#run enrichment
enriched_up <- clusterProfiler::enricher(
	gene = up_genes,
	TERM2GENE = gene_sets
)

enriched_down <- clusterProfiler::enricher(
	gene = down_genes,
	TERM2GENE = gene_sets
)

enriched_all <- clusterProfiler::enricher(
	gene = res_gene_filt$Target.Gene,
	TERM2GENE = gene_sets
)

# check
dim(enriched_up)
dim(enriched_down)
dim(enriched_all) 

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
  write.csv(all_df, paste0(prefix, ".csv"), row.names = FALSE)

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
      write.csv(up_df, paste0(prefix, ".csv"), row.names = FALSE)

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
      write.csv(down_df, paste0(prefix, ".csv"), row.names = FALSE)

      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- down_df$ID[down_df$padj < 0.05]

    } else {
      enrich_list[["NASH_vs_NORMAL_DOWN"]] <- character(0)
    }
  
#Pathway presence/absence between contrasts
all_pathways <-c(enrich_list[["NASH_vs_NORMAL_UP"]], enrich_list[["NASH_vs_NORMAL_DOWN"]], enrich_list[["NASH_vs_NORMAL_ALL"]])

pathway_matrix <- data.frame(
	Pathway = all_pathways, 
	UP = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_UP"]]), 
	DOWN = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_DOWN"]]), 
	ALL = as.integer(all_pathways %in% enrich_list[["NASH_vs_NORMAL_ALL"]])
)

head(pathway_matrix)

write.csv(pathway_matrix, "EN_NASH_vs_NORMAL_pathway_presence_absence.csv", row.names = FALSE)

# dotplot
png("plots/dotplot_enrich.png", width=1200, height=1200, res=150)
dotplot(enriched_all)
dev.off()

png("plots/barplot_up.png", width=1200, height=1200, res=150)
barplot(enriched_all, title = "GO Biological Process Enriched")
dev.off()

# 2. GSEA
#set the organism
organism="org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#prepare input
res_gene_unique <- aggregate(t ~ Target.Gene, data = res_gene, FUN = max)
original_gene_list <- res_gene_unique$t
names(original_gene_list) <- res_gene_unique$Target.Gene
gene_list <- na.omit(original_gene_list)
gene_list <- sort(gene_list, decreasing =TRUE)

#gene set enrichment
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

gse_MF <- gseGO(
	geneList=gene_list, 
	ont="MF", 
	keyType="SYMBOL", 
	nPerm=1000, 
	minGSSize =3, 
	maxGSSize=800, 
	pvalueCutoff=0.05, 
	verbose=TRUE, 
	OrgDb = organism, 
	pAdjustMethod="none"
)

gse_all <- gseGO(
	geneList=gene_list, 
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

# GSEA plot
png("plots/GSEA.png", width = 1200, height = 800)
gseaplot(
	gse, 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

png("plots/GSEA_MF.png", width = 1200, height = 800)
gseaplot(
	gse_MF, 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

png("plots/GSEA_all.png", width = 1200, height = 800)
gseaplot(
	gse_all, 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

png("plots/ridgeplot.png", width = 1200, height = 800)
ridgeplot(gse, showCategory=15) +  
labs(x = "enrichment distribution") + 
theme(
	axis.text.y = element_text(size =9), 
	axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14)
)
dev.off()

png("plots/ridgeplot_MF.png", width = 1200, height = 800)
ridgeplot(gse_MF) +  
labs(x = "enrichment distribution") + 
theme(
	axis.text.y = element_text(size = 10), 
	axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14)
)
dev.off()

png("plots/ridgeplot_all.png", width = 1200, height = 800)
ridgeplot(gse_all, showCategory = 15) +  
labs(x = "enrichment distribution") + 
theme(
	axis.text.y = element_text(size = 9), 
	axis.text.x = element_text(size = 10), axis.title.x = element_text(size = 14)
)
dev.off()

# dotplot
png("plots/dotplot_2.png", width=3000, height=3000, res=300)
dotplot(
	gse, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10
)
dev.off()

png("plots/dotplot_MF.png", width=3000, height=3000, res=300)
dotplot(
	gse_MF, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10
)
dev.off()

png("plots/dotplot_all.png", width=3000, height=3000, res=300)
dotplot(
	gse_all, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10
)
dev.off()

# 3. KEGG
Human.GRCh38.p13.annot <- read.table("/ref/Human.GRCh38.p13.annot.tsv", header=T, sep="\t", fill= T, comment.char ="#")
res_gene_entrez <- merge(res_gene_filt, Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], by.x ="Target.Gene", by.y = "Symbol", all.x =T)
res_gene_entrez_filt <- res_gene_entrez[!is.na(res_gene_entrez$GeneID), ] # filter NA values

#https://www.kegg.jp/pathway/hsa04932
#pathway=hsa04932

path <- keggGet("hsa04932") #get the pathway
nafld_genes<- path[[1]]$GENE #extract the genes
nafld_entrez <- nafld_genes[seq(1, length(nafld_genes), 2)] #extract the entrez
nafld_entrez <- sub(";.*", "", nafld_entrez) #remove ;

#overlap to see if the gene I obtained are involved in the pathology
gene_list <- unique(na.omit(res_gene_entrez$GeneID)) #extract the gene obtained 
overlap <- intersect(gene_list, nafld_entrez)
length(overlap)

png("plots/barplot_kegg.png", width = 1200, height = 800)
barplot(c(
	Overlap = length(overlap),
	Not_overlap = length(gene_list) - length(overlap))
)
dev.off()

png("plots/venn_kegg.png", width = 800, height = 800)
venn.plot <- venn.diagram(
	x = list(
		miRNA = gene_list, 
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
	main = "Overlap miRNA vs KEGG NAFLD", 
	main.cex=2
)
grid.newpage()
grid.draw(venn.plot)
dev.off()

overlap_symbols <- Human.GRCh38.p13.annot[Human.GRCh38.p13.annot$GeneID %in% overlap, c("GeneID", "Symbol", "Description")]
overlap_symbols <- overlap_symbols[order(overlap_symbols$Symbol), ]
write.csv(overlap_symbols, "NASH_overlap_genes.csv", row.names = FALSE, quote=F)
