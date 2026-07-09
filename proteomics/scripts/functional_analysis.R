#!/bin/bash

# Proteomics functional analysis
# Dataset:PXD026717

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
library(KEGGREST)

# 1. GSEA
# set the organism
organism="org.Hs.eg.db"
#BiocManager::install(organism, character.only = TRUE)
library(organism, character.only = TRUE)

#symbol <- bitr(sig_proteins, fromType = "UNIPROT", toType = "SYMBOL", OrgDb = organism) #convert the Uniprot name in a symbol
symbol_all <- bitr(
	res_3_clean$Uniprot, 
	fromType = "UNIPROT", 
	toType = "SYMBOL", 
	OrgDb = organism
) 

res_gene <- merge(
	res_sig, 
	symbol, 
	by.x="Uniprot", 
	by.y="UNIPROT", 
	all.x=TRUE
)
res_gene_all <- merge(
	res_clean, 
	symbol_all, 
	by.x="Uniprot", 
	by.y="UNIPROT", 
	all.x=TRUE
)

res_gene_clean <- res_gene[!is.na(res_gene$SYMBOL), ] #remove NA
res_gene_all_clean <- res_gene_all[!is.na(res_gene_all$SYMBOL), ] #remove NA

res_gene_sel<- res_gene_clean[, c("SYMBOL", "logFC")]
res_gene_all_sel <- res_gene_all_clean[, c("SYMBOL", "logFC")]

write.table(res_gene_sel, file="proteomic_symbol_logFC.txt", sep="\t", quote=F, col.names=T, row.names=F)

# prepare input
res_gene_unique <- aggregate(
	t ~ SYMBOL, 
	data = res_gene_all_clean, 
	FUN = max
)

original_gene_list <- res_gene_unique$t
names(original_gene_list) <- res_gene_unique$SYMBOL
gene_list <- na.omit(original_gene_list)
gene_list <- sort(
	gene_list, 
	decreasing =TRUE
)

## gene set enrichment
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

png("plots/GSEA.png", width = 1200, height = 800)
gseaplot(
	gse , 
	by="all", 
	title=gse$Description[1], 
	geneSetID=1
)
dev.off()

## GSEA random
set.seed(56)
gene_list_random <- sample(gene_list)
names(gene_list_random) <- names(gene_list)
gene_list_random <- sort(
	gene_list_random, 
	decreasing=TRUE
)

## gene set enrichment
gse_random <- gseGO(
	geneList=gene_list_random, 
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
png("plots/GSEA_random.png", width = 1200, height = 800)
gseaplot(
	gse_random, 
	by="all", 
	title=gse_random$Description[1], 
	geneSetID=1
)
dev.off()

# Ridgeplot
png("plots/ridgeplot.png", width = 1200, height = 800)
ridgeplot(gse) +  
labs(x = "enrichment distribution") + 
theme(
	axis.text.y = element_text(size = 8), 
	axis.text.x = element_text(size = 10), 
	axis.title.x = element_text(size = 14)
)
dev.off()

# Dotplot
png("plots/dotplot.png", width=3000, height=3000, res=300)
dotplot(
	gse, 
	showCategory=15, 
	orderBy="NES", 
	font.size=10,
	title = "Biological Process Enriched"
)
dev.off()

# Create the pathway matrix with up and down proteins
all_pathways <- gse@result
all_pathways$Direction <- ifelse(all_pathways$NES > 0, "UP", "DOWN")
up_genes <- all_pathways$Description[all_pathways$NES > 0]
down_genes <- all_pathways$Description[all_pathways$NES < 0]
all_pathways<- unique(c(up_genes, down_genes))

pathway_matrix <- data.frame(
	Pathway = all_pathways, 
	UP = as.integer(all_pathways %in% up_genes), 
	DOWN = as.integer(all_pathways %in% down_genes)
)

write.csv(pathway_matrix, "results/EN_NASH_vs_HEALTHY_pathway_presence_absence.csv", row.names = FALSE)

# 2. ORA
log2fc_list <- list()
results_list <- list()
enrich_list <- list()

heat_colors <- colorRampPalette(c("green", "black", "red"))(150)
gene_sets <- clusterProfiler::read.gmt("/media/matteo/SSD5/lucrezia/ref/c5.go.bp.v2024.1.Hs.symbols.gmt")

results_list[["NASH_vs_NORMAL"]] <- res_gene_sel
log2fc_list[["NASH_vs_NORMAL"]] <- res_gene_sel$logFC
names(log2fc_list[["NASH_vs_NORMAL"]]) <- res_gene_sel$SYMBOL

# split genes by direction
up_genes <- res_gene_sel$SYMBOL[res_gene_sel$logFC > 0]
down_genes <- res_gene_sel$SYMBOL[res_gene_sel$logFC < 0]

## check
length(up_genes)
length(down_genes)

# run enrichment
enriched_up <- clusterProfiler::enricher(
	gene = up_genes,
	TERM2GENE = gene_sets
)

enriched_down <- clusterProfiler::enricher(
	gene = down_genes,
	TERM2GENE = gene_sets
)

enriched_all <- clusterProfiler::enricher(
	gene = res_gene_sel$SYMBOL,
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

Up pathtway
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

Down pathway
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
  
png("plots/dotplot_enrich.png", width=1200, height=1200, res=150)
dotplot(enriched_all)
dev.off()

# 3. KEGG
Human.GRCh38.p13.annot <- read.table(
	"/ref/Human.GRCh38.p13.annot.tsv", 
	header=T, 
	sep="\t", 
	fill= T, 
	comment.char ="#"
)

res_gene_entrez <- merge(
	res_gene_sel, 
	Human.GRCh38.p13.annot[, c("Symbol", "GeneID")], 
	by.x ="SYMBOL", 
	by.y = "Symbol", 
	all.x =T
)

res_gene_entrez_filt <- res_gene_entrez[!is.na(res_gene_entrez$GeneID), ] # filter NA values

#https://www.kegg.jp/pathway/hsa04932
#pathway=hsa04932

path <- keggGet("hsa04932") #get the pathway
nafld_genes<- path[[1]]$GENE #extract the genes

nafld_entrez <- nafld_genes[seq(1, length(nafld_genes), 2)] #extract the entrez
nafld_entrez <- sub(";.*", "", nafld_entrez) #remove ;

## overlap to see if the gene I obtained are involved in the pathology
gene_list <- unique(na.omit(res_gene_entrez$GeneID)) #extract the gene obtained
overlap <- intersect(
	gene_list, 
	nafld_entrez
)
length(overlap)

# repeat the analysis with gene symbol instead of entrez to see if there are overlapping genes
res_gene_symbol <- res_gene_sel$SYMBOL
res_gene_symbol <- res_gene_symbol[!is.na(res_gene_symbol)]
nafld_symbols <- nafld_genes[seq(2, length(nafld_genes), 2)]
nafld_symbols <- sub(".*\\[", "", nafld_symbols)  
nafld_symbols <- sub("\\].*", "", nafld_symbols)  
gene_list <- unique(na.omit(res_gene_symbol))
overlap <- intersect(gene_list, nafld_symbols)
length(overlap)
