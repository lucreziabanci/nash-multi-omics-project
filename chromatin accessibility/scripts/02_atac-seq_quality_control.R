#!/bin/bash

# ATAC-seq quality control
# Dataset:PRJNA725028

# 1. Quality control

bamfiles<-list.files(
	pattern="*_sorted.bam$", 
	full.names=TRUE
	)

filelabels<-gsub("_sorted.bam", "", basename(bamfiles))

dir.create("fragment_size", showWarnings=FALSE, recursive=TRUE)

for (i in 1:length(bamfiles)) {
	print(paste("Plotting", filelabels[i], "fragment sizes"))
	tiff(paste0("fragment_size/", filelabels[i], "_frag_size",".tiff"), 600,600)
	ATACseqQC::fragSizeDist(bamfiles[i], filelabels[i])
	dev.off()
}

tags <- c("AS", "XN", "XM", "XO", "XG", "NM", "MD", "YS", "YT")

#check the levels
samtools view -H mapped/SRR14321495.rmDup.sorted.bam | grep "@SQ"

seqlvls <- paste0("chr", c(1:22, "X", "Y"))

#create an object with genome info
seq_obj <- GenomeInfoDb::seqinfo(TxDb.Hsapiens.UCSC.hg19.knownGene::TxDb.Hsapiens.UCSC.hg19.knownGene)

#create GRanges: object with chromosomes and ready to work with functions
which <- as(seq_obj[seqlvls], "GRanges")

#object to make calculation in parallel using 6 core
bpparam <- BiocParallel::MulticoreParam(workers=20)

gal <- list()
for (i in 1:length(bamfiles)) {
    print(paste("Processing", filelabels[i]))
    gal[[i]] <- ATACseqQC::readBamFile(bamfiles[i], tag=tags,
                          which=which, asMates=T, bigFile=T, BPPARAM=bpparam)
    out_bam <- paste0(filelabels[i], "_shift.bam")
    objs <- ATACseqQC::shiftGAlignmentsList(gal[[i]], outbam=paste0(out_bam))
}

ann <- read.delim("ref/genes_ID.annot.txt", header=FALSE, sep="\t")

names(ann)=c("chr", "type", "start", "end","strand", "gene_id")

genes <- GenomicRanges::makeGRangesFromDataFrame(ann) #transform the data.frame in a GRanges

print(paste("Plotting TSS enrichment"))
tsse <- ATACseqQC::TSSEscore(objs, genes) #calculate enrichment

dir.create("TSSEnrichment", showWarnings=FALSE)

png("tss_enrichment.png", 600, 600)

plot(
	100*(-9:10-.5), 
	tsse$values, 
	type="b", 
	xlab="distance to TSS", 
	ylab="aggregate TSS score"
)
dev.off()
