# NASH Multi-omics Project
NASH multi-omics project based on the analysis of transcriptomics, chromatin accessibility, epigenomics, proteomics and miRNA datasets

## Transcriptomics Analysis
RNA-seq differential expression and pathway enrichment analysis in NASH.

This repository contains a bioinformatics workflow for the analysis of RNA-seq data from liver samples:
- SRA download and FASTQ file generation;
- Transcript quantification using Salmon;
- PCA analysis;
- Differential expression analysis with DESeq2;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

This pipeline has been tested on the following dataset:

- **GEO accession**: GSE126848
- **Organism**: Homo sapiens

### Usage

Overview of the pipeline:

1) `RNA-seq_download_and_quantification.sh`: download SRA files and perform transcript quantification using Salmon.
2) `RNA-seq_quantification_and_pca.R`: import Salmon quantification files and perform PCA analysis.
3) `RNA-seq_differential_analysis_and_volcano.R`: perform differential expression analysis using DESeq2 and generate volcano plots.
4) `RNA-seq_enrichement_and_kegg_analysis.R`: perform GSEA, ORA, and KEGG pathway analyses.

### Dependencies

The following annotation files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:

- `c5.go.bp.v2024.1.Hs.symbols.gmt`
- `Human.GRCh38.p13.annot.tsv`

The following R packages are required:

```
tximport (v 1.38.2)
DESeq2 (v 1.50.0)
biomaRt (v 2.66.2)
clusterProfiler (4.18.4)
enrichplot (v 1.30.5)
ggplot2 (v 4.0.0)
KEGGREST (1.50.0)
```

### Output
The pipeline generates:
- PCA plots
- Volcano plots
- Differentially expressed gene tables
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

## Chromatin Accessibility Analysis
ATAC-seq differential expression and pathway enrichment analysis in NASH.

This repository contains a bioinformatics workflow for the analysis of ATAC-seq data from liver samples:
- SRA download and FASTQ and BAM file generation;
- Quality control;
- Peak calling;
- Differential expression analysis with DESeq2;
- PCA analysis;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

This pipeline has been tested on the following dataset:
- **GEO accession**: PRJNA725028
- **Organism**: Homo sapiens

Usage
Overview of the pipeline:
1) `ATAC-seq_preprocessing_and_mapping.sh`: download SRA files, map with Bowtie2 and process BAM files.
2) `ATAC-seq_quality_control.R`: import processed BAM files and perform quality control.
3) `ATAC-seq_peak_calling.sh`: call the peak.
4) `ATAC-seq_differential_analysis_pca_and_volcano_plot.R`: perform differential expression analysis using DESeq2 and generate pca and volcano plots.
5) `ATAC-seq_enrichment_and_kegg_analysis.R`: perform GSEA, ORA, and KEGG pathway analyses.

Dependencies
The following annotation files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:

- `TxDb.Hsapiens.UCSC.hg19.knownGene`
- `genes_ID.annot.txt`
- `EnsDb.Hsapiens.v75`
- `c5.go.bp.v2024.1.Hs.symbols.gmt`
- `Human.GRCh38.p13.annot.tsv`

The following R packages are required:

```
ATACseqQC (v1.34.0)
rtracklayer (v1.70.1)
GenomicRanges (v1.62.0)
ChIPpeakAnno (v3.44.0)
csaw (v1.44.0)
DESeq2 (v 1.50.0)
clusterProfiler (4.18.4)
enrichplot (v 1.30.5)
ggplot2 (v 4.0.0)
KEGGREST (1.50.0)
```

### Output
The pipeline generates:
- PCA plots
- Volcano plots
- Differentially expressed gene tables
- TSS enrichment
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

# Author

Lucrezia Banci, Lorenzo Casbarra, and Matteo Ramazzotti
