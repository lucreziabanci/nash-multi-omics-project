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
tximport (v1.38.2)
DESeq2 (v1.50.0)
biomaRt (v2.66.2)
clusterProfiler (v4.18.4)
enrichplot (v1.30.5)
ggplot2 (v4.0.0)
KEGGREST (v1.50.0)
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
DESeq2 (v1.50.0)
clusterProfiler (v4.18.4)
enrichplot (v1.30.5)
ggplot2 (v4.0.0)
KEGGREST (v1.50.0)
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

## Epigenomics Analysis
Array DNA methylation differential expression and pathway enrichment analysis in NASH using GPL13534 array.

This repository contains a bioinformatics workflow for the analysis of array data from liver samples:
- Download GEO file;
- Calculate beta and M values;
- Quality control using mean detection p-value and density;
- PCA analysis;
- Differential expression analysis with limma;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

This pipeline has been tested on the following dataset:

- **GEO accession**: GSE48325
- **Organism**: Homo sapiens

### Usage

Overview of the pipeline:

1) `download_beta_and_M_values.R`: download GEO files, apply quality filter and calculate beta and M values.
2) `quality_control_and_PCA.R`: perform quality control and PCA analysis.
3) `differential_analysis_and_volcano_plot.R`: perform differential expression analysis using limma and generate volcano plots.
4) `functional_and_enrichment_analysis.R`: perform GSEA, ORA, and KEGG pathway analyses.

### Dependencies

The following annotation files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:

- `c5.go.bp.v2024.1.Hs.symbols.gmt`
- `Human.GRCh38.p13.annot.tsv`

The following R packages are required:

```
GEOquery (v2.78.0)
limma (v3.66.0)
ggplot2 (v4.0.0)
Biobase (v2.70.0)
minfi (v1.56.0)
IlluminaHumanMethylation450kmanifest (v0.4.0)
lluminaHumanMethylation450kanno.ilmn12.hg19 (v0.6.1)
clusterProfiler (v4.18.4)
enrichplot (v1.30.5)
org.Hs.eg.db (v3.22.0)
dplyr (v1.1.4)
VennDiagram (v1.7.3)
KEGGREST (v1.50.0)
```

### Output
The pipeline generates:
- PCA plots
- Volcano plots
- Differentially expressed gene tables
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

## Proteomics Analysis
LC/MS proteomics differential expression and pathway enrichment analysis in NASH using MaxQuant (v2.7.5.0) data.

This repository contains a bioinformatics workflow for the analysis of array data from liver samples:
- Upload proteinGroups.txt file;
- Differential expression analysis with limma;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

This pipeline has been tested on the following dataset:

- **PRIDE accession**: PXD026717
- **Organism**: Homo sapiens

### Usage

Overview of the pipeline:

1) `data_upload_and_differential_analysis.R`: upload data from MaxQuant analysis, perform differential analysis using limma and generate volcano plot.
2) `functional_analysis.R`: perform ORA, GSEA and KEGG analysis.

### Dependencies

The following annotation files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:

- `c5.go.bp.v2024.1.Hs.symbols.gmt`
- `Human.GRCh38.p13.annot.tsv`

The following R packages are required:

```
limma (v3.66.0)
ggplot2 (v4.0.0)
org.Hs.eg.db (v3.22.0)
clusterProfiler (v4.18.4)
KEGGREST (v1.50.0)
```

### Output
The pipeline generates:
- Volcano plots
- Differentially expressed gene tables
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

## miRNA Analysis
Array miRNA differential expression and pathway enrichment analysis in NASH using GPL16384 array data.

This repository contains a bioinformatics workflow for the analysis of array data from liver samples:
- Download data from GEOquery;
- Differential expression analysis with limma;
- Volcano plot visualization;
- Gene Set Enrichment Analysis (GSEA);
- Over-Representation Analysis (ORA);
- KEGG pathway overlap analysis.

This pipeline has been tested on the following dataset:

- **GEO accession**: GSE59492
- **Organism**: Homo sapiens

### Usage

Overview of the pipeline:

1) `download_and_PCA.R`: download data from GEO and perform PCA.
2) `differential_analysis_and_volcano_plot.R`: perform differential analysis using limma and generate volcano plot.
2) `functional_analysis.R`: perform ORA, GSEA and KEGG analysis.

### Dependencies

The following annotation files should be manually downloaded and placed inside the metadata/ directory before running enrichment analyses:

- `c5.go.bp.v2024.1.Hs.symbols.gmt`
- `Human.GRCh38.p13.annot.tsv`

The following R packages are required:

```
GEOquery (v2.78.0)
Biobase (v2.70.0)
limma (v3.66.0)
umap (v0.2.10.0)
ggplot2 (v4.0.0)
clusterProfiler (v4.18.4)
enrichplot (v1.30.5)
org.Hs.eg.db (v3.22.0)
KEGGREST (v1.50.0)
VennDiagram (v1.7.3)

```

### Output
The pipeline generates:
- PCA plots
- Volcano plots
- Differentially expressed gene tables
- GSEA enrichment plots
- ORA pathway analyses
- KEGG overlap analyses

# Author

Lucrezia Banci, Lorenzo Casbarra, and Matteo Ramazzotti
