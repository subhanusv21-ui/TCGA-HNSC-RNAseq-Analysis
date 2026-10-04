# TCGA-HNSC RNA-seq Analysis

## Project Overview

This project analyzes RNA-seq gene expression data from The Cancer Genome Atlas (TCGA) for head and neck squamous cell carcinoma (HNSC).

The aim is to identify genes that are differentially expressed between tumour and normal tissue and investigate the biological processes and pathways associated with these changes.

## Research Question

**Which genes are significantly differentially expressed between HNSC tumour and normal tissue, and what biological processes and pathways do these genes regulate?**

## Dataset

- **Cancer type:** Head and Neck Squamous Cell Carcinoma (TCGA-HNSC)
- **Data type:** RNA-seq gene expression
- **Data source:** The Cancer Genome Atlas (TCGA)
- **Data format:** STAR-Counts
- **Tumour samples:** Primary Tumour
- **Normal samples:** Solid Tissue Normal

## Analysis Workflow

1. Retrieved TCGA-HNSC RNA-seq data.
2. Examined sample types and clinical information.
3. Removed duplicate tumour-normal sample barcodes.
4. Filtered low-expression genes.
5. Performed variance-stabilizing transformation.
6. Visualized sample variation using PCA.
7. Performed differential expression analysis using DESeq2.
8. Generated volcano plots and heatmaps.
9. Performed Gene Ontology (GO) enrichment analysis.
10. Performed pathway enrichment analysis.
11. Examined log2 fold-change patterns.
12. Interpreted the biological processes and pathways associated with differentially expressed genes.

## Tools Used

- R
- RStudio
- Bioconductor
- TCGAbiolinks
- DESeq2
- clusterProfiler
- survival
- survminer
- org.Hs.eg.db
- ggplot2

## Key Findings

Differential expression analysis showed substantial transcriptional differences between HNSC tumour and normal tissue.

Functional enrichment analysis highlighted biological processes related to extracellular matrix organization, cell division and immune-associated functions.

Pathway analysis also highlighted pathways including ECM-receptor interaction and integrin signalling.

These results provide a basis for further investigation of immune-associated transcriptional patterns in HNSC.

## Visualizations

The `figures/` folder contains the main visual outputs from the analysis, including:

- PCA plot
- Volcano plot
- Differentially expressed gene heatmap
- GO enrichment bar plot
- GO enrichment dot plot
- Pathway enrichment plots
- log2 fold-change plot

## Project Status

This is an ongoing bioinformatics portfolio project.

The current stage focuses on tumour-versus-normal differential expression and functional enrichment. The next stage will investigate immune-associated genes and their relationship with patient survival.

## Repository Structure

```text
TCGA-HNSC-RNAseq-Analysis/
│
├── README.md
│
└── figures/
    ├── README.md
    ├── volcano_plot.png
    ├── GO_enrichment_barplot.png
    ├── GO_enrichment_dotplot.png
    ├── DEG_heatmap.png
    ├── pathway_enrichment_1.png
    ├── pathway_enrichment_2.png
    └── log2FC_plot.png
```

## Author

**Subhanu SV**  
MSc Biotechnology  
St. Joseph's University, Bengaluru
