# TCGA-HNSC RNA-seq Analysis

## Project Overview

This project analyzes RNA-seq gene expression data from The Cancer Genome Atlas (TCGA) for head and neck squamous cell carcinoma (HNSC).

The analysis compares tumour and normal tissue to identify differentially expressed genes and investigate the biological processes and pathways associated with these transcriptional changes.

The project was developed as a beginner-to-intermediate bioinformatics portfolio project using R and Bioconductor.

## Research Question

**Which genes are differentially expressed between HNSC tumour and normal tissue, and what biological processes and pathways are associated with these changes?**

## Dataset

- **Cancer type:** Head and Neck Squamous Cell Carcinoma (TCGA-HNSC)
- **Data type:** RNA-seq gene expression
- **Data source:** The Cancer Genome Atlas (TCGA)
- **Data format:** STAR-Counts
- **Tumour samples:** Primary Tumour
- **Normal samples:** Solid Tissue Normal

## Analysis Workflow

1. Retrieved TCGA-HNSC RNA-seq data using `TCGAbiolinks`.
2. Examined sample types and sample information.
3. Removed duplicate tumour-normal sample barcodes.
4. Filtered genes with low expression.
5. Performed variance-stabilizing transformation.
6. Visualized sample-level variation using principal component analysis (PCA).
7. Performed differential expression analysis using DESeq2.
8. Annotated Ensembl gene identifiers with gene symbols and gene descriptions.
9. Identified significantly differentially expressed genes.
10. Performed Gene Ontology (GO) enrichment analysis.
11. Performed KEGG pathway enrichment analysis.
12. Generated a heatmap of differentially expressed genes.
13. Visualized representative gene log2 fold-change values.
14. Saved analysis results as CSV files for reproducibility.

## Differential Expression Results

Using an adjusted p-value cutoff of **< 0.05** and an absolute log2 fold-change cutoff of **≥ 1**, the analysis identified:

- **612 significant differentially expressed genes**
- **274 upregulated genes**
- **338 downregulated genes**

This indicates substantial transcriptional differences between HNSC tumour and normal tissue.

## Functional Enrichment Findings

GO enrichment analysis highlighted biological processes associated with:

- Extracellular matrix organization
- Cell division
- Immune-associated biological processes
- B-cell mediated and adaptive immune functions

Pathway analysis also highlighted pathways including:

- ECM-receptor interaction
- Integrin-related signalling

These findings suggest changes in extracellular matrix organization, cellular proliferation and immune-associated processes within HNSC tumour tissue.

## Visualizations

The `figures/` directory contains the main visual outputs from the analysis:

- GO enrichment bar plot
- GO enrichment dot plot
- Differentially expressed gene heatmap
- Pathway enrichment plots
- Representative gene log2 fold-change plot

## Tools and Packages

### Programming and Analysis

- R
- RStudio
- Bioconductor

### R Packages

- TCGAbiolinks
- DESeq2
- clusterProfiler
- AnnotationDbi
- org.Hs.eg.db
- survival
- survminer
- ggplot2

## Repository Structure

```text
TCGA-HNSC-RNAseq-Analysis/
│
├── README.md
├── TCGA_HNSC_analysis.R
├── TCGA_HNSC_analysis.Rproj
│
├── figures/
│   ├── GO_enrichment_barplot.png
│   ├── GO_enrichment_dotplot.png
│   ├── DEG_heatmap.png
│   ├── pathway_enrichment_1.png
│   ├── pathway_enrichment_2.png
│   └── log2FC_plot.png
│
└── results/
    ├── HNSC_DEG_results.csv
    ├── HNSC_GO_upregulated.csv
    ├── HNSC_GO_downregulated.csv
    ├── HNSC_KEGG_upregulated.csv
    └── HNSC_KEGG_downregulated.csv
```

## Project Status

### Completed

- TCGA-HNSC RNA-seq data retrieval
- Sample processing and filtering
- PCA
- Differential expression analysis
- Gene annotation
- GO enrichment analysis
- KEGG pathway analysis
- Visualization
- Export of analysis results

### Planned Next Step

The next stage of the project will investigate **immune-associated gene expression patterns and their relationship with patient survival** using clinical and survival data from TCGA-HNSC.

## Author

**Subhanu SV**  
MSc Biotechnology  
St. Joseph's University, Bengaluru
