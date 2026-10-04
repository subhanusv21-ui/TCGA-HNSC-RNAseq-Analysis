# ============================================================
# TCGA-HNSC RNA-seq Analysis
# ============================================================
#
# Project:
# Identification of immune-associated transcriptomic
# signatures associated with patient survival in HNSC
#
# Current analysis:
# Tumour vs Normal Differential Gene Expression
# and Functional Enrichment Analysis
#
# Dataset: TCGA-HNSC
# Data type: RNA-seq STAR-Counts
# Author: Subhanu SV
# ============================================================


# ------------------------------------------------------------
# 1. Load required packages
# ------------------------------------------------------------

library(TCGAbiolinks)
library(SummarizedExperiment)
library(DESeq2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(ggplot2)
library(pheatmap)


# ------------------------------------------------------------
# 2. Retrieve TCGA-HNSC RNA-seq data
# ------------------------------------------------------------

query_hnsc <- GDCquery(
  project = "TCGA-HNSC",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query_hnsc)

hnsc_data <- GDCprepare(query_hnsc)


# ------------------------------------------------------------
# 3. Examine sample information
# ------------------------------------------------------------

sample_info <- colData(hnsc_data)

table(sample_info$sample_type)


# ------------------------------------------------------------
# 4. Remove duplicate tumour-normal sample barcodes
# ------------------------------------------------------------

sample_barcodes <- substr(colnames(hnsc_data), 1, 12)

duplicate_barcodes <- duplicated(sample_barcodes) |
  duplicated(sample_barcodes, fromLast = TRUE)

sum(duplicate_barcodes)

hnsc_data <- hnsc_data[, !duplicate_barcodes]

table(colData(hnsc_data)$sample_type)


# ------------------------------------------------------------
# 5. Define tumour and normal conditions
# ------------------------------------------------------------

condition <- ifelse(
  colData(hnsc_data)$sample_type == "Primary Tumor",
  "Tumor",
  "Normal"
)

condition <- factor(
  condition,
  levels = c("Normal", "Tumor")
)

colData(hnsc_data)$condition <- condition

table(colData(hnsc_data)$condition)


# ------------------------------------------------------------
# 6. Create DESeq2 dataset
# ------------------------------------------------------------

dds <- DESeqDataSet(
  hnsc_data,
  design = ~ condition
)

dds


# ------------------------------------------------------------
# 7. Filter low-expression genes
# ------------------------------------------------------------

keep <- rowSums(counts(dds) >= 10) >= 2

dds <- dds[keep, ]

dim(dds)


# ------------------------------------------------------------
# 8. Perform differential expression analysis
# ------------------------------------------------------------

dds <- DESeq(dds)

res <- results(
  dds,
  contrast = c("condition", "Tumor", "Normal")
)

summary(res)


# ------------------------------------------------------------
# 9. Variance-stabilizing transformation
# ------------------------------------------------------------

vsd <- vst(dds, blind = FALSE)


# ------------------------------------------------------------
# 10. Principal Component Analysis
# ------------------------------------------------------------

pca_plot <- plotPCA(
  vsd,
  intgroup = "condition"
) +
  ggtitle("PCA of TCGA-HNSC RNA-seq Samples") +
  theme_classic()

pca_plot


# ------------------------------------------------------------
# 11. Annotate differential expression results
# ------------------------------------------------------------

gene_annot <- AnnotationDbi::select(
  org.Hs.eg.db,
  keys = rownames(res),
  columns = c("SYMBOL", "GENENAME"),
  keytype = "ENSEMBL"
)

res_df <- as.data.frame(res)

res_df$ENSEMBL <- rownames(res_df)

res_annotated <- merge(
  gene_annot,
  res_df,
  by = "ENSEMBL"
)

res_annotated <- res_annotated[
  !duplicated(res_annotated$ENSEMBL),
]


# ------------------------------------------------------------
# 12. Identify differentially expressed genes
# ------------------------------------------------------------

deg <- subset(
  res_annotated,
  padj < 0.05 &
    abs(log2FoldChange) >= 1
)

up_genes <- subset(
  deg,
  log2FoldChange >= 1
)

down_genes <- subset(
  deg,
  log2FoldChange <= -1
)

n_total <- nrow(deg)
n_up <- nrow(up_genes)
n_down <- nrow(down_genes)

n_total
n_up
n_down


# ------------------------------------------------------------
# 13. Rank significant genes
# ------------------------------------------------------------

up_genes <- up_genes[
  order(up_genes$padj),
]

down_genes <- down_genes[
  order(down_genes$padj),
]

head(
  up_genes[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "GENENAME",
      "log2FoldChange",
      "pvalue",
      "padj"
    )
  ],
  20
)

head(
  down_genes[
    ,
    c(
      "ENSEMBL",
      "SYMBOL",
      "GENENAME",
      "log2FoldChange",
      "pvalue",
      "padj"
    )
  ],
  20
)


# ------------------------------------------------------------
# 14. Gene Ontology enrichment - Upregulated genes
# ------------------------------------------------------------

up_symbols <- up_genes$SYMBOL

up_symbols <- up_symbols[
  !is.na(up_symbols) &
    up_symbols != ""
]

up_conversion <- bitr(
  up_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_up <- enrichGO(
  gene = up_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_up_simple <- simplify(
  ego_up,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


# ------------------------------------------------------------
# 15. Gene Ontology enrichment - Downregulated genes
# ------------------------------------------------------------

down_symbols <- down_genes$SYMBOL

down_symbols <- down_symbols[
  !is.na(down_symbols) &
    down_symbols != ""
]

down_conversion <- bitr(
  down_symbols,
  fromType = "SYMBOL",
  toType = "ENTREZID",
  OrgDb = org.Hs.eg.db
)

ego_down <- enrichGO(
  gene = down_conversion$ENTREZID,
  OrgDb = org.Hs.eg.db,
  keyType = "ENTREZID",
  ont = "BP",
  pAdjustMethod = "BH",
  pvalueCutoff = 0.05,
  qvalueCutoff = 0.05,
  readable = TRUE
)

ego_down_simple <- simplify(
  ego_down,
  cutoff = 0.7,
  by = "p.adjust",
  select_fun = min
)


# ------------------------------------------------------------
# 16. KEGG pathway enrichment - Upregulated genes
# ------------------------------------------------------------

ekegg_up <- enrichKEGG(
  gene = up_conversion$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)


# ------------------------------------------------------------
# 17. KEGG pathway enrichment - Downregulated genes
# ------------------------------------------------------------

ekegg_down <- enrichKEGG(
  gene = down_conversion$ENTREZID,
  organism = "hsa",
  pvalueCutoff = 0.05
)


# ------------------------------------------------------------
# 18. Heatmap of top differentially expressed genes
# ------------------------------------------------------------

top_genes <- head(
  order(res$padj),
  20
)

top_gene_ids <- rownames(res)[top_genes]

mat <- assay(vsd)[
  top_gene_ids,
]

annotation_col <- as.data.frame(
  colData(dds)[
    ,
    "condition",
    drop = FALSE
  ]
)

pheatmap(
  mat,
  scale = "row",
  annotation_col = annotation_col,
  show_rownames = TRUE,
  fontsize_row = 8,
  main = "Top 20 Differentially Expressed Genes"
)


# ------------------------------------------------------------
# 19. Representative differentially expressed genes
# ------------------------------------------------------------

key_up_genes <- c(
  "MMP11",
  "MMP9",
  "MMP13",
  "COL1A1",
  "COL3A1",
  "POSTN",
  "LOXL2",
  "MKI67",
  "TOP2A",
  "CDK1",
  "CCNB1",
  "PLK1",
  "AURKA",
  "UBE2C",
  "KIF2C"
)

key_down_genes <- c(
  "MYH11",
  "ACTA2",
  "DES",
  "MYL9",
  "CNN1",
  "TPM2",
  "PLN",
  "CCL28",
  "CST4",
  "MPO",
  "GPD1",
  "PLIN1"
)

key_genes <- c(
  key_up_genes,
  key_down_genes
)

key_table <- res_annotated[
  res_annotated$SYMBOL %in% key_genes,
  c(
    "SYMBOL",
    "GENENAME",
    "log2FoldChange",
    "padj"
  )
]

key_table <- key_table[
  match(key_genes, key_table$SYMBOL),
]

key_table


# ------------------------------------------------------------
# 20. Representative gene log2 fold-change plot
# ------------------------------------------------------------

ggplot(
  key_table,
  aes(
    x = reorder(SYMBOL, log2FoldChange),
    y = log2FoldChange
  )
) +
  geom_col() +
  coord_flip() +
  geom_hline(
    yintercept = 0,
    linetype = "dashed"
  ) +
  labs(
    title = "Representative Differentially Expressed Genes in HNSC",
    x = "Gene",
    y = "log2 Fold Change (Tumor vs Normal)"
  ) +
  theme_classic()


# ------------------------------------------------------------
# 21. Save important results
# ------------------------------------------------------------

dir.create(
  "results",
  showWarnings = FALSE
)

dir.create(
  "figures",
  showWarnings = FALSE
)

write.csv(
  as.data.frame(res_annotated),
  "results/HNSC_differential_expression_results.csv",
  row.names = FALSE
)

write.csv(
  key_table,
  "results/HNSC_representative_genes.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_up_simple),
  "results/HNSC_GO_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ego_down_simple),
  "results/HNSC_GO_downregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg_up),
  "results/HNSC_KEGG_upregulated.csv",
  row.names = FALSE
)

write.csv(
  as.data.frame(ekegg_down),
  "results/HNSC_KEGG_downregulated.csv",
  row.names = FALSE
)


# ------------------------------------------------------------
# 22. Save PCA plot
# ------------------------------------------------------------

ggsave(
  "figures/PCA_plot.png",
  plot = pca_plot,
  width = 7,
  height = 5,
  dpi = 300
)


# ============================================================
# End of current analysis
# ============================================================
#
# Next stage:
# Identify immune-associated genes and investigate their
# relationship with patient survival.
# ============================================================