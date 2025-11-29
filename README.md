# Understudied-genes-analysis
#### Summary

This project analyzes candidate gene expression across single-cell RNA-seq datasets. It includes two complementary workflows: a detailed comparison of two tumor-associated scRNA-seq datasets, and a scalable pipeline capable of integrating ten or more datasets. Using Seurat, the analysis identifies cell types, detects cluster-specific marker genes, compares tumor microenvironments across studies, and evaluates how candidate genes behave across diverse cellular contexts. The multi-dataset workflow automates preprocessing, batch integration, visualization, and pseudo-bulk summarization, enabling high-throughput exploration of hundreds of genes across large scRNA-seq collections.

## scRNA-seq Candidate Gene Analysis
### Part 1: Comparative Analysis of Selected Datasets
This repository demonstrates analysis and comparison of two scRNA-seq datasets, GSE72056 and GSE115978, both profiling tumor and immune/microenvironment cell populations. The workflow includes:

- Subsetting and preprocessing individual Seurat objects.
- Identification of cluster-specific marker genes.
- Cell type annotation based on top markers.
- Comparison of shared and unique cell populations between datasets.
This approach highlights similarities and differences in tumor microenvironments and enables systematic evaluation of candidate gene expression.

### Part 2: Scaling Analysis Across ≥10 scRNA-seq Datasets 
#### (candidate_gene_analysis_multi_scRNAseq_datasets.R)
To handle 10 or more scRNA-seq datasets efficiently:

- Store datasets in a structured directory with consistent file naming.
- Preprocess each dataset individually (normalization, variable feature selection).
- Integrate datasets in batches using Seurat’s FindIntegrationAnchors to reduce memory usage.
- Merge all integrated datasets into a single Seurat object for downstream visualization, clustering, and candidate gene analysis.
- Automate candidate gene expression plots and generate pseudo-bulk summaries to compare hundreds of genes across multiple datasets.
