01 Preprocessed Data - README

This notebook is primarily concerned with the preprocessing of single-cell RNA sequencing (scRNA-seq) and antibody-derived tags (ADTs) datasets from three different experiments (exp100, exp105, and exp106). The preprocessing includes quality control (QC), filtering, and exploratory data analysis.

Dependencies

Scanpy
Matplotlib
NumPy
Summary of Steps

Data Loading
Datasets from three experiments are loaded into Scanpy's AnnData objects for further analysis.

Quality Control and Filtering
Filter cells based on gene count: Cells that have less than a certain number of genes are filtered out to remove low-quality cells.
exp100: Minimum 700 genes
exp105: Minimum 500 genes
exp106: Minimum 1000 genes
Annotate mitochondrial genes: Mitochondrial genes are annotated as 'mt'.
Calculate QC metrics: Additional quality control metrics such as total_counts and pct_counts_mt (percent of counts in mitochondrial genes) are calculated for each dataset.
Data Visualization
Violin Plots: Plots are generated to visualize the distribution of QC metrics such as the number of genes, total counts, and percentage of mitochondrial genes.
Scatter Plots: Scatter plots are used to explore relationships between various metrics, including total_counts and pct_counts_mt.
Further Filtering Based on QC Metrics
Cells are further filtered based on QC metrics. For example, cells with more than 50,000 total counts or 10% counts in mitochondrial genes are removed in exp100.

Copying Raw Counts
Original count data is copied to a new layer named "counts" for future reference.

Cleaning Up Doublets
Manual doublet cleanup is performed based on mutual exclusivity of heavy chain isotypes. Thresholds are set for different immunoglobulin genes, and cells violating these thresholds are removed.

Saving Processed Data
Finally, the cleaned-up dataset is saved to disk as exp106.h5.

Pearson Correlations
Pearson correlation is calculated between different metrics such as IGHG1 vs IgG_ADT and IGHM vs IgG_ADT for exploratory analysis.

