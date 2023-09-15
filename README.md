# Single-cell RNA Sequencing (scRNA-seq) and Antibody-derived Tags (ADTs) Data Analysis Pipeline

This repository contains two Jupyter Notebooks that detail the workflow for analyzing scRNA-seq and ADTs datasets. The data originates from three separate experiments: **exp100**, **exp105**, and **exp106**.

## Common Dependencies

- **Scanpy**: For data manipulation and analysis
- **Matplotlib**: For data visualization
- **NumPy**: For numerical operations
- **Seaborn**: For advanced data visualization

---

## 01 Preprocessed Data

The primary aim of this notebook is to preprocess the scRNA-seq and ADTs datasets. The steps are designed to ensure that the data is clean and of high quality, ready for downstream analyses.

### Summary of Steps

#### Data Loading

- Datasets from **exp100**, **exp105**, and **exp106** are loaded into Scanpy's AnnData objects for efficient manipulation and visualization.

#### Quality Control and Filtering

1. **Filter cells based on gene count**: To eliminate low-quality cells, a gene count threshold is applied to each dataset.
    - **exp100**: Minimum 700 genes
    - **exp105**: Minimum 500 genes
    - **exp106**: Minimum 1000 genes
2. **Annotate mitochondrial genes**: These are marked as 'mt' for easy identification and further analyses.
3. **Calculate QC metrics**: Metrics like `total_counts` and `pct_counts_mt` are calculated. These are crucial for ensuring the quality of the cells in the dataset.
4. **Further filtering**: Cells exceeding certain QC metric thresholds are filtered out. For instance, in **exp100**, cells with more than 50,000 total counts or with 10% mitochondrial counts are removed.

#### Data Visualization

- **Violin Plots**: These help in understanding the distribution of QC metrics like gene counts, total counts, and mitochondrial percentage.
- **Scatter Plots**: Used to explore relationships between various QC metrics, enabling more informed filtering decisions.

#### Additional Steps

- **Copying raw counts**: The original count data is saved in a new layer named "counts" for backtracking or additional analyses.
- **Manual doublet cleanup**: Doublets or cell multiplets are identified and removed based on the mutually exclusive expression of heavy chain isotypes.
- **Saving processed data**: The final AnnData object is saved as `exp106.h5` for future use.
- **Pearson Correlations**: Correlations between different metrics like **IGHG1 vs IgG_ADT** and **IGHM vs IgG_ADT** are calculated for exploratory analysis.

---

## 02 Immunoglobulin Annotation and Visualization

The focus of this notebook is to annotate and visualize immunoglobulin (Ig) genes and classes in the datasets.

### Summary of Steps

#### Preparing Data for Annotation

1. **Load preprocessed dataset**: The cleaned dataset from the previous notebook is imported.
2. **Logarithmic scaling of data**: Data is log-scaled to make it more conducive for machine learning algorithms.

#### Immunoglobulin (Ig) Annotation

1. **Heavy Chain Annotation**: The dataset is annotated for the presence of different heavy chains like IgG, IgA, etc.
2. **Light Chain Annotation**: Similar to heavy chains, light chains are also annotated.

#### Labeling

- **Assigning Ig labels**: Cells are labeled based on the presence of specific Ig classes and light and heavy chains.

#### Visualization

1. **UMAP Plotting**: UMAP plots are generated to visualize the clusters of cells based on Ig classes.
2. **Scatter Plots of Ig Subclasses**: These are used to analyze the distribution of different Ig subclasses across cells.

#### Subsetting Data

- **Filtering**: Cells can be filtered out based on their expression of specific Ig classes or subclasses or under certain conditions like `negA` or `negM`.

#### Further Explorations

1. **Data log scaling**: Additional log scaling is performed for more advanced analyses.
2. **Density Plots**: These are used to understand the distribution of cells across different Ig classes.

#### Merged Blocks for Code Efficiency

- A combined block for Ig annotation steps
- A combined block for visualization steps

---

By modularizing each notebook into distinct sections, the analysis becomes more interpretable and manageable, aiding in debugging and modification.


