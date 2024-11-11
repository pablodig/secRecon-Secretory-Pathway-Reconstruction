# secRecon: Secretory Pathway Reconstruction

This repository, **secRecon**, contains all the necessary code to generate the figures and analyses presented in the  paper:

This paper introduces secRecon, a comprehensive reconstruction of the mammalian secretory pathway, highlighting its utility in contextualizing omics data and uncovering insights into secretory phenotypes. The repository includes multiple modules to facilitate the reconstruction and analysis of the secretory pathway, and it is structured into several folders, each addressing specific aspects of the network reconstruction and analysis process.

## 1 - Network Reconstruction

This folder, **01 - Network Reconstruction**, contains all the essential used to map data from different databases and datasets into secRecon and the Functional and PPI networkl topology generation.

### Contents

1. **1.1 - Feature\_Extraction.ipynb**: This Jupyter notebook extracts relevant features from different databases and includes it into secRecon

2. **1.2 - Functional\_Ontology\_Network.ipynb**: This notebook is responsible for constructing the **Functional Ontology Network**. It uses the processed gene and process data to build a network that describes relationships between genes based on their biological function and ontology.&#x20;

3. **1.3 - PPI\_Ontology\_Network.ipynb**: The third notebook in the series focuses on integrating **Protein-Protein Interactions (PPI)** with ontology data to create a comprehensive **PPI Ontology Network**.&#x20;

## 2 - CHO vs Plasma

This folder contains notebooks and scripts used to perform various analyses on the CHO vs Plasma multi-omics datasets, aiming to extract meaningful insights using secRecon.

### Contents

1. **2.1 - preprocess\_multiomics.ipynb**: This Jupyter notebook focuses on preprocessing multi-omics data, which includes transcriptomics, proteomics, and other omics data types.

2. **2.2 - correlation\_hclust.Rmd**: This R Markdown file performs hierarchical clustering of correlation data from multi-omics datasets. The clustering helps in identifying patterns across genes and proteins, providing insights into their co-regulation and functional interactions.

3. **2.3 - gsva\_limma.Rmd**: This R Markdown file applies Gene Set Variation Analysis (GSVA) followed by differential analysis using `limma`. The goal is to identify significantly varying pathways across different conditions in the secretory pathway.

4. **2.4 - gsva\_secRecon.Rmd**: This file focuses on applying GSVA specifically to the reconstructed secretory pathway (secRecon). The analysis helps in understanding how specific pathways behave under various experimental conditions.

5. **2.5 - Cho\_vs\_plasma\_network\_analysis.ipynb**: This Jupyter notebook generates and compares network features between CHO cells and plasma cells. It aims to highlight significant differences that could be leveraged to enhance protein production or understand cell-specific characteristics of the secretory pathway.

## 3 - Sec-seq Analysis

This folder, contains scripts used to perform analyses with a Sec-seq dataset using secRecon, such as correlation studies, dimension reduction, and dominance analysis.

### Contents

1. **0a\_extract\_ppi\_network.Rmd**: This R Markdown file extracts the Protein-Protein Interaction (PPI) network data. This serves as the initial step to prepare the dataset for further correlation analysis and network integration.

2. **0b\_Preprocess-DimensionReduction-Pseudotime.ipynb**: This Jupyter notebook preprocesses the omics data and performs dimension reduction and pseudotime analysis. This helps in identifying key states and transitions in the secretory pathway.

3. **1\_score-secRecon-activity.ipynb**: This notebook scores the activity of secRecon under different conditions. It aims to quantify pathway activity and identify conditions that promote or suppress secretory functions.

4. **2\_Dominance-Analysis.ipynb**: This notebook performs dominance analysis to determine the impact of different genes and conditions on the secretory pathway. This type of analysis is useful for identifying major regulators and their influence on pathway outcomes.

5. **3\_Correlation-DEG.ipynb**: This notebook identifies differentially expressed genes (DEGs) and performs correlation analysis. The results help in understanding co-regulation and identifying key gene modules that respond to experimental perturbations.

6. **4\_IgGpopulation\_correlation.Rmd**: This R Markdown file performs correlation analysis focused specifically on IgG populations. It aims to explore relationships within IgG production and how these relate to other functional components in the secretory pathway.

## 4 - secRecon\_TF\_enrichment

This folder contains notebooks and scripts for analyzing transcription factor enrichment in the secRecon components. The aim is to identify transcription factors that may play a regulatory role in the observed pathway dynamics.

### Contents

1. **secRecon\_TF\_enrichment.ipynb**: This Jupyter notebook performs transcription factor enrichment analysis on secRecon. It identifies key transcription factors that may be regulating the secretory pathway and visualizes their potential impact on network activity.

## 5 - Gtex\_correlation

This folder contains scripts and data for characterizing the secretory pathway using GTEx data. The aim is to explore gene expression in various tissues and gain insights into how the secretory pathway is regulated across different biological contexts.

### Contents

1. **0\_gtex\_secrecon\_char.Rmd**: This R Markdown file uses GTEx data to characterize the reconstructed secretory pathway. It performs tissue-specific analysis, identifying differential expression patterns and linking these patterns to secretory pathway activity.

## Usage Instructions

1. Clone this repository.
2. Navigate to the relevant folder for the analysis you wish to perform.
3. Make sure you have the appropriate software installed (e.g., Jupyter Notebook, R) along with the required Python and R packages. You can install the Python dependencies via:
   ```bash
   pip install -r requirements.txt
   ```
4. Run the Jupyter notebooks or R Markdown files in the specified order to ensure proper data processing and analysis.

## Unified Requirements

- Python 3.8+
- R version 4.0+
- Jupyter Notebook
- Required Python packages (install via `pip` or conda): `pandas`, `numpy`, `networkx`, `matplotlib`, `scikit-learn`, `scanpy`, `seaborn`
- Required R packages: `limma`, `GSVA`, `ggplot2`, `dplyr`, `corrplot`, `EnhancedVolcano`

---
