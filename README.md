<!DOCTYPE html>
<html>
<head>
    <title>01 Preprocessed Data - README</title>
</head>
<body>
    
    <h1>01 Preprocessed Data - README</h1>

    <p>This notebook is primarily concerned with the preprocessing of single-cell RNA sequencing (scRNA-seq) and antibody-derived tags (ADTs) datasets from three different experiments (<code>exp100</code>, <code>exp105</code>, and <code>exp106</code>). The preprocessing includes quality control (QC), filtering, and exploratory data analysis.</p>

    <h2>Dependencies</h2>
    <ul>
        <li>Scanpy</li>
        <li>Matplotlib</li>
        <li>NumPy</li>
    </ul>

    <h2>Summary of Steps</h2>

    <h3>Data Loading</h3>
    <p>Datasets from three experiments are loaded into Scanpy's AnnData objects for further analysis.</p>

    <h3>Quality Control and Filtering</h3>
    <ol>
        <li><strong>Filter cells based on gene count:</strong> Cells that have less than a certain number of genes are filtered out to remove low-quality cells.
            <ul>
                <li><code>exp100</code>: Minimum 700 genes</li>
                <li><code>exp105</code>: Minimum 500 genes</li>
                <li><code>exp106</code>: Minimum 1000 genes</li>
            </ul>
        </li>
        <li><strong>Annotate mitochondrial genes:</strong> Mitochondrial genes are annotated as 'mt'.</li>
        <li><strong>Calculate QC metrics:</strong> Additional quality control metrics such as <code>total_counts</code> and <code>pct_counts_mt</code> (percent of counts in mitochondrial genes) are calculated for each dataset.</li>
    </ol>

    <h3>Data Visualization</h3>
    <ol>
        <li><strong>Violin Plots:</strong> Plots are generated to visualize the distribution of QC metrics such as the number of genes, total counts, and percentage of mitochondrial genes.
            <p><img src="path/to/your/violin_plot_figure.png" alt="Violin Plots"></p>
        </li>
        <li><strong>Scatter Plots:</strong> Scatter plots are used to explore relationships between various metrics, including <code>total_counts</code> and <code>pct_counts_mt</code>.
            <p><img src="path/to/your/scatter_plot_figure.png" alt="Scatter Plots"></p>
        </li>
    </ol>

    <h3>Further Filtering Based on QC Metrics</h3>
    <p>Cells are further filtered based on QC metrics. For example, cells with more than 50,000 total counts or 10% counts in mitochondrial genes are removed in <code>exp100</code>.</p>

    <h3>Copying Raw Counts</h3>
    <p>Original count data is copied to a new layer named "counts" for future reference.</p>

    <h3>Cleaning Up Doublets</h3>
    <p>Manual doublet cleanup is performed based on mutual exclusivity of heavy chain isotypes. Thresholds are set for different immunoglobulin genes, and cells violating these thresholds are removed.</p>

    <h3>Saving Processed Data</h3>
    <p>Finally, the cleaned-up dataset is saved to disk as <code>exp106.h5</code>.</p>

    <h3>Pearson Correlations</h3>
    <p>Pearson correlation is calculated between different metrics such as <code>IGHG1</code> vs <code>IgG_ADT</code> and <code>IGHM</code> vs <code>IgG_ADT</code> for exploratory analysis.</p>

</body>
</html>