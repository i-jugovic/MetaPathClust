# Metabolic Pathway Clustering (MetaPathClust)
**MetaPathClust** is a python script for enzyme sequence and metabolic pathway comparison and clustering using query cover, percent identity, and E-values. Performs PCA, Elbow method, Silhouette analysis, K-means, and hierarchical clustering to group species metabolically, integrates taxonomic grouping, and outputs visualizations with a short summary statistics.

# Overview - MetaPathClust
**Purpose.** Given a semicolon-delimited CSV of enzyme search metrics (per pathway P1–P8: query coverage, % identity, E-value) across species, the script:
- **Loads & cleans** feature columns (coverage, identity, E-value).
- **Transforms** E-values (−log10), standardizes features, and **handles missing values**.
- **Assigns taxonomic groups** by name heuristics for colored plots and special markers.
- Runs **PCA** (full and three feature subsets) with variance plots.
- Performs **cluster model selection** (Elbow via WCSS + Silhouette) and then:
  - **Hierarchical clustering** (Ward) with PCA scatter + dendrogram.
  - **K-means** at the selected k with a **detailed silhouette plot** and PCA view.
- **Highlights** the selected species (in our example *Gluconacetobacter brunescens* Hr-1-5) if present.
- **Writes figures** (PNG/SVG) and prints a concise **summary report** (feature counts, group/cluster distributions, optimal k, silhouette score, explained variance, file list).

# System requirements
## Runtime
- Python ≥ 3.9 (3.10+ recommended)

## Python packages
- `pandas`, `numpy`
- `matplotlib`, `seaborn`
- `scikit-learn` (`sklearn`)
- `scipy`
- `kneed`
- (optional) `warnings` is stdlib

Make sure the packages are installed.

# Input data layout
An example .csv file is uploaded (named `Raw enzyme data.csv`). The delimiter is the semicolon `;`. The first required column is the `Species` column. For pathways there are enzymes from P1 to P8 in our example. `_Cov` is Query Cover in percentages, `ID` is Percent Identity also in percentages, `_E` is E-value in the scientific notation (Note: E-values may include zeros (treated as very significant), scientific notation, or blanks.).

# Outputs
- Figures (PNG + SVG):
  - `pca_full_dataset.`, `pca_p1_p2_p7_p8.`, `pca_p4_p5_p6_p7_p8.`, `pca_p1_p2_p4_p5_p6.`
  - `cluster_optimization_hierarchical.`, `hierarchical_clustering_optimal.`
  - `cluster_optimization_kmeans.`, `silhouette_analysis_detailed.`
- Console summary with optimal k, silhouette score, distributions, and file list.

# Citing

