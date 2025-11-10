# Metabolic Pathway Clustering (MetaPathClust)
**MetaPathClust** is a python script for enzyme sequence and metabolic pathway comparison and clustering using query cover, percent identity, and E-values. Performs PCA, Elbow method, Silhouette analysis, K-means, and hierarchical clustering to group species metabolically, integrates taxonomic grouping, and outputs visualizations with a short summary statistics.

# Overview - MetaPathClust
**Purpose:** Given a semicolon-delimited CSV of enzyme search metrics (per pathway P1–P8: query coverage, % identity, E-value) across species, the script:
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
An example .csv file is uploaded (named `Raw enzyme data.csv`). The delimiter is the semicolon `;`. The first required column is the `Species` column. For pathways there are enzymes from P1 to P8 in our example. `_Cov` is Query Coverage in percentages, `ID` is Percent Identity also in percentages, `_E` is E-value in the scientific notation (Note: E-values may include zeros (treated as very significant), scientific notation, or blanks.).

# Outputs
- Figures (PNG + SVG):
  - `pca_full_dataset.`, `pca_p1_p2_p7_p8.`, `pca_p4_p5_p6_p7_p8.`, `pca_p1_p2_p4_p5_p6.`
  - `cluster_optimization_hierarchical.`, `hierarchical_clustering_optimal.`
  - `cluster_optimization_kmeans.`, `silhouette_analysis_detailed.`
- Console summary with optimal k, silhouette score, distributions, and file list.

# Detailed explanatory notes on methods used in the script
This Python script performs multivariate statistical analysis on enzyme data from various bacterial and plant species, incorporating coverage, identity percentages, and E-values from protein homology searches (P1-P8 proteins). 

Firstly, it preprocesses data. The code begins by loading enzyme data from CSV file containing the data from the previously described table. The scientific notation of E-value is converted to float. E-value of 0 is the most significant, presenting the best match, while E-value of 1 is the least significant, presenting the worst match. The zero E-values are replaced with 1×10<sup>-300</sup> to enable -log<sub>10</sub> transformation, in which lower E-values become higher transformed values. Missing E-values are filled with 1. After that the program uses StandardScaler for zero mean and unit variance.

Secondly, the script automatically assigns species to taxonomic groups based on genus names and uses special markers to identify key species of interest, previously mentioned.

Thirdly, analysis takes its place. The program performs Principal Component Analysis (PCA) on multiple feature subsets: full dataset which includes all data, subset of P1, P2, P7 and P8 for betalain biosynthesis, subset of P4, P5, P6, P7 and P8 for eumelanin biosynthesis, and subset of P1, P2, P4, P5 and P6 for the intersection of betalain and melanins biosynthesis. The data is visualized with scatter plots colored by taxonomic groups; special markers are the highlighted species of interest and the explained variance plots show PC contribution.

Finally, to determine the metabolic clusters regardless of data taxonomy we used two complementary methods to find optimal cluster numbers – Elbow method and Silhouette analysis. Elbow method focuses on compactness (minimizing With-in-Cluster Sum of Squares) but does not consider how distinct clusters are from each other, while Silhouette analysis balances compactness and separation to find clusters that are both tight and far apart. Therefore, in case of different optimal k-number of clusters determined by each method, the result from Silhouette analysis is preferred. With this known optimal number of clusters we performed hierarchical clustering with optimal cut-off line and K-means clustering with detailed silhouette analysis with pre-sample coefficients and cluster centroids.

For statistical validation, we employed multiple complementary metrics. The silhouette coefficient (ranging from –1 to +1) was used to assess cluster cohesion and separation, providing a measure of how well-defined the identified clusters were. The explained variance ratios from the PCA were examined to quantify the degree of in-formation retained in the reduced-dimensionality space. Additionally, cross-method validation was performed by comparing the optimal number of clusters (k) suggested by the elbow method with that indicated by the silhouette analysis, ensuring robustness in cluster selection.

The output is automatically generated in high-resolution PNG (300 DPI) and vec-tor SVG formats. The plots are separated for each analysis type, and the script gives a comprehensive summary statistics.

# Usage
This methodology is suited for or can be partly used for:
- **comparative enzymology** (identifying metabolic similarities across species),
- **phylogenetic analysis** (comparing taxonomic vs. functional clustering),
- **biomarker discovery** (finding enzyme profiles that distinguish species groups),
- **evolutionary studies** (understanding enzyme conservation patterns) and
- **synthetic biology** (finding optimal enzyme sources for pathway engineering).

# Citing
Karničnik, B.; Jugović, I.; Accetto, T.; Fanedl, L.; Avguštin, G.; Trček, J. *Gluconacetobacter brunescens* sp. nov., a Novel Acetic Acid Bacterium Isolated from Pear Vinegar, Producing a Water-Soluble Brown Pigment. Microorganisms 2025, 13, x. https://doi.org/10.3390/xxxxx

