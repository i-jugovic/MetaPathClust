"""
MetaPathClust

Python script for enzyme sequence and metabolic pathway comparison and clustering using query cover, percent identity, and E-values.
Performs PCA, Elbow method, Silhouette analysis, K-means, and hierarchical clustering to group species metabolically, integrates taxonomic grouping,
and outputs visualizations with a short summary statistics.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from sklearn.cluster import KMeans
from sklearn.metrics import silhouette_score, silhouette_samples
from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from scipy.spatial.distance import pdist
from kneed import KneeLocator
import warnings
warnings.filterwarnings('ignore')

# Set up plotting parameters
plt.rcParams['figure.figsize'] = (12, 8)
plt.rcParams['font.size'] = 10
sns.set_style("whitegrid")

# Load and prepare data
def load_and_prepare_data():
    """Load the CSV file and prepare data for analysis"""
    # Read the CSV file
    data = pd.read_csv('Raw enzyme data.csv', delimiter=';')

    # Clean column names - remove spaces and special characters
    data.columns = data.columns.str.strip()

    # Extract coverage, identity percentage, and E-value columns for analysis
    # P1_Cov (%), P1_ID (%), P1_E, P2_Cov (%), P2_ID (%), P2_E-, etc.
    feature_cols = []
    for i in range(1, 9):  # P1 to P8
        cov_col = f'P{i}_Cov (%)'
        id_col = f'P{i}_ID (%)'
        eval_col = f'P{i}_E'

        if cov_col in data.columns:
            feature_cols.append(cov_col)
        if id_col in data.columns:
            feature_cols.append(id_col)
        if eval_col in data.columns:
            feature_cols.append(eval_col)

    # Create feature matrix
    X = data[feature_cols].copy()

    # Handle missing values and convert to numeric
    for col in X.columns:
        if '_E' in col:
            # Convert E-values from scientific notation to float
            X[col] = pd.to_numeric(X[col], errors='coerce')
            # E-value of 0 is the BEST (most significant), 1 is worst
            # Replace 0 E-values with a very small number to enable log transformation
            X[col] = X[col].replace(0, 1e-300)
            # Take negative log10 of E-values (better E-values become higher values)
            X[col] = -np.log10(X[col])
            # Handle any remaining infinite values
            X[col] = X[col].replace([np.inf, -np.inf], 300)  # Cap at reasonable value
        else:
            X[col] = pd.to_numeric(X[col], errors='coerce')

    # Fill remaining missing values with appropriate defaults
    for col in X.columns:
        if 'E-value' in col:
            # Missing E-values get value 1 (worst), which becomes 0 after -log10 transformation
            X[col] = X[col].fillna(0)  # 0 means worst significance after transformation
        else:
            X[col] = X[col].fillna(0)  # 0 for coverage and identity

    return data, X, feature_cols

# Define organism groups
def assign_groups(species_names):
    """Assign each species to its taxonomic group"""
    groups = []
    special_markers = []

    for species in species_names:
        species_lower = species.lower()

        # Assign groups based on genus/species names
        if 'novacetimonas' in species_lower:
            groups.append('Novacetimonas')
            if 'hansenii' in species_lower:
                special_markers.append('1')
            else:
                special_markers.append('')
        elif 'komagataeibacter' in species_lower:
            groups.append('Komagataeibacter')
            if 'xylinus' in species_lower:
                special_markers.append('2')
            else:
                special_markers.append('')
        elif 'gluconacetobacter' in species_lower:
            groups.append('Gluconacetobacter')
            if 'brunescens' in species_lower:
                special_markers.append('3')
            elif 'tumulicola' in species_lower:
                special_markers.append('4')
            elif 'tumulisoli' in species_lower:
                special_markers.append('5')
            elif 'diazotrophicus' in species_lower:
                special_markers.append('6')
            elif 'liquefaciens' in species_lower:
                special_markers.append('7')
            else:
                special_markers.append('')
        elif 'streptomyces' in species_lower:
            groups.append('Streptomyces')
            if 'scabiei' in species_lower:
                special_markers.append('8')
            else:
                special_markers.append('')
        elif 'pseudomonas' in species_lower:
            groups.append('Pseudomonas')
            if 'putida' in species_lower:
                special_markers.append('9')
            elif 'aeruginosa' in species_lower:
                special_markers.append('10')
            else:
                special_markers.append('')
        elif 'bacillus' in species_lower:
            groups.append('Bacillus')
            special_markers.append('')
        elif 'escherichia' in species_lower:
            groups.append('Escherichia coli')
            special_markers.append('')
        elif 'acinetobacter' in species_lower:
            groups.append('Acinetobacter')
            if 'baumannii' in species_lower:
                special_markers.append('11')
            else:
                special_markers.append('')
        elif 'silene' in species_lower:
            groups.append('Silene latifolia')
            special_markers.append('')
        elif 'spinacia' in species_lower:
            groups.append('Spinacia oleracea')
            special_markers.append('')
        elif 'beta vulgaris' in species_lower:
            groups.append('Beta vulgaris')
            special_markers.append('')
        elif 'mirabilis' in species_lower:
            groups.append('Mirabilis jalapa')
            special_markers.append('')
        elif 'portulaca' in species_lower:
            groups.append('Portulaca')
            special_markers.append('')
        elif 'amaranthus' in species_lower:
            groups.append('Amaranthus')
            special_markers.append('')
        elif 'bougainvillea' in species_lower:
            groups.append('Bougainvillea')
            special_markers.append('')
        else:
            groups.append('Other')
            special_markers.append('')

    return groups, special_markers

# Optimal cluster determination functions
def find_optimal_clusters_elbow(X_scaled, max_clusters=15):
    """Find optimal number of clusters using elbow method"""
    inertias = []
    K_range = range(2, min(max_clusters + 1, len(X_scaled)))

    for k in K_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        kmeans.fit(X_scaled)
        inertias.append(kmeans.inertia_)

    # Use KneeLocator to find the elbow
    try:
        knee_locator = KneeLocator(K_range, inertias, curve='convex', direction='decreasing')
        optimal_k_elbow = knee_locator.elbow
    except:
        # Fallback: find the point with maximum second derivative
        if len(inertias) >= 3:
            second_derivatives = []
            for i in range(1, len(inertias) - 1):
                second_deriv = inertias[i-1] - 2*inertias[i] + inertias[i+1]
                second_derivatives.append(second_deriv)
            optimal_k_elbow = K_range[np.argmax(second_derivatives) + 1]
        else:
            optimal_k_elbow = 3  # Default fallback

    return optimal_k_elbow, K_range, inertias

def find_optimal_clusters_silhouette(X_scaled, max_clusters=15):
    """Find optimal number of clusters using silhouette analysis"""
    silhouette_scores = []
    K_range = range(2, min(max_clusters + 1, len(X_scaled)))

    for k in K_range:
        kmeans = KMeans(n_clusters=k, random_state=42, n_init=10)
        cluster_labels = kmeans.fit_predict(X_scaled)
        silhouette_avg = silhouette_score(X_scaled, cluster_labels)
        silhouette_scores.append(silhouette_avg)

    optimal_k_silhouette = K_range[np.argmax(silhouette_scores)]

    return optimal_k_silhouette, K_range, silhouette_scores

def plot_cluster_optimization(X_scaled, title_suffix=""):
    """Plot both elbow method and silhouette analysis results"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Elbow method
    optimal_k_elbow, K_range_elbow, inertias = find_optimal_clusters_elbow(X_scaled)
    ax1.plot(K_range_elbow, inertias, 'bo-', linewidth=2, markersize=8)
    ax1.axvline(x=optimal_k_elbow, color='red', linestyle='--', linewidth=2,
                label=f'Optimal k = {optimal_k_elbow}')
    ax1.set_xlabel('Number of Clusters (k)')
    ax1.set_ylabel('Within-Cluster Sum of Squares (WCSS)')
    ax1.set_title(f'Elbow Method')
    ax1.legend()
    ax1.grid(True, alpha=0.3)

    # Silhouette analysis
    optimal_k_silhouette, K_range_sil, silhouette_scores = find_optimal_clusters_silhouette(X_scaled)
    ax2.plot(K_range_sil, silhouette_scores, 'go-', linewidth=2, markersize=8)
    ax2.axvline(x=optimal_k_silhouette, color='red', linestyle='--', linewidth=2,
                label=f'Optimal k = {optimal_k_silhouette}')
    ax2.set_xlabel('Number of Clusters (k)')
    ax2.set_ylabel('Average Silhouette Score')
    ax2.set_title(f'Silhouette Analysis')
    ax2.legend()
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, optimal_k_elbow, optimal_k_silhouette

def plot_silhouette_detailed(X_scaled, n_clusters, title_suffix=""):
    """Create detailed silhouette plot for a specific number of clusters"""
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(16, 6))

    # Perform K-means clustering
    kmeans = KMeans(n_clusters=n_clusters, random_state=42, n_init=10)
    cluster_labels = kmeans.fit_predict(X_scaled)

    # Calculate silhouette scores
    silhouette_avg = silhouette_score(X_scaled, cluster_labels)
    sample_silhouette_values = silhouette_samples(X_scaled, cluster_labels)

    # Silhouette plot
    y_lower = 10
    for i in range(n_clusters):
        ith_cluster_silhouette_values = sample_silhouette_values[cluster_labels == i]
        ith_cluster_silhouette_values.sort()

        size_cluster_i = ith_cluster_silhouette_values.shape[0]
        y_upper = y_lower + size_cluster_i

        color = plt.cm.nipy_spectral(float(i) / n_clusters)
        ax1.fill_betweenx(np.arange(y_lower, y_upper),
                         0, ith_cluster_silhouette_values,
                         facecolor=color, edgecolor=color, alpha=0.7)

        ax1.text(-0.05, y_lower + 0.5 * size_cluster_i, str(i+1))
        y_lower = y_upper + 10

    ax1.set_xlabel('Silhouette Coefficient Values')
    ax1.set_ylabel('Metabolic Cluster Label')
    ax1.set_title(f'Silhouette Plot (k={n_clusters})')
    ax1.axvline(x=silhouette_avg, color="red", linestyle="--",
                label=f'Average Score = {silhouette_avg:.3f}')
    ax1.legend()

    # PCA visualization of clusters
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_scaled)

    colors = plt.cm.nipy_spectral(np.linspace(0, 1, n_clusters))
    for i in range(n_clusters):
        mask = cluster_labels == i
        ax2.scatter(X_pca[mask, 0], X_pca[mask, 1],
                   c=[colors[i]], label=f'Metabolic Cluster {i+1}', alpha=0.7, s=60)

    # Plot cluster centers
    centers_pca = pca.transform(kmeans.cluster_centers_)
    ax2.scatter(centers_pca[:, 0], centers_pca[:, 1],
               c='red', marker='x', s=200, linewidths=3, label='Centroids')

    ax2.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax2.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax2.set_title(f'K-means Clustering (k={n_clusters})')
    ax2.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, cluster_labels, silhouette_avg

# PCA analysis function
def perform_pca_analysis(X, groups, special_markers, species_names, feature_subset=None, title_suffix=""):
    """Perform PCA analysis and create visualization"""

    if feature_subset is not None:
        # Select specific features
        X_subset = X[feature_subset].copy()
        title_suffix = f" - {title_suffix}"
    else:
        X_subset = X.copy()

    # Standardize the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X_subset)

    # Perform PCA
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)

    # Create color map for groups
    unique_groups = list(set(groups))
    colors = plt.cm.tab20(np.linspace(0, 1, len(unique_groups)))
    group_colors = dict(zip(unique_groups, colors))

    # Create the plot
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # PCA scatter plot
    for group in unique_groups:
        mask = np.array(groups) == group
        ax1.scatter(X_pca[mask, 0], X_pca[mask, 1],
                   c=[group_colors[group]], label=group, alpha=0.7, s=60)

    # Highlight special markers
    for i, marker in enumerate(special_markers):
        if marker:
            ax1.scatter(X_pca[i, 0], X_pca[i, 1],
                       c='red', s=150, marker='*',
                       edgecolors='black', linewidth=1)
            ax1.annotate(marker, (X_pca[i, 0], X_pca[i, 1]),
                        xytext=(0, 10), textcoords='offset points',
                        fontsize=8, ha='center')

    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax1.set_title(f'PCA Analysis of Enzyme Data{title_suffix}')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left', prop={'style': 'italic'})
    ax1.grid(True, alpha=0.3)

    # Explained variance plot
    cumvar = np.cumsum(pca.explained_variance_ratio_)
    ax2.bar(range(1, len(cumvar)+1), pca.explained_variance_ratio_, alpha=0.7)
    ax2.plot(range(1, len(cumvar)+1), cumvar, 'ro-', linewidth=2)
    ax2.set_xlabel('Principal Component')
    ax2.set_ylabel('Explained Variance Ratio')
    ax2.set_title('PCA Explained Variance')
    ax2.grid(True, alpha=0.3)

    plt.tight_layout()

    return fig, pca, X_scaled

# Hierarchical clustering analysis
def perform_hierarchical_clustering(X, species_names):
    """Perform hierarchical clustering analysis with optimal cluster determination"""

    # Standardize the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Find optimal number of clusters
    print("Finding optimal number of clusters...")

    # Create cluster optimization plots
    fig_opt, optimal_k_elbow, optimal_k_silhouette = plot_cluster_optimization(X_scaled)

    # Use the average of both methods, with preference for silhouette if they differ significantly
    if abs(optimal_k_elbow - optimal_k_silhouette) <= 1:
        optimal_k = int(np.mean([optimal_k_elbow, optimal_k_silhouette]))
    else:
        optimal_k = optimal_k_silhouette  # Prefer silhouette method

    print(f"Elbow method suggests: {optimal_k_elbow} clusters")
    print(f"Silhouette analysis suggests: {optimal_k_silhouette} clusters")
    print(f"Using: {optimal_k} clusters")

    # Perform hierarchical clustering
    linkage_matrix = linkage(X_scaled, method='ward')
    clusters = fcluster(linkage_matrix, optimal_k, criterion='maxclust')

    # Create cluster labels
    cluster_labels = [f'Metabolic Cluster {i}' for i in clusters]

    # Perform PCA for visualization
    pca = PCA()
    X_pca = pca.fit_transform(X_scaled)

    # Create visualization
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

    # PCA plot with hierarchical clusters
    colors = plt.cm.tab10(np.linspace(0, 1, optimal_k))
    for i in range(1, optimal_k + 1):
        mask = clusters == i
        ax1.scatter(X_pca[mask, 0], X_pca[mask, 1],
                   c=[colors[i-1]], label=f'Metabolic Cluster {i}',
                   alpha=0.7, s=60)

    # Highlight Gluconacetobacter brunescens Hr-1-5
    brunescens_idx = [i for i, name in enumerate(species_names)
                     if 'brunescens hr-1-5' in name.lower()]
    if brunescens_idx:
        idx = brunescens_idx[0]
        ax1.scatter(X_pca[idx, 0], X_pca[idx, 1],
                   c='red', s=200, marker='*',
                   edgecolors='black', linewidth=2)
        ax1.annotate('G. brunescens Hr-1-5', (X_pca[idx, 0], X_pca[idx, 1]),
                    xytext=(-70, 10), textcoords='offset points',
                    fontsize=10, ha='left', style='italic', weight='bold')

    ax1.set_xlabel(f'PC1 ({pca.explained_variance_ratio_[0]:.1%} variance)')
    ax1.set_ylabel(f'PC2 ({pca.explained_variance_ratio_[1]:.1%} variance)')
    ax1.set_title(f'Hierarchical Clustering Analysis (k={optimal_k})')
    ax1.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
    ax1.grid(True, alpha=0.3)

    # Dendrogram
    dendrogram(linkage_matrix, ax=ax2, truncate_mode='level', p=10)
    ax2.set_title('Hierarchical Clustering Dendrogram')
    ax2.set_xlabel('Sample Index')
    ax2.set_ylabel('Distance')
    ax2.axhline(y=linkage_matrix[-optimal_k+1, 2], color='red', linestyle='--',
                label=f'Cut for {optimal_k} clusters')
    ax2.legend()

    plt.tight_layout()

    return fig, fig_opt, clusters, cluster_labels, optimal_k

# K-means clustering analysis
def perform_kmeans_clustering(X, species_names):
    """Perform K-means clustering analysis with optimal cluster determination"""

    # Standardize the features
    scaler = StandardScaler()
    X_scaled = scaler.fit_transform(X)

    # Find optimal number of clusters
    fig_opt, optimal_k_elbow, optimal_k_silhouette = plot_cluster_optimization(X_scaled, " - K-means")

    # Use the average of both methods, with preference for silhouette if they differ significantly
    if abs(optimal_k_elbow - optimal_k_silhouette) <= 1:
        optimal_k = int(np.mean([optimal_k_elbow, optimal_k_silhouette]))
    else:
        optimal_k = optimal_k_silhouette  # Prefer silhouette method

    # Create detailed silhouette plot
    fig_sil, cluster_labels, silhouette_avg = plot_silhouette_detailed(X_scaled, optimal_k, " - K-means")

    return fig_opt, fig_sil, cluster_labels, optimal_k, silhouette_avg

# Main analysis function
def main_analysis():
    """Perform the complete analysis with optimal clustering"""
    print("Loading and preparing data...")
    data, X, feature_cols = load_and_prepare_data()

    # Get species names and assign groups
    species_names = data['Species'].tolist()
    groups, special_markers = assign_groups(species_names)

    print(f"Loaded data with {len(species_names)} species and {len(feature_cols)} features")
    print(f"Features included: {feature_cols}")
    print(f"Unique groups: {set(groups)}")

    # 1. Full dataset analysis
    print("\n1. Performing full dataset PCA analysis...")
    fig1, pca1, X_scaled1 = perform_pca_analysis(X, groups, special_markers, species_names)
    fig1.savefig('pca_full_dataset.png', dpi=300, bbox_inches='tight')
    fig1.savefig('pca_full_dataset.svg', bbox_inches='tight')
    plt.show()

    # 2. P1, P2, P7, P8 analysis
    print("\n2. Performing P1, P2, P7, P8 analysis...")
    p1278_features = [col for col in feature_cols if any(p in col for p in ['P1_', 'P2_', 'P7_', 'P8_'])]
    fig2, pca2, X_scaled2 = perform_pca_analysis(X, groups, special_markers, species_names,
                                                p1278_features, "P1, P2, P7, P8")
    fig2.savefig('pca_p1_p2_p7_p8.png', dpi=300, bbox_inches='tight')
    fig2.savefig('pca_p1_p2_p7_p8.svg', bbox_inches='tight')
    plt.show()

    # 3. P4, P5, P6, P7, P8 analysis
    print("\n3. Performing P4, P5, P6, P7, P8 analysis...")
    p45678_features = [col for col in feature_cols if any(p in col for p in ['P4_', 'P5_', 'P6_', 'P7_', 'P8_'])]
    fig3, pca3, X_scaled3 = perform_pca_analysis(X, groups, special_markers, species_names,
                                                p45678_features, "P4, P5, P6, P7, P8")
    fig3.savefig('pca_p4_p5_p6_p7_p8.png', dpi=300, bbox_inches='tight')
    fig3.savefig('pca_p4_p5_p6_p7_p8.svg', bbox_inches='tight')
    plt.show()

    # 4. P1, P2, P4, P5, P6 analysis
    print("\n4. Performing P1, P2, P4, P5, P6 analysis...")
    p12456_features = [col for col in feature_cols if any(p in col for p in ['P1_', 'P2_', 'P4_', 'P5_', 'P6_'])]
    fig4, pca4, X_scaled4 = perform_pca_analysis(X, groups, special_markers, species_names,
                                                p12456_features, "P1, P2, P4, P5, P6")
    fig4.savefig('pca_p1_p2_p4_p5_p6.png', dpi=300, bbox_inches='tight')
    fig4.savefig('pca_p1_p2_p4_p5_p6.svg', bbox_inches='tight')
    plt.show()

    # 5. Hierarchical clustering analysis with optimal cluster determination
    print("\n5. Performing hierarchical clustering analysis with optimal cluster determination...")
    fig5, fig5_opt, clusters_hier, cluster_labels_hier, optimal_k_hier = perform_hierarchical_clustering(X, species_names)
    fig5_opt.savefig('cluster_optimization_hierarchical.png', dpi=300, bbox_inches='tight')
    fig5_opt.savefig('cluster_optimization_hierarchical.svg', bbox_inches='tight')
    fig5.savefig('hierarchical_clustering_optimal.png', dpi=300, bbox_inches='tight')
    fig5.savefig('hierarchical_clustering_optimal.svg', bbox_inches='tight')
    plt.show()
    plt.show()

    # 6. K-means clustering analysis with optimal cluster determination
    print("\n6. Performing K-means clustering analysis with optimal cluster determination...")
    fig6_opt, fig6_sil, cluster_labels_kmeans, optimal_k_kmeans, silhouette_avg = perform_kmeans_clustering(X, species_names)
    fig6_opt.savefig('cluster_optimization_kmeans.png', dpi=300, bbox_inches='tight')
    fig6_opt.savefig('cluster_optimization_kmeans.svg', bbox_inches='tight')
    fig6_sil.savefig('silhouette_analysis_detailed.png', dpi=300, bbox_inches='tight')
    fig6_sil.savefig('silhouette_analysis_detailed.svg', bbox_inches='tight')
    plt.show()
    plt.show()

    # Create summary report
    print("\n" + "="*70)
    print("ANALYSIS SUMMARY WITH OPTIMAL CLUSTERING")
    print("="*70)

    # Count features by type
    cov_features = [col for col in feature_cols if 'Cov' in col]
    id_features = [col for col in feature_cols if 'ID' in col]
    eval_features = [col for col in feature_cols if 'E' in col]

    print(f"\nFeature breakdown:")
    print(f"  Coverage features: {len(cov_features)}")
    print(f"  Identity features: {len(id_features)}")
    print(f"  E-value features: {len(eval_features)} (log-transformed)")
    print(f"  Total features: {len(feature_cols)}")

    # Find Gluconacetobacter brunescens Hr-1-5
    brunescens_idx = [i for i, name in enumerate(species_names)
                     if 'brunescens hr-1-5' in name.lower()]

    if brunescens_idx:
        idx = brunescens_idx[0]
        print(f"\nGluconacetobacter brunescens Hr-1-5 found at index {idx}")
        print(f"Assigned to taxonomic group: {groups[idx]}")
        print(f"Hierarchical cluster: {cluster_labels_hier[idx]}")
        print(f"K-means cluster: Cluster {cluster_labels_kmeans[idx]}")

        # Show PC coordinates for each analysis
        print(f"\nPCA Coordinates:")
        print(f"Full dataset - PC1: {pca1.transform(X_scaled1)[idx, 0]:.3f}, PC2: {pca1.transform(X_scaled1)[idx, 1]:.3f}")

    print(f"\nClustering Results:")
    print(f"Optimal clusters (Hierarchical): {optimal_k_hier}")
    print(f"Optimal clusters (K-means): {optimal_k_kmeans}")
    print(f"Average silhouette score: {silhouette_avg:.3f}")

    print(f"\nTotal variance explained by first 2 PCs (full dataset): {pca1.explained_variance_ratio_[:2].sum():.1%}")

    # Group distribution
    print(f"\nTaxonomic group distribution:")
    for group in set(groups):
        count = groups.count(group)
        print(f"  {group}: {count} species")

    # Cluster distribution for hierarchical clustering
    print(f"\nHierarchical cluster distribution:")
    for i in range(1, optimal_k_hier + 1):
        count = list(clusters_hier).count(i)
        print(f"  Metabolic Cluster {i}: {count} species")

    # Show E-value transformation info
    print(f"\nE-value Processing:")
    print("  - E-values: 0 = best match (most significant), 1 = worst match")
    print("  - Zero E-values replaced with 1e-300 for log transformation")
    print("  - Applied -log10 transformation (better matches get higher values)")
    print("  - Missing E-values filled with worst significance (becomes 0 after transformation)")
    print("  - Infinite values capped at 300 for numerical stability")

    print("\n" + "="*70)
    print("Files saved:")
    print("- pca_full_dataset.png/svg")
    print("- pca_p1_p2_p7_p8.png/svg")
    print("- pca_p4_p5_p6_p7_p8.png/svg")
    print("- pca_p1_p2_p4_p5_p6.png/svg")
    print("- cluster_optimization_hierarchical.png/svg")
    print("- hierarchical_clustering_optimal.png/svg")
    print("- cluster_optimization_kmeans.png/svg")
    print("- silhouette_analysis_detailed.png/svg")
    print("="*70)

# Run the analysis
if __name__ == "__main__":
    main_analysis()
