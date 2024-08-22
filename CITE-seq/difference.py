import pandas as pd
import matplotlib.pyplot as plt
import scanpy as sc
import os

def distribution_stackedbar(rna, sample_col='sample', celltype_col='celltype_V1',figsize=(8, 6)):
    """
    Plots a stacked barplot showing the percentage of each cell type in each sample.
    
    Parameters:
    - rna: AnnData object containing the observations (`obs`) and uns (`uns`) data.
    - sample_col: The column in `rna.obs` that contains sample identifiers.
    - celltype_col: The column in `rna.obs` that contains cell type annotations.
    - figsize: A tuple specifying the figure size.
    """
    
 
    df = pd.DataFrame(rna.obs)

    # Calculate the percentage 
    count_df = df.groupby([sample_col, celltype_col]).size().unstack(fill_value=0)
    percent_df = count_df.div(count_df.sum(axis=1), axis=0) * 100
    file_path = f"/data/msun/projects/Stephen/PDAC_scRNA/02_difference/cell_distribution/celltype_by{sample_col}.csv"
    if os.path.exists(file_path):
        file_path = file_path.replace(".csv", "_1.csv")
        for i in range(2,10):
            if not os.path.exists(file_path):
                break
            file_path = file_path.replace(f"_{i-1}.csv", f"_{i}.csv")
    count_df.to_csv(file_path, index=False)
    print(f"DataFrame saved to {file_path}")
    
    # Extract the colors for cell types
    color_key = celltype_col + "_colors"
    cell_type_colors = rna.uns[color_key]
    cell_types = rna.obs[celltype_col].cat.categories
    color_map = dict(zip(cell_types, cell_type_colors))

    # Calculate total cell counts per sample
    total_cells_per_sample = count_df.sum(axis=1)
    x_labels = [f'{sample}\n(n={total_cells})' for sample, total_cells in total_cells_per_sample.items()]

    fig, ax = plt.subplots(figsize=figsize)
    percent_df.plot(kind='bar', stacked=True, color=[color_map[ct] for ct in percent_df.columns], ax=ax)
    ax.set_xticklabels(x_labels, rotation=0)

    # Customize the plot
    plt.title('Percentage of Cell Types in Each Sample')
    plt.xlabel('Sample')
    plt.ylabel('Percentage')
    plt.legend(title='Cell Type', bbox_to_anchor=(1.05, 1), loc='upper left')
    plt.tight_layout()

    plt.show()
    plt.close()
    

def distribution_barplot(adata, cell_type_col, condition_col,figsize=(14, 7), bar_width=0.3):
    """
    Plots the percentage of each cell type under different conditions with custom colors.
    
    Parameters:
    adata (AnnData): The AnnData object containing the single-cell data.
    cell_type_col (str): The column name in adata.obs that contains cell type information.
    condition_col (str): The column name in adata.obs that contains condition information.
    condition_colors (dict): A dictionary mapping conditions to their respective colors.
    
    Returns:
    None
    """
    # Calculate total cells per condition
    total_cells_per_condition = adata.obs.groupby(condition_col).size()
    
    # Calculate cell counts per cell type and condition
    cell_counts = adata.obs.groupby([cell_type_col, condition_col]).size().unstack(fill_value=0)
    
    # Calculate percentages
    cell_type_percentages = cell_counts.div(total_cells_per_condition, axis=1) * 100
    
    # Plotting
    fig, ax = plt.subplots(figsize= figsize) 
    colors = ["#E8AECB","#263260","#D6C685","#FFFFFF"]
    index = np.arange(len(cell_type_percentages))

    for i, condition in enumerate(cell_type_percentages.columns):
        ax.bar(index + i * bar_width, cell_type_percentages[condition], 
               bar_width, label=condition, color=colors[i])

    ax.set_xlabel('Cell Type')
    ax.set_ylabel('Percentage of cell number by total cell number under one condition')
    ax.set_title('Cell Type Distribution Across Conditions')
    ax.set_xticks(index + bar_width)
    ax.set_xticklabels(cell_type_percentages.index, rotation=45, ha='right')
    ax.legend(title='Condition')

    plt.tight_layout()
    plt.show()

def split_plots(rna, celltype_col='celltype_V1', sample_col='sample', n_cols=3, figsize=(15, 12), ptsize=26):
    """
    Plots UMAPs for each condition in a grid layout with a single shared legend.
    
    Parameters:
    - rna: AnnData object containing the data.
    - celltype_col: The column in `rna.obs` used for coloring (default is 'celltype_V1').
    - sample_col: The column in `rna.obs` that contains sample identifiers (default is 'sample').
    - n_cols: Number of columns for the grid layout (default is 3).
    - figsize: Size of the entire figure (default is (15, 12)).
    - ptsize: Size of the points in the UMAP plot.
    """
    # Get the unique conditions
    conditions = rna.obs[sample_col].unique()
    
    # Determine the number of rows needed based on the number of conditions
    n_rows = -(-len(conditions) // n_cols)  # This is a ceiling division
    
    # Create a figure with subplots
    fig, axes = plt.subplots(n_rows, n_cols, figsize=figsize)
    axes = axes.flatten()  # Flatten the axes array for easy iteration

    # Plot UMAP for each condition without individual legends
    for i, cond in enumerate(conditions):
        ax = axes[i]
        sc.pl.umap(rna[rna.obs[sample_col] == cond], 
                   title=f'UMAP - {cond}', 
                   color=[celltype_col], 
                   size=ptsize, 
                   ax=ax, 
                   show=False)  

    # Hide any unused axes
    for j in range(i + 1, len(axes)):
        axes[j].axis('off')
    
    plt.tight_layout()  
    plt.show()
    plt.close()


import scanpy as sc
import pandas as pd
import numpy as np
from itertools import combinations

def findDEG(adata, cluster_key, condition_key, method='wilcoxon',layers = "rna"):
    """
    Perform differential expression analysis between conditions for each cluster.

    Parameters:
    adata (AnnData): The annotated data matrix.
    cluster_key (str): The key in adata.obs that indicates the cluster assignments.
    condition_key (str): The key in adata.obs that indicates the conditions (e.g., treatment, control).
    method (str): The method to use for differential expression testing (e.g., 't-test', 'wilcoxon').

    Returns:
    dict: A dictionary where keys are cluster names and values are DataFrames with differential expression results.
    """

    # Get the unique conditions in the dataset
    conditions = sorted(adata.obs[condition_key].unique())
    
    # Iterate over each pair of conditions
    for cond1, cond2 in combinations(conditions, 2):
        print(f"Processing comparison: {cond1} vs {cond2}")

        # Create output directories if they don't exist
        output_dir = f"/data/msun/projects/Stephen/PDAC_scRNA/02_difference/DEG/{layers}_{cond1}_vs_{cond2}"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
        output_dir0 = f"/data/msun/projects/Stephen/PDAC_scRNA/02_difference/volcano/{layers}_{cond1}_vs_{cond2}"
        if not os.path.exists(output_dir0):
            os.makedirs(output_dir0)

        # Iterate over each cluster
        for cluster in adata.obs[cluster_key].unique():
            print(f"  Processing cluster: {cluster}")

            try:
                # Subset the data to the current cluster
                adata_cluster = adata[adata.obs[cluster_key] == cluster]
                adata_pair = adata_cluster[adata_cluster.obs[condition_key].isin([cond1, cond2])]
                # Only set adata_pair.X if layers is not "prot"
                if layers != "prot":
                    adata_pair.X = adata_pair.layers['log1p_norm'].toarray()
                
                # Perform differential expression testing
                sc.tl.rank_genes_groups(adata_pair, groupby=condition_key, groups=[cond1], reference=cond2, method=method)
                
                # Store the results in a DataFrame
                de_results = sc.get.rank_genes_groups_df(adata_pair, group=cond1)
                de_results['cluster'] = cluster
                de_results['condition_1'] = cond1
                de_results['condition_2'] = cond2
    
                # Define the filename
                file_name = f"{output_dir}/{cluster}.csv"
                
                # Write the results to a CSV file
                de_results.to_csv(file_name, index=False)
                print(f"Results saved to {file_name}")
    
                # Create and save a volcano plot
                plot_volcano(de_results, cluster, cond1, cond2)

            except ValueError as e:
                print(f"\033[91mSkipping cluster {cluster} for comparison {cond1} vs {cond2} due to error: {e}\033[0m")


def plot_volcano(de_results, cluster, cond1, cond2, layers = "rna"):
    """
    Create and save a volcano plot for the differential expression results.

    Parameters:
    de_results (DataFrame): The differential expression results DataFrame.
    cluster (str): The cluster identifier.
    cond1 (str): The first condition.
    cond2 (str): The second condition.
    output_dir (str): The directory where the plots will be saved.
    """

    plt.figure(figsize=(10, 6))

    # Use log2 fold change and p-value to create a volcano plot
    log_fold_change = de_results['logfoldchanges']
    p_values = de_results['pvals_adj']

    # Plot the data
    plt.scatter(log_fold_change, -np.log10(p_values), c='gray', alpha=0.7)

    # Highlight significant genes
    significant_pos = (de_results['pvals_adj'] < 0.05) & (de_results['logfoldchanges'] > 1)
    significant_neg = (de_results['pvals_adj'] < 0.05) & (de_results['logfoldchanges'] < -1)
    plt.scatter(log_fold_change[significant_pos], -np.log10(p_values[significant_pos]), c='red')
    plt.scatter(log_fold_change[significant_neg], -np.log10(p_values[significant_neg]), c='blue')

    # Label top 10 most significant genes
    de_results['rank_metric'] = -np.log10(de_results['pvals_adj']) * np.abs(de_results['logfoldchanges'])
    top_pos = de_results.loc[significant_pos].nlargest(10, 'rank_metric')
    top_neg = de_results.loc[significant_neg].nlargest(10, 'rank_metric')
    top_genes = pd.concat([top_pos, top_neg])
    for i, row in top_genes.iterrows():
        plt.text(row['logfoldchanges'], -np.log10(row['pvals_adj']), row['names'], fontsize=8, ha='right')

    
    # Add labels and title
    plt.title(f"Volcano Plot: {cluster}, {cond1} vs {cond2}")
    plt.xlabel('Log2 Fold Change')
    plt.ylabel('-Log10(p-value)')

    # Save the plot
    output_dir0 = f"/data/msun/projects/Stephen/PDAC_scRNA/02_difference/volcano/{layers}_{cond1}_vs_{cond2}"
    plot_file_name = f"{output_dir0}/{cluster}.pdf"
    plt.savefig(plot_file_name)
    plt.close()

    print(f"Volcano plot saved to {plot_file_name}")




def compare_deg_counts(adata, cluster_key, condition_key, cond1, cond2, pval_threshold=0.05, logfc_threshold=1.0,layers = "rna"):
    """
    Compare the number of DEGs across different cell types for a single pair of conditions and visualize it with a bar plot.

    Parameters:
    adata (AnnData): The annotated data matrix.
    cluster_key (str): The key in adata.obs that indicates the cluster assignments.
    condition_key (str): The key in adata.obs that indicates the conditions (e.g., treatment, control).
    cond1 (str): The first condition in the comparison.
    cond2 (str): The second condition in the comparison.
    pval_threshold (float): The p-value threshold for considering a gene as differentially expressed.
    logfc_threshold (float): The log fold change threshold for considering a gene as differentially expressed.

    Returns:
    DataFrame: A DataFrame with the number of DEGs per cell type.
    """

    deg_counts = []

    # Iterate over each cluster
    for cluster in adata.obs[cluster_key].unique():
        print(f"Processing cluster: {cluster}")

        try: 
            # Subset the data to the current cluster
            adata_cluster = adata[adata.obs[cluster_key] == cluster]
            # Only set adata_pair.X if layers is not "prot"
            if layers != "prot":
                adata_cluster.X = adata_cluster.layers['log1p_norm'].toarray()
            
            # Perform differential expression analysis between the two conditions
            sc.tl.rank_genes_groups(adata_cluster, groupby=condition_key, groups=[cond1], reference=cond2, method='wilcoxon')
            
            # Extract the results for the comparison
            de_results = sc.get.rank_genes_groups_df(adata_cluster, group=cond1)
            
            # Filter based on p-value and log fold change thresholds
            significant_degs = de_results[(de_results['pvals_adj'] < pval_threshold) & (abs(de_results['logfoldchanges']) > logfc_threshold)]
            num_degs = significant_degs.shape[0]
            
            # Append the number of DEGs for this cluster
            deg_counts.append((cluster, num_degs))
        except ValueError as e:
            print(f"\033[91mSkipping cluster {cluster} for comparison {cond1} vs {cond2} due to error: {e}\033[0m")
    
    # Create a DataFrame from the results
    deg_counts_df = pd.DataFrame(deg_counts, columns=['Cell Type', 'DEG Count'])

    # Plot the results as a bar plot
    plt.figure(figsize=(5, 4))
    bars = plt.barh(deg_counts_df['Cell Type'], deg_counts_df['DEG Count'], color='skyblue')
    plt.xlabel('Number of DEGs')
    plt.ylabel('Cell Type')
    plt.title(f'Number of DEGs in {cond1} vs {cond2} across Cell Types')
    plt.gca().invert_yaxis()  # Invert y-axis to have the largest bar at the top
    
    # Add numbers on top of the bars
    for bar in bars:
        plt.text(bar.get_width() + 0.1, bar.get_y() + bar.get_height()/2,
                 f'{int(bar.get_width())}', va='center', fontsize=10, color='black')

    plt.show()

    return deg_counts_df
