import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import plotly.graph_objects as go
import matplotlib.patches as mpatches
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap
from matplotlib.lines import Line2D

def downsample_data(X, labels, subset_fraction=0.1):
    """
    Downsamples the data to a specified fraction.
    
    Parameters:
    - X (array-like): The feature matrix.
    - labels (dict): A dictionary of label arrays (e.g., {"acdc": acdc_labels, "manual": manual_labels}).
    - subset_fraction (float): Fraction of the data to keep (e.g., 0.1 for 10%).

    Returns:
    - X_subset (array): Subsampled feature matrix.
    - labels_subset (dict): Subsampled labels.
    """
    subset_size = int(X.shape[0] * subset_fraction)
    random_indices = np.random.choice(X.shape[0], subset_size, replace=False)
    X_subset = X[random_indices]
    labels_subset = {key: label[random_indices] for key, label in labels.items()}
    return X_subset, labels_subset

def compute_coordinates(X, n_pca_components=10):
    """
    Computes UMAP and t-SNE coordinates after PCA dimensionality reduction.

    Parameters:
    - X (array): Input data.
    - n_pca_components (int): Number of PCA components for dimensionality reduction.

    Returns:
    - umap_coords (array): UMAP coordinates.
    - tsne_coords (array): t-SNE coordinates.
    """
    pca = PCA(n_components=n_pca_components)
    X_pca = pca.fit_transform(X)
    umap_coords = umap.UMAP().fit_transform(X_pca)
    tsne_coords = TSNE(n_components=2, random_state=42).fit_transform(X_pca)
    return umap_coords, tsne_coords

def scatter_with_legend(ax, coords, labels, unique_labels, title, cmap="tab10"):
    """
    Creates a scatter plot with a legend.

    Parameters:
    - ax (Axes): Matplotlib Axes object.
    - coords (array): Coordinates for the scatter plot.
    - labels (array): Labels corresponding to the points.
    - unique_labels (array): Unique labels for creating the legend.
    - title (str): Title for the plot.
    - cmap (str): Colormap to use for the plot.
    """
    scatter = ax.scatter(
        coords[:, 0], coords[:, 1],
        c=[unique_labels.tolist().index(lbl) for lbl in labels],
        cmap=cmap, alpha=0.5, s=10
    )
    legend_labels = [Line2D([0], [0], marker='o', color='w', 
                            markerfacecolor=scatter.cmap(scatter.norm(i)), markersize=8)
                     for i, _ in enumerate(unique_labels)]
    ax.legend(
        legend_labels, unique_labels, title="Cell Types",
        loc="upper left", bbox_to_anchor=(1.05, 1)
    )
    ax.set_title(title)

def plot_dim_reduction(X, labels, results_dir=None, file_prefix=None):
    """
    Generates UMAP and t-SNE scatter plots with legends.

    Parameters:
    - X (array): Feature matrix.
    - labels (dict): Dictionary with label arrays (e.g., {"acdc": acdc_labels, "manual": manual_labels}).
    - results_dir (str, optional): Directory to save plots. If None, plots are displayed but not saved.
    - file_prefix (str, optional): Prefix for saved file names.
    """
    # Dimensionality reduction
    umap_coords, tsne_coords = compute_coordinates(X)

    # Plot UMAP
    fig, axes = plt.subplots(1, 2, figsize=(22, 8))
    scatter_with_legend(axes[0], umap_coords, labels['acdc'], np.unique(labels['acdc']), "UMAP of ACDC Marker Model")
    scatter_with_legend(axes[1], umap_coords, labels['manual'], np.unique(labels['manual']), "UMAP of Manual Annotation")
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if results_dir and file_prefix:
        plt.savefig(f"{results_dir}/{file_prefix}_umap.png")
    plt.show()
    plt.close()

    # Plot t-SNE
    fig, axes = plt.subplots(1, 2, figsize=(22, 8))
    scatter_with_legend(axes[0], tsne_coords, labels['acdc'], np.unique(labels['acdc']), "t-SNE of ACDC Marker Model")
    scatter_with_legend(axes[1], tsne_coords, labels['manual'], np.unique(labels['manual']), "t-SNE of Manual Annotation")
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    if results_dir and file_prefix:
        plt.savefig(f"{results_dir}/{file_prefix}_tsne.png")
    plt.show()
    plt.close()


# --- Compute Thresholds ---
def compute_thresholds(filtered_df):
    robust_scaled_df = filtered_df.drop('cell_type', axis=1).copy()
    thresholds = {}
    for marker in robust_scaled_df.columns:
        median = np.median(robust_scaled_df[marker])
        iqr = np.percentile(robust_scaled_df[marker], 75) - np.percentile(robust_scaled_df[marker], 25)
        robust_scaled_df[marker] = (robust_scaled_df[marker] - median) / iqr
        thresholds[marker] = np.percentile(robust_scaled_df[marker], 75)
    return robust_scaled_df, thresholds

# --- Plot Histograms ---
def plot_histograms(robust_scaled_df, thresholds, results_dir = "./"):
    fig, axes = plt.subplots(nrows=(len(robust_scaled_df.columns) // 3 + 1), ncols=3, figsize=(15, 4 * (len(robust_scaled_df.columns) // 3 + 1)))
    axes = axes.flatten()

    for idx, marker in enumerate(robust_scaled_df.columns):
        sns.histplot(robust_scaled_df[marker], bins=50, kde=True, color='skyblue', ax=axes[idx])
        axes[idx].axvline(thresholds[marker], color='red', linestyle='--', linewidth=2, label=f'Threshold: {thresholds[marker]:.2f}')
        axes[idx].set_title(f"{marker} Distribution")
        axes[idx].legend()

    for i in range(len(robust_scaled_df.columns), len(axes)):
        fig.delaxes(axes[i])

    plt.tight_layout()

    # Save or display the plot
    if results_dir:
        plt.savefig(os.path.join(results_dir, "marker_threshold_histograms.png"))
    plt.show()
    plt.close()
    

# --- Plot Heatmap ---
def plot_heatmap(X0, ct_index0, y0, marker_table, mode="marker_model", results_dir=None):
    """
    Plots a heatmap for either the marker model or manual annotations.
    
    Parameters:
    - X0 (ndarray): Data matrix.
    - ct_index0 (ndarray): Cell type index matrix from ACDC computation (used for marker model heatmap).
    - y0 (ndarray): Ground truth cell type annotations (used for manual annotations heatmap).
    - marker_table (DataFrame): Marker table with marker names as columns and cell types as rows.
    - mode (str): Either "marker_model" or "manual_annotations".
    - results_dir (str, optional): Directory to save the plot. If None, the plot is displayed but not saved.
    """
    if mode == "marker_model":
        avg_marker_model = np.vstack([np.mean(X0[ct_index0[:, i] == 1], axis=0) for i in range(len(marker_table.index))])
        feature_matrix = pd.DataFrame(avg_marker_model, columns=marker_table.columns, index=marker_table.index)
        title = "Average Heatmap from Marker Model"
        file_name = "marker_model_heatmap.png"
    elif mode == "manual_annotations":
        unique_cell_types = np.unique(y0)
        avg_manual = np.vstack([np.mean(X0[y0 == ct, :], axis=0) for ct in unique_cell_types])
        feature_matrix = pd.DataFrame(avg_manual, columns=marker_table.columns, index=unique_cell_types)
        title = "Average Heatmap from Manual Annotations"
        file_name = "manual_annotation_heatmap.png"
    else:
        raise ValueError("Invalid mode. Use 'marker_model' or 'manual_annotations'.")
    
    # Plotting the heatmap
    fig = plt.figure(figsize=(6, 4))
    sns.heatmap(feature_matrix, vmin=0.0, vmax=2.0, cmap="YlGnBu", center=0)
    plt.xticks(rotation=90)
    plt.title(title)
    
    # Save or display the plot
    if results_dir:
        plt.savefig(f"{results_dir}/{file_name}")
    plt.show()
    plt.close()
    
    
# --- Pie Chart ---
def plot_pie_chart(data, title, results_dir=None, file_name=None):
    """
    Plots a pie chart for the given data.
    
    Parameters:
    - data (Series or list): Data for plotting, typically a Pandas Series or list of cell types.
    - title (str): Title for the pie chart.
    - results_dir (str, optional): Directory to save the plot. If None, the plot is displayed but not saved.
    - file_name (str, optional): Name of the file to save the plot as. Required if results_dir is specified.
    """
    ax = pd.Series(data).value_counts().plot.pie(
        autopct="%1.1f%%",
        startangle=90,
        labels=pd.Series(data).value_counts().index,
        wedgeprops={'edgecolor': 'black'},
        figsize=(8, 8)
    )

    # Rotate the labels to avoid overlap
    for label in ax.get_xticklabels():
        label.set_rotation(45)

    plt.title(title)

    # Save or display the plot
    if results_dir and file_name:
        plt.savefig(f"{results_dir}/{file_name}")
    plt.show()
    plt.close()
    
    
    
def plot_cell_type_abundance_or_counts(predicted_cell_types, y0, plot_type='abundance', results_dir='./'):
    """
    Plot cell type distribution for ACDC marker model vs manual annotations.

    Parameters:
    - predicted_cell_types: List or array of predicted cell types (from marker model)
    - y0: List or array of manual annotations
    - plot_type: Either 'abundance' or 'counts'. 'abundance' shows the normalized proportions,
                 'counts' shows the actual counts.
    - results_dir: Directory where the plots will be saved.
    """
    
    # Calculate abundance (normalized counts) or cell counts based on the plot_type argument
    if plot_type == 'abundance':
        marker_abundance = pd.Series(predicted_cell_types).value_counts(normalize=True)
        manual_abundance = pd.Series(y0).value_counts(normalize=True)
    elif plot_type == 'counts':
        marker_abundance = pd.Series(predicted_cell_types).value_counts()
        manual_abundance = pd.Series(y0).value_counts()
    else:
        raise ValueError("plot_type must be either 'abundance' or 'counts'")

    # Prepare data for plotting
    abundance_df = pd.DataFrame({
        "ACDC Marker Model": marker_abundance,
        "Manual Annotation": manual_abundance
    }).fillna(0)

    # Create subplots: two plots side by side
    fig, axes = plt.subplots(1, 2, figsize=(15, 6))

    # Plot for ACDC Marker Model
    axes[0].bar(marker_abundance.index, marker_abundance.values, color='b', alpha=0.7)
    axes[0].set_title(f'ACDC Marker Model {plot_type.capitalize()}', fontsize=14)
    axes[0].set_xlabel('Cell Type', fontsize=12)
    axes[0].set_ylabel(f'{plot_type.capitalize()}', fontsize=12)
    axes[0].tick_params(axis='x', rotation=45)  # Rotate labels by 45 degrees

    # Plot for Manual Annotation
    axes[1].bar(manual_abundance.index, manual_abundance.values, color='g', alpha=0.7)
    axes[1].set_title(f'Manual Annotation {plot_type.capitalize()}', fontsize=14)
    axes[1].set_xlabel('Cell Type', fontsize=12)
    axes[1].set_ylabel(f'{plot_type.capitalize()}', fontsize=12)
    axes[1].tick_params(axis='x', rotation=45)  # Rotate labels by 45 degrees

    # Set tight layout to prevent overlap
    plt.tight_layout()
    
    # Save or display the plot
    if results_dir:
        plot_filename = f"side_by_side_{plot_type}_plot.png"
        plt.savefig(os.path.join(results_dir, plot_filename))
    plt.show()
    plt.close()
    

    
def create_sankey_data(data, source_col, target_col):
    """
    Prepares data for the Sankey diagram by grouping and mapping transitions.
    
    Parameters:
    - data (DataFrame): The input data containing source and target columns.
    - source_col (str): Name of the source column.
    - target_col (str): Name of the target column.
    
    Returns:
    - transitions (DataFrame): DataFrame with source, target, and count columns.
    - all_labels (list): List of unique labels used in the diagram.
    """
    # Group and count transitions
    transitions = data.groupby([source_col, target_col]).size().reset_index(name='count')

    # Get unique labels for source and target
    all_labels = pd.concat([transitions[source_col], transitions[target_col]]).unique()
    label_dict = {label: i for i, label in enumerate(all_labels)}

    # Map source and target labels to indices
    transitions['source'] = transitions[source_col].map(label_dict)
    transitions['target'] = transitions[target_col].map(label_dict)

    return transitions, all_labels

def plot_sankey(transitions, all_labels, title, results_dir='./'):
    """
    Plots a Sankey diagram using the given transitions data and labels.
    
    Parameters:
    - transitions (DataFrame): DataFrame containing source, target, and count columns.
    - all_labels (list): List of unique labels used in the diagram.
    - title (str): Title for the Sankey diagram.
    - output_path (str, optional): Path to save the plot as an HTML file. If None, the plot is displayed.
    """
    # Create the Sankey diagram
    fig = go.Figure(go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=list(all_labels)
        ),
        link=dict(
            source=transitions['source'],
            target=transitions['target'],
            value=transitions['count']
        )
    ))

    fig.update_layout(title_text=title, font_size=10)

    # Export to HTML or display the plot
    if results_dir:
        fig.write_html(results_dir)
        print(f"Plot saved as {results_dir}. Open this file to view the diagram.")
    else:
        fig.show()