import os
import pandas as pd
import scanpy as sc
import numpy as np
import time
import anndata as ad
import sys
sys.path.append('/dss/dsshome1/0F/di93quv/Systems_biomedicine/acdc/')
from ACDC.cell_type_annotation import *
from ACDC.random_walk_classifier import *
from ACDC.visualization import *
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import umap.umap_ as umap
import phenograph

# --- Setup directories ---
results_dir = "./results/"
os.makedirs(results_dir, exist_ok=True)
log_file = os.path.join(results_dir, "output_log.txt")

def log(message):
    with open(log_file, "a") as f:
        f.write(f"{message}\n")

# --- Load marker table ---
log("Loading marker table...")
marker_table = pd.read_csv("/dss/dsshome1/0F/di93quv/Systems_biomedicine/acdc/data/Manual_Gating_Noemi/cell_types.csv")
marker_table.set_index('Cell type', inplace=True)

# --- Load and preprocess data ---
log("Loading data...")
adata = sc.read_h5ad("/dss/dsshome1/0F/di93quv/Systems_biomedicine/acdc/data/Manual_Gating_Noemi/sce_scaled.h5ad")
#adata.X = adata.layers['scaled_exprs'] #scaled
adata.X = adata.layers['exprs'] # uncsclaed
expr_df = pd.DataFrame(adata.X, index=adata.obs_names, columns=adata.var_names)
filtered_df = pd.concat([adata.obs, expr_df], axis=1)
filtered_df = filtered_df.drop(columns=['sample_id', 'group', 'patient_id', 'timepoint', 'measurement_day', 'condition', 'cluster_id'])
columns = [col for col in filtered_df.columns if col != 'cell_type'] + ['cell_type']
filtered_df = filtered_df[columns]
filtered_df.to_csv(os.path.join(results_dir, "filtered_expression_with_metadata.csv"), index=False)

# --- Compute thresholds ---
robust_scaled_df, thresholds = compute_thresholds(filtered_df)

# Plot histograms
log("Plotting histograms...")
plot_histograms(robust_scaled_df, thresholds)

# --- Marker model computation ---
log("Starting ACDC...")
start_time = time.time()
X0 = robust_scaled_df.to_numpy()
y0 = filtered_df['cell_type'].values
mk_model = compute_marker_model(robust_scaled_df, marker_table, 0.0)
score0 = get_score_mat(X0, [], marker_table, [], mk_model)
score0 = np.concatenate([score0, 1.0 - score0.max(axis=1)[:, np.newaxis]], axis=1)
thres = 0.4
ct_index0 = get_unique_index(X0, score0, marker_table, thres)
res_c = get_landmarks(X0, score0, ct_index0, marker_table.index.tolist(), phenograph, thres)
landmark_mat, landmark_label = output_feature_matrix(res_c, marker_table.index.tolist())
end_time = time.time()
log(f"ACDC runtime: {end_time - start_time:.2f} seconds")

# --- Create a new AnnData object ---
log(" Create a new AnnData object...")
filtered_data = filtered_df.drop(columns=['cell_type']).to_numpy()
filtered_obs = filtered_df[['cell_type']].copy()
filtered_var = adata.var
filtered_adata = ad.AnnData(X=filtered_data, obs=filtered_obs, var=filtered_var)

# --- Assign ACDC predictions ---
log(" Assign ACDC predictions...")
predicted_cell_types = [marker_table.index[i] if i < len(marker_table.index) else 'unknown' for i in ct_index0.argmax(axis=1)]
filtered_adata.obs['ACDC_cell_type'] = predicted_cell_types
filtered_adata.obs.to_csv(
    os.path.join(results_dir, "filtered_adata_obs.csv"), 
    index=False  
)

# --- Heatmaps ---
# Plot and save the heatmap for the marker model
log(" Plotting heatmaps...")
plot_heatmap(X0, ct_index0, y0, marker_table, mode="marker_model", results_dir=results_dir)

# Plot and save the heatmap for manual annotations
plot_heatmap(X0, ct_index0, y0, marker_table, mode="manual_annotations", results_dir=results_dir)


# --- Downsample data for visualization ---
# For UMAP/t-SNE, subset to around 100k cells (1% cells)
log(" Plotting UMAP and TSNE...")
subset_fraction = min(100000 / 16800000, 1.0)
X_filtered = filtered_adata.X
X_subset, labels_subset = downsample_data(
    X_filtered, 
    labels={"acdc": filtered_adata.obs["ACDC_cell_type"], "manual": filtered_adata.obs["cell_type"]}, 
    subset_fraction=subset_fraction
)

# Generate UMAP and t-SNE plots
plot_dim_reduction(X_subset, labels_subset, results_dir=results_dir, file_prefix="dim_reduction_filtered_unscaled")

# Generate pie charts and cell type abundance plots
log(" Plotting pie chart and abundance plots...")
plot_pie_chart(data=predicted_cell_types, title="ACDC Marker Model Cell Type Distribution", results_dir=results_dir, file_name="pie_acdc_model_unscaled.png")
plot_pie_chart(data=y0, title="Manual Cell Type Distribution", results_dir=results_dir, file_name="pie_manual_unscaled.png")
plot_cell_type_abundance_or_counts(predicted_cell_types, y0, plot_type='abundance', results_dir=results_dir)

# Generate Sankey plot
log(" Plotting sankey plot...")
transitions, all_labels = create_sankey_data(filtered_df, 'cell_type', 'ACDC_cell_type')
plot_sankey(transitions, all_labels, "Sankey Diagram: Cell Type vs ACDC Cell Type", os.path.join(results_dir, "sankey_diagram_unscaled.html"))