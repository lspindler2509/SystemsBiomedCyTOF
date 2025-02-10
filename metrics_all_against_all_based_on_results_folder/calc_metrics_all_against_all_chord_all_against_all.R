library(tibble)
library(data.table)
library(fossil) # ARI
library(yardstick) # balanced accuracy
library(MLmetrics) # F1 score
library(mclust)
library(irr)

library(viridis)
library(circlize)
library(SingleCellExperiment)
library(RColorBrewer)

# Input parameters
mapper_path <- "~/mapping_cell_types_labels.csv"     # File with uniform cell type labels as mapper
folder_path <- "~/results_celltypes_all"            # Folder with the result file(s) for the different tools/methods/etc.
file_pattern <- "*.csv"                             # Pattern of the wanted file(s) (! not the complete filename)
out_path <- "metrics_all_against_all_original_with_weighted__intersection.csv"           # File to save the results of the metrics

# Functions

# Mapping the cell labels to the uniform labels
trans_cell_types <- function(vec, dt, col) {
  trans_vec <- dt$label[match(vec, dt[[col]])]
  trans_vec <- factor(trans_vec, levels = unique(dt$label)[unique(dt$label) != ""]) # Remove empty levels
  return(trans_vec)
}

create_chord_diagram <- function(sce_cells, pred_cells, output_path, cell_type_dt) {
  message("Generating frequency matrix...")
  df <- data.frame("manual_gating" = sce_cells, "gating_tool" = pred_cells)
  freq_matrix <- as.matrix(table(df$manual_gating, df$gating_tool))

  # Farben für Labels definieren
  unique_labels <- unique(as.vector(cell_type_dt$label))
  label_colors <- c("#ff6db6", "#004949", "#db6d00",  "#B2DF8A", "#FDB462", "#490092", "#009999", "#8f4e00", "#ffdf4d","#b66dff")
  names(label_colors) <- unique_labels
  
  # Plot starten
  message("Creating Chord Diagram: ", output_path)
  png(output_path, width = 3000, height = 3000, res = 300)  # Erhöhte Auflösung
  # Chord Diagram
  chordDiagram(
    freq_matrix,
    directional = 1,
    diffHeight = mm_h(8),
    self.link = 1,
    grid.col = label_colors,
    target.prop.height = mm_h(4)
  )
  
  # Titel mit Dateinamen
  formatted_title <- gsub("_", " ", basename(output_path))
  formatted_title <- gsub(".png", "", formatted_title)
  title(paste("Chord Diagram:", formatted_title))
  message("Chord Diagram saved: ", output_path)
  dev.off()
}

# Calculating the metrics for pairwise comparisons
calc_pairwise_metrics <- function(method1_vec, method2_vec, method1, method2) {
  # Adjusted Rand Index
  if (all.equal(method1_vec, method2_vec) == TRUE) {
    return(data.table(
      "method1" = method1,
      "method2" = method2,
      "ACC" = 1,
      "ARI" = 1,
      "F1" = 1,
      "kappa" = 1,
      "F1_weighted" = 1,
      "ACC_weighted" = 1
    ))
  }

  # Finden der Schnittmenge der tatsächlichen Werte beider Vektoren
  common_values <- intersect(unique(method1_vec), unique(method2_vec))
  print(common_values)
  
  # Filtere die Vektoren, so dass nur die gemeinsamen Werte verbleiben
  method1_vec_filtered <-  factor(method1_vec, levels = common_values)
  method2_vec_filtered <-  factor(method2_vec, levels = common_values)

  kappa_result <- kappa2(cbind(method1_vec, method2_vec))$value
  # Balanced Accuracy
  acc <- bal_accuracy_vec(method1_vec_filtered, method2_vec_filtered, na.rm = TRUE)
  acc_weighted <- bal_accuracy_vec(method1_vec_filtered, method2_vec_filtered, na.rm = TRUE, estimator = "macro_weighted")

  ari_cluster <- adjustedRandIndex(method1_vec, method2_vec)

  f1_yardstick <- f_meas(data.frame(truth = method1_vec_filtered, estimate = method2_vec_filtered), truth = "truth", estimate = "estimate", na.rm = TRUE)
  f1_yardstick_weighted <- f_meas(data.frame(truth = method1_vec_filtered, estimate = method2_vec_filtered), truth = "truth", estimate = "estimate", na.rm = TRUE, estimator = "macro_weighted")

  
  return(data.table("method1" = method1, "method2" = method2, "ACC" = acc, "ARI" = ari_cluster, "F1" = f1_yardstick$.estimate, "kappa" = kappa_result, "F1_weighted" = f1_yardstick_weighted$.estimate, "ACC_weighted" = acc_weighted))
}

# Extracting the result vectors from the file(s)
filenames <- list.files(folder_path, pattern = file_pattern, full.names = TRUE)
results_list <- lapply(filenames, function(file) {
  lines <- readLines(file)
  
  # Remove empty lines, if any
  lines <- lines[lines != ""]
  
  # Remove quotes from lines
  lines <- gsub('"', '', lines)
  return(lines)
})
names(results_list) <- unlist(lapply(filenames, function(file) sub(file_pattern, '', basename(file))))

# Translating the cell labels to uniform names
cell_type_dt <- fread(mapper_path, na.strings = c(NA_character_, "")) # Reading in the mapper data table
pred_cells <- lapply(results_list, trans_cell_types, cell_type_dt, "matrix_label")

# Performing all-vs-all comparisons using expand.grid
methods <- names(pred_cells)
print(methods)

# Create pairwise comparisons and calculate metrics (both directions)
all_combinations <- expand.grid(method1 = methods, method2 = methods)
print(all_combinations)

# Calculate metrics for both directions
all_metrics <- apply(all_combinations, 1, function(pair) {
  method1 <- pair[1]
  method2 <- pair[2]
  
  cat("Calculating metrics for comparison:", method1, "vs", method2, "\n")
  
  # Two methods to compare
  vec1 <- pred_cells[[method1]]
  vec2 <- pred_cells[[method2]]
  
  # Compute metrics for both directions
  metrics1 <- calc_pairwise_metrics(vec1, vec2, method1, method2)
  create_chord_diagram(vec1, vec2, paste0("~/all_against_all_chord/", method1,"_", method2, ".png"), cell_type_dt)

  print(metrics1)
  
  return(list(metrics1))
})

print(all_metrics)

# Flatten the list and remove NULL values
all_metrics_flat <- unlist(all_metrics, recursive = FALSE)
all_metrics_flat <- all_metrics_flat[!sapply(all_metrics_flat, is.null)]

# Combine all results into a single data table
comp_metrics <- rbindlist(all_metrics_flat)

# Save the metrics results in a file
fwrite(comp_metrics, file = out_path)

cat("Metrics calculation complete. Results saved to:", out_path, "\n")
