library(data.table)
library(viridis)
library(circlize)

# --- Function: Transform cell types ---
transform_cell_types <- function(vec, dt, col) {
  transformed_vec <- dt$label[match(vec, dt[[col]])]
  transformed_vec <- factor(transformed_vec, levels = unique(dt$label))
  return(transformed_vec)
}

# --- Function: Create Chord Diagram ---
create_chord_diagram <- function(sce_cells, pred_cells, output_path) {
  message("Generating frequency matrix...")
  df <- data.frame("manual_gating" = sce_cells, "gating_tool" = pred_cells)
  freq_matrix <- as.data.frame.matrix(table(df$manual_gating, df$gating_tool))
  
  message("Creating Chord Diagram: ", output_path)
  png(output_path, width = 3000, height = 3000, res = 300)  # Increase dimensions for clarity
  circos.clear()
  circos.par(gap.degree = 5, track.margin = c(0.01, 0.01), points.overflow.warning = FALSE)
  
  chordDiagram(
    freq_matrix,
    grid.col = viridis::viridis(nrow(freq_matrix) + ncol(freq_matrix)),
    transparency = 0.4, 
    annotationTrack = "grid", 
    preAllocateTracks = list(track.height = 0.05),
    cex = 0.6  # Reduce font size for entity labels
  )
  title("Chord Diagram: Manual Gating vs Gating Tool")
  dev.off()
  message("Chord Diagram saved: ", output_path)
}

# --- Main Function ---
main <- function() {
  message("Reading mapping file...")
  mapper_path <- "~/mapping_cell_types_labels.csv"
  cell_type_dt <- fread(mapper_path, na.strings = c(NA_character_, ""))
  
  message("Loading SCE data...")
  sce <- readRDS("/nfs/proj/collab_anke_hannover/Manual_Gating_Noemi/sce.rds")
  sce_cells <- transform_cell_types(sce$cell_type, cell_type_dt, "sce_label")
  
  message("Loading results from CSV files...")
  folder_path <- "~/results_celltypes"
  file_pattern <- "*.csv"
  filenames <- list.files(folder_path, pattern = file_pattern, full.names = TRUE)
  
  results_list <- lapply(filenames, function(file) {
    lines <- readLines(file)
    lines <- lines[lines != ""]
    lines <- gsub('"', '', lines)
    return(lines)
  })
  names(results_list) <- unlist(lapply(filenames, function(file) sub(file_pattern, '', basename(file))))
  
  pred_cells_list <- lapply(results_list, transform_cell_types, cell_type_dt, "matrix_label")
  
  message("Creating output directory for chord plots...")
  dir.create("~/chord_plots", showWarnings = FALSE)

  lapply(names(pred_cells_list), function(name) {
    message("Processing comparison: ", name)
    pred_cells <- pred_cells_list[[name]]
    output_path_chord <- paste0("~/chord_plots/chord_", name, ".png")
    message("Generating diagram for: ", name)
    create_chord_diagram(sce_cells, pred_cells, output_path_chord)
  })
  message("All chord diagrams generated!")
}

# --- Run the script ---
main()
