library(data.table)
library(viridis)
library(circlize)
library(SingleCellExperiment)
library(RColorBrewer)


# --- Function: Transform cell types ---
transform_cell_types <- function(vec, dt, col) {
  transformed_vec <- dt$label[match(vec, dt[[col]])]
  transformed_vec <- factor(transformed_vec, levels = unique(dt$label))
  return(transformed_vec)
}

# --- Function: Create Chord Diagram ---
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

  # Hier wird die Legende nach dem Schließen des Diagramms erstellt
  legend_output_path <- paste0(output_path,"_legend.png")
  
  message("Creating legend: ", legend_output_path)
  # Erstelle ein separates Bild für die Legende
  png(legend_output_path, width = 3000, height = 3000, res = 300)  
  # Generiere eine neue Plotfläche
  plot.new()
  # Legende erstellen
  legend("topright", 
         legend = unique_labels, 
         fill = label_colors, 
         title = "Cell Type", 
         cex = 1.2, 
         border = "black", 
         box.lwd = 2)  # Anpassungen der Legende
  dev.off()
  
  message("Legend saved: ", legend_output_path)
}



# --- Main Function ---
main <- function() {
  message("Reading mapping file...")
  mapper_path <- "~/mapping_cell_types_labels.csv"
  cell_type_dt <- fread(mapper_path, na.strings = c(NA_character_, ""))
  
  message("Loading SCE data...")
  sce <- readRDS("/nfs/proj/collab_anke_hannover/Manual_Gating_Noemi/sce.rds")

  sce_cells <- transform_cell_types(sce$cell_type, cell_type_dt, "sce_label")

  output_file <- "manual_sce_unscaled_original.csv"

  # Schreibe sce$cell_type in die Datei
  writeLines(as.character(sce_cells), con = output_file)

  message("Cell types erfolgreich in ", output_file, " geschrieben.")
  
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
    output_path_chord <- paste0("~/chord_plots/", name, "_dark.png")
    message("Generating diagram for: ", name)
    create_chord_diagram(sce_cells, pred_cells, output_path_chord, cell_type_dt)
  })
  message("All chord diagrams generated!")
}

# --- Run the script ---
main()
