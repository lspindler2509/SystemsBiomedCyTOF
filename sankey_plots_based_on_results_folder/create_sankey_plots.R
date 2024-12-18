library(data.table)
library(ggplot2)
library(ggsankey)
library(dplyr)
library(tidyr)
library(viridis)
library(stringr)

trans_cell_types <- function(vec, dt, col) {
  trans_vec <- dt$label[ match(vec, dt[[col]]) ]
  trans_vec <- factor(trans_vec, levels = unique(dt$label))
  return(trans_vec)
}

prepare_sankey_data <- function(sce_cells, pred_cells) {
  # Erstelle Dataframe mit manual_gating und prediction
  df <- data.frame("manual_gating" = sce_cells, "gating_tool" = pred_cells)
  
  # Konvertiere in long format für Sankey-Plot
  df_long <- df %>% make_long(manual_gating, gating_tool)
  
  # Berechne Prozentsätze je Node
  df_long_pct <- df_long %>% 
    group_by(node) %>%
    tally() %>%
    mutate(pct = n / sum(n))
  
  # Merge Prozentsätze zurück in den Long DataFrame
  df_plot <- merge(df_long, df_long_pct, by = 'node', all.x = TRUE)
  
  return(df_plot)
}

# --- Funktion 3: Sankey-Plot erstellen ---
create_sankey_plot <- function(df_plot, output_path) {
  p <- ggplot(df_plot, aes(x = x,
                           next_x = next_x,
                           node = node,
                           next_node = next_node,
                           fill = factor(node),
                           label = paste0(node, " n=", n, '(', round(pct * 100, 1), '%)' ))) +
    geom_sankey(flow.alpha = 0.5, color = "gray40", show.legend = TRUE) +
    geom_sankey_label(size = 2, color = "black", fill = "white", hjust = -0.1) +
    theme_bw() +
    theme(legend.position = "none",
          axis.title = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank(),
          plot.background = element_rect(fill = "white", color = NA)) + # Weißer Hintergrund
    labs(fill = 'Nodes') +
    scale_fill_viridis_d(option = "inferno")
  
  # Speichere Plot als PNG
  ggsave(output_path, p, width = 8, height = 6, dpi = 300, bg = "white")
}

main <- function() {
  # Cell type Mapper einlesen
  mapper_path <- "~/mapping_cell_types_labels.csv"
  cell_type_dt <- fread(mapper_path, na.strings = c(NA_character_, ""))
  
  # Manual Gating Daten aus SCE
  sce <- readRDS("/nfs/proj/collab_anke_hannover/Manual_Gating_Noemi/sce.rds")
  sce_cells <- trans_cell_types(sce$cell_type, cell_type_dt, "sce_label")
  
  # Ergebnisse aus Dateien einlesen
  folder_path <- "~/results_celltypes"
  file_pattern <- "*.csv"
  filenames <- list.files(folder_path, pattern = file_pattern, full.names = TRUE)
  
  # Ergebnisse transformieren und für jeden Vergleich Sankey-Plot erstellen
  results_list <- lapply(filenames, function(file) {
    lines <- readLines(file)
    lines <- lines[lines != ""]
    lines <- gsub('"', '', lines) # Entferne Anführungszeichen
    return(lines)
  })
  names(results_list) <- unlist(lapply(filenames, function(file) sub(file_pattern, '', basename(file))))
  
  pred_cells_list <- lapply(results_list, trans_cell_types, cell_type_dt, "matrix_label")
  
  dir.create("~/sankey_plots", showWarnings = FALSE) # Ordner erstellen, falls nicht vorhanden

  lapply(names(pred_cells_list), function(name) {
    message("Current col: ", name)
    pred_cells <- pred_cells_list[[name]]
    message("before prepare!")
    df_plot <- prepare_sankey_data(sce_cells, pred_cells)
    message("after prepare!")
    output_path <- paste0("~/sankey_plots/sankey_", name, ".png")
    message("before plot creation!")
    create_sankey_plot(df_plot, output_path)
    message("Plot gespeichert: ", output_path)
  })
}

# --- Starte das Skript ---
main()
