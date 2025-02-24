library(data.table)
library(tidyr)
library(dplyr)
library(pheatmap)
library(RColorBrewer)


filepath <- "create_heatmaps/metrics_all_against_all_original_with_weighted_no_intersection.csv" # file with computed metric scores
metrics_all <- fread(filepath)

# modify metrics df
metrics_all2 <- as.data.frame(metrics_all) %>%
  separate(method1, into=c("tool1", NA, "data1", "matrix1"), sep="_") %>%
  separate(method2, into=c("tool2", NA, "data2", "matrix2"), sep="_") %>%
  unite("data_matrix", c(data1, matrix1), sep = "_") %>%
  unite("data_matrix2", c(data2, matrix2), sep = "_") %>%
  filter(data_matrix == data_matrix2) %>%
  select(- data_matrix2) %>%
  unique()

create_heatmap_for_metric <- function(filt, df, metric, dir, clust=FALSE) {
  # define legend range and colors by the possible range of the metric
  if (metric %in% c("ACC", "F1", "F1_weighted",  "ACC_weighted")) {
    cols <- colorRampPalette(c("white", "red", "darkred"))(1000)
    breaks <- seq(0,1,length.out=1001);
  } else{
    cols <- colorRampPalette(c("darkblue", "blue", "white", "red", "darkred"))(1000)
    breaks <- seq(-1,1,length.out=1001)
  }
  
  # filter the metric df to only contain one set of input parameters (filt) and one metric
  metrics <- df %>%
    filter(data_matrix == filt) %>%
    select(tool1, tool2, matches(paste0("^",metric, "$"))) %>%
    unique() %>%
    as.data.table() %>%
    dcast(formula = ... ~tool1, value.var = metric) %>%
    tibble::column_to_rownames(var="tool2")
  
  # create the heatmap
  png(paste0(dir, "/", paste(metric, filt, sep = "_"), ".png"))
  pheatmap(as.matrix(metrics), cluster_rows = clust, cluster_cols = clust, color=cols, breaks = breaks, clustering_method = "ward.D2", display_numbers = TRUE, number_format = "%.2f", number_color = "white", fontsize_number = 16, fontsize_col = 20, fontsize_row = 20)
  # main = paste(metric, gsub("_", " ", filt))
  dev.off()
}

metric_list <- c("ACC", "F1", "ARI", "kappa", "F1_weighted", "ACC_weighted")
metric_list <- c("F1", "ARI", "kappa", "F1_weighted", "ACC_weighted")

# create output directories
dir.create(file.path("create_heatmaps/heatmaps_new_clustered"))
dir.create(file.path("create_heatmaps/heatmaps_new"))

# save heatmaps (clustered and not clustered) in files
lapply(metric_list, function(metric) lapply(unique(metrics_all2$data_matrix), create_heatmap_for_metric, metrics_all2, metric, "create_heatmaps/heatmaps_new_clustered", TRUE))
lapply(metric_list, function(metric) lapply(unique(metrics_all2$data_matrix), create_heatmap_for_metric, metrics_all2, metric, "create_heatmaps/heatmaps_new"))
