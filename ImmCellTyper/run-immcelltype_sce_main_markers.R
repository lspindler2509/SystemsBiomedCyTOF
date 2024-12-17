library(ImmCellTyper)
library(CATALYST)
library(readxl)
library(flowCore)
library(SingleCellExperiment)
library(cowplot)
library(ggplot2)
library(reshape2)
library(dplyr)
library(caret)
library(rstatix)
library(Rphenograph)
library(readr)
library(RColorBrewer)
library(igraph)
library(pheatmap)
library(CytofRUV)

count <- 1

run_kmeans2 <- function (data) 
{
  k.means <- kmeans(data, 2, iter.max = 500)
  pos <- which(k.means$centers == max(k.means$centers))
  neg <- which(k.means$centers == min(k.means$centers))
  print(pos)
  print(neg)
  print(k.means$size)
  results <- k.means$cluster

  print(count)
  results[results == pos] <- "+"
  results[results == neg] <- "-"

  plot_data <- data.frame(
    Index = seq_along(data),
    Value = data,
    Cluster = factor(results)
  )

  directory <- "systems_biomed/plots/sce/eval_raw_main_markers/"

  p <- ggplot(plot_data, aes(x = Value)) +
  geom_histogram(aes(fill = Cluster), position = "dodge", bins = 30, alpha = 0.6) +
  labs(
    title = paste("K-means Clustering for", count),
    x = "Value",
    y = "Density",
    fill = "Cluster"
  ) +
  theme_minimal()
  ggsave(paste0(directory, "pos_neg_colored",count,".png"), plot = p, bg = "white")

  p <- ggplot(plot_data, aes(x = Value)) +
  geom_density() +
  labs(
    title = paste("K-means Clustering for", count),
    x = "Value",
    y = "Density",
    fill = "Cluster"
  ) +
  theme_minimal()
  ggsave(paste0(directory, "pos_neg_",count,".png"), plot = p, bg = "white")

  assign("count", count + 1, envir = .GlobalEnv)
  return(results)
}

binaryClass2 <- function (data = NULL, class.file = NULL, scale = FALSE) 
{
  if (is.null(data)) {
    return("Please provide data. Use load_data")
  }
  if (is.null(class.file)) {
    return("Please provide a classification file.")
  }
  else {
    data <- data[, colSums(data != 0, na.rm = TRUE) > 0]
    if (scale == TRUE) {
      data <- data.frame(scale(data))
    }
    else if (scale == FALSE) {
      data <- data
    }
    else {
      return("Please choose either TRUE or FALSE for scale.")
    }
    kmeans.results <- data.frame(apply(data, 2, run_kmeans2))
    colnames(kmeans.results) <- colnames(data)
    kmeans.limits <- data.frame(matrix(nrow = 0, ncol = 3))
    colnames(kmeans.limits) <- c("pos_min", "Q1", "Q3")
    for (marker in colnames(kmeans.results)) {
      if (length(unique(kmeans.results[, marker])) == 
        1 & unique(kmeans.results[, marker])[1] == "-") {
        pos_min <- 999
        Q1 <- 0
        Q3 <- 0
      }
      else {
        pos_min <- as.numeric(min(data[kmeans.results[, 
          marker] == "+", marker]))
        Q1 <- as.numeric(quantile(data[kmeans.results[, 
          marker] == "+", marker], 0.25))
        Q3 <- as.numeric(quantile(data[kmeans.results[, 
          marker] == "+", marker], 0.75))
      }
      kmeans.results[, marker] <- as.character(kmeans.results[, 
        marker])
      kmeans.results[data[, marker] <= Q1 & kmeans.results[, 
        marker] == "+", marker] <- "L"
      kmeans.results[data[, marker] > Q1 & data[, marker] <= 
        Q3, marker] <- "I"
      kmeans.results[data[, marker] > Q3, marker] <- "H"
      kmeans.limits[marker, ] <- c(pos_min, Q1, Q3)
    }
    kmeans.results[, "Number"] <- rownames(kmeans.results)
    class.results <- data.frame(matrix(nrow = dim(data)[1], 
      ncol = 1))
    colnames(class.results) <- c("Cell.Type")
    class.cat <- read.table(class.file, header = TRUE, sep = ",", 
      row.names = 1, check.names = FALSE)
    class.cat <- apply(class.cat, c(1, 2), function(x) {
      gsub(" ", "", x)
    })
    class.cat[class.cat == ""] <- "A"
    for (cell.type in rownames(class.cat)) {
      temp.row <- class.cat[cell.type, , drop = FALSE]
      temp.row <- temp.row[, temp.row != "A", drop = FALSE]
      temp.row.fine <- temp.row[, temp.row != "+", drop = FALSE]
      temp.row.all <- temp.row[, temp.row == "+", drop = FALSE]
      rows.required <- merge(kmeans.results, temp.row.fine)
      if (length(colnames(temp.row.all)) == 1) {
        rows.required <- rows.required[rows.required[, 
          colnames(temp.row.all)] %in% c("L", "I", "H"), 
          ][, "Number"]
      }
      else {
        rows.required <- rows.required[apply(rows.required[, 
          colnames(temp.row.all)], 1, function(x) {
          all(x %in% c("L", "I", "H"))
        }), ][, "Number"]
      }
      class.results[rows.required, "Cell.Type"] <- cell.type
    }
    class.results[is.na(class.results[, "Cell.Type"]), "Cell.Type"] <- "Unclassified"
    return(class.results)
  }
}


sce <- readRDS("/nfs/proj/collab_anke_hannover/Manual_Gating_Noemi/sce.rds")
class_dir <- 'classification_matrix_main_markers.csv'
types <- read_csv(class_dir)
print(types)
  
name <- "unscaled"
directory <- "systems_biomed/plots/sce/eval_results_raw_main_markers/"
# Schritte für jede Iteration
plot <- displayMarkers(sce, types, color_by = 'sample_id', ncol = 3)
ggsave(paste0(directory,"displayMarkers_excl_", name, ".png"), plot = plot, create.dir = TRUE)
print(rownames(sce))

print(assayNames(sce))
exprs_matrix <- assay(sce, 'exprs')
exprs_matrix[exprs_matrix < 0] <- 0
assay(sce, 'exprs') <- exprs_matrix
exprs_matrix <- NULL

exprs <- t(assay(sce, 'exprs'))
binary.results <- binaryClass2(exprs, class_dir)
print(class(binary.results))
write.csv(binary.results, file = "results_sce_main_markers.csv", row.names = FALSE)
sce$cluster_id <- unlist(binary.results)

plot <- plotbcFreq(sce, binary.results)
ggsave(paste0(directory, "plotbcFreq_", name, ".png"), plot = plot, bg = "white")


plot <- plotbcHeatmap(sce, binary.results, remove.unclass = F)
ggsave(paste0(directory, "plotbcHeatmap_", name, ".png"), plot = plot)

set.seed(1234)
sce <- runDR(sce, dr = "UMAP", cells = 2000, features = "type")
sce <- runDR(sce, dr = "TSNE", cells = 2000, features = "type")

plot <- plotDR(sce, "UMAP", color_by = "cluster_id")
ggsave(paste0(directory, "plotDRumap_", name, ".png"), plot = plot, bg = "white")

plot <- plotDR(sce, "TSNE", color_by = "cluster_id")
ggsave(paste0(directory, "plotDRtsne_", name, ".png"), plot = plot, bg = "white")

plot <- plotDR(sce, "UMAP", color_by = "cluster_id", facet_by = 'condition')
ggsave(paste0(directory, "plotDRUMAP_facet_", name, ".png"), plot = plot, bg = "white")

plot <- plotDR(sce, "TSNE", color_by = "cluster_id", facet_by = 'condition')
ggsave(paste0(directory, "plotDRTSNE_facet_", name, ".png"), plot = plot, bg = "white")

plot <- plotFreq(sce, type = 'stacked', group_by = 'condition')
ggsave(paste0(directory, "stacked_plot_", name, ".png"), plot = plot, bg = "white")

plot <- plotFreq(sce, type = 'box', group_by = 'condition')
ggsave(paste0(directory, "boxplot_", name, ".png"), plot = plot, bg = "white")

plot <- plotStateMarkerExprs(sce, group = 'condition')
ggsave(paste0(directory, "plotStateMarkerExprs_", name, ".png"), plot = plot, bg = "white")