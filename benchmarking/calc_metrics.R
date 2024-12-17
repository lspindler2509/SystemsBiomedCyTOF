library(tibble)
library(data.table)

library(fossil) # ARI
library(yardstick) # balanced accuracy
library(MLmetrics) # F1 score



# input parameters

sce_path <- "data/sce_cleaned.rds"     # path to sce file
mapper_path <- "cell_label_mapper.csv"     # file with uniform cell type labels as mapper
folder_path <- "results"     # folder with the result file(s) for the different tools/methods/etc.
file_pattern <- "*.csv"     # pattern of the wanted file(s) (! not the complete filename)

out_path <- "metrics.csv"     # file to save the results of the metrics



# functions

# mapping the cell labels to the uniform labels
trans_cell_types <- function(vec, dt, col) {
  trans_vec <- dt$label[ match(vec, dt[[col]]) ]
  trans_vec <- factor(trans_vec, levels = unique(dt$label)[unique(dt$label) != "conventional CD4 cells"]) # without filtering balanced accuracy can't be computed
  return(trans_vec)
}

# calculating the metrics
calc_metrics <- function(truth, estimate, method) {
  # adjusted rand index
  ari <- adj.rand.index(as.numeric(truth), as.numeric(estimate))
  # balanced accuracy
  acc <- bal_accuracy_vec(truth, estimate)
  # F1 score
  f1 <- F1_Score(truth, estimate)
  
  return(data.table("method"=method, "ARI"=ari, "ACC"=acc, "F1"=f1))
}



# extracting the result vectors from the file(s)
filenames <- list.files(folder_path, pattern=file_pattern, full.names=TRUE)
results_list <- lapply(filenames, scan, character(), quote="", sep="\n")
names(results_list) <- unlist(lapply(filenames, function(file) sub(file_pattern, '', basename(file))))


# translating the cell labels to uniform names
sce <- readRDS(sce_path)
cell_type_dt <- fread(mapper_path, na.strings = c(NA_character_, "")) # reading in the mapper data table

sce_cells <- trans_cell_types(sce$cell_type, cell_type_dt, "sce_label")
pred_cells <- lapply(results_list, trans_cell_types, cell_type_dt, "matrix_label")


# calculating the metrics (adjusted rand index, balanced accuracy, F1 score)
metrics <- lapply(names(pred_cells), function(method) calc_metrics(sce_cells, pred_cells[[method]], method))
comp_metrics <- rbindlist(metrics)


# saving the metrics results in a file
fwrite(comp_metrics, file = out_path)
