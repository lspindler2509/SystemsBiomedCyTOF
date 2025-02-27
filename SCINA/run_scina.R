library(SCINA)
library(SingleCellExperiment)
library(data.table)

# functions
signatures_from_classification_matrix <- function(classification_matrix) {
  class_mat <- as.data.frame(t(classification_matrix))
  
  mat_list <- as.list(class_mat)
  
  list_high <- lapply(mat_list, function(x) rownames(class_mat)[x == '1'])
  list_low <- lapply(mat_list, function(x) rownames(class_mat)[x == '-1'])
  
  
  list_low <- lapply(list_low, function(x) {
    if(!identical(x, character(0))) {
      paste0("low_", x)
    }
  })
  
  signatures <- lapply(colnames(class_mat), function(x) c(list_high[[x]], list_low[[x]]))
  names(signatures) <- colnames(class_mat)
  
  return(signatures)
}



sce <- readRDS("data/sce.rds")
sce_scaled <- readRDS("data/sce_scaled.rds")

# transformed counts
exp <- assay(sce, "exprs")
exp_scaled <- assay(sce_scaled, "scaled_exprs")



orig_markers <- read.csv("data/cell_types_original_markers.csv", row.names=1)
colnames(orig_markers) <- sub("\\.", "-", colnames(orig_markers))
signs_orig_markers <- signatures_from_classification_matrix(orig_markers)

main_markers <- read.csv("data/cell_types_main_markers.csv", row.names=1)
colnames(main_markers) <- sub("\\.", "-", colnames(main_markers))
signs_main_markers <- signatures_from_classification_matrix(main_markers)

print(signs_orig_markers)



func <- function(exp_ma, signs) {
  print("scina start")
  s <- Sys.time()
  results = SCINA(exp_ma, signs, rm_overlap=FALSE)
  print("done after (%.2f)", Sys.time()-s)
  
  return(results)
}

exp_list <- list("exp" = exp, "exp_scaled" = exp_scaled)
results_orig <- lapply(exp_list, func, signs_orig_markers)

results_main_scaled <- func(exp_list$exp_scaled, signs_main_markers)
fwrite(list(results_main_scaled$cell_labels), file = "results_new/scina_sce_scaled_main_markers.csv")
write.table(results_main_scaled$probabilities, "results_new/scina_scaled_main_probabilities.csv", quote=FALSE, sep=",", col.names = FALSE)


print("writing to files")

fwrite(list(results_orig$exp$cell_labels), file = "results_new/scina_sce_unscaled_original_markers.csv")
fwrite(list(results_orig$exp_scaled$cell_labels), file = "results_new/scina_sce_scaled_original_markers.csv")

write.table(results_orig$exp$probabilities, "results_new/scina_unscaled_original_probabilities.csv", quote=FALSE, sep=",", col.names = FALSE)
write.table(results_orig$exp_scaled$probabilities, "results_new/scina_scaled_original_probabilities.csv", quote=FALSE, sep=",", col.names = FALSE)





print("Done")


