library(SingleCellExperiment)
library(scater)
library(sceasy)
library(reticulate)
library(zellkonverter)
library(dplyr)


# Load the .rds file
sce <- readRDS("~/Systems_biomedicine/acdc/data/Manual_Gating_Noemi/sce_scaled.rds")

# General structure of the object
sce


# Check the metadata
colData(sce)

# Check the scaled expression matrix
assayNames(sce)    # available assays
counts <- assay(sce, "counts")  
head(counts)

expr <- assay(sce, "exprs") 
head(expr)

scaled_expr <- assay(sce, "scaled_exprs") 
head(scaled_expr)

# Check the row (gene/marker) and column (cell) annotations
rowData(sce)  # Marker information
colData(sce)  # Cell metadata

# Summary of the dataset
summary(sce)


# Save as h5ad using zellkonverter
zellkonverter::writeH5AD(
  sce,
  "~/Systems_biomedicine/acdc/data/Manual_Gating_Noemi/sce_scaled.h5ad"
)




# Convert colData to a dataframe
metadata <- as.data.frame(colData(sce))

# Overview of unique values in key columns
cat("Number of unique patients:", length(unique(metadata$patient_id)), "\n") #20
cat("Number of unique timepoints:", length(unique(metadata$timepoint)), "\n") #3
cat("Number of unique groups:", length(unique(metadata$group)), "\n") # 6

# Summarize counts per patient, timepoint, and group
patient_summary <- metadata %>% count(patient_id)
timepoint_summary <- metadata %>% count(timepoint)
group_summary <- metadata %>% count(group)

# Print summaries
print("Patients overview:")
print(patient_summary)

print("Timepoints overview:")
print(timepoint_summary)

print("Groups overview:")
print(group_summary)


celltype_summary <- metadata %>% count(cell_type)
print("Cell types overview:")
print(celltype_summary)

