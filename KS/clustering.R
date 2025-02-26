library(data.table)
library(SingleCellExperiment)
library(CATALYST)


print("load data ...")
sce <- readRDS("/nfs/proj/collab_anke_hannover/sce_scaled.rds")


print("clustering ...")
sce <- cluster(sce, features = "type", 
               xdim = 20, ydim = 20, maxK = 50, 
               verbose = FALSE, seed = 1)

print(sce)

print("export files ...")
meta50 <- cluster_codes(sce)$meta50
soms_for_columns <- colData(sce)$cluster_id
columns_meta50 <- meta50[soms_for_columns]

write.csv(data.frame("patient_id" = colData(sce)$patient_id) , "CyTOF_LS/data/patient_id.csv")
write.csv(data.table( "meta50" = columns_meta50) , "CyTOF_LS/data/meta50.csv")
write.csv(assay(sce , "scaled_exprs") , "CyTOF_LS/data/markerXcells.csv") 

print("done.")

