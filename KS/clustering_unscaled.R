library(data.table)
library(SingleCellExperiment)
library(CATALYST)


print("load data ...")
sce <- readRDS("CyTOF_LS/Manual_Gatin_Noemi/sce.rds")

print("clustering ...")
sce <- cluster(sce, features = "type", 
               xdim = 20, ydim = 20, maxK = 300, 
               verbose = FALSE, seed = 1)

print(sce)



print("export files ...")
meta300 <- cluster_codes(sce)$meta300
soms_for_columns <- colData(sce)$cluster_id
columns_meta300 <- meta300[soms_for_columns]

write.csv(data.frame("patient_id" = colData(sce)$patient_id) , "CyTOF_LS/data/300/patient_id_unscaled.csv")
write.csv(data.table( "meta300" = columns_meta300) , "CyTOF_LS/data/300/meta300_unscaled.csv")
write.csv(assay(sce , "exprs") , "CyTOF_LS/data/300/markerXcells_unscaled.csv")
write.table(DataFrame("manual_gating" = colData(sce)$cell_type) , "CyTOF_LS/data/300/manual_gating_unscaled.csv" , quote = F , row.names = F ,sep = "\t")

print("done.")