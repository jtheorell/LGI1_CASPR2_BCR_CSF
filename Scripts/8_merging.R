#Here, we will add a lot of BCR and TCR-related information to the experiment
library(SingleCellExperiment)
library(scran)
library(scater)

#First, the import of the data from previously
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_GLMPCA_UMAP.rds")
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]
#We start by removing the few cells that have insufficiently useful transcriptomes
BCR_H <- BCR_H[which(BCR_H$CELL %in% colnames(csfSce)),]
BCR_H_ordered <- BCR_H[match(colnames(csfSce), BCR_H$CELL),]

identical(BCR_H_ordered$CELL, colnames(csfSce))
#TRUE

#So now, we can start integrating straight. 
colData(csfSce) <- cbind(colData(csfSce), BCR_H_ordered[,c("V_CALL","JUNCTION",
                                                          "JUNCTION_LENGTH",
                                                          "ISOTYPE","HamClone",
                                                          "light_type","Clonal",
                                                          "All_mutations",
                                                          "Non_silent_mutations",
                                                          "Specific", "SubClone")])

saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_5_BCR.rds")

