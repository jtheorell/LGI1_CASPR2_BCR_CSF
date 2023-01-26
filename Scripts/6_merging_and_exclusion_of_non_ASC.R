#Here, we will add a lot of BCR and TCR-related information to the experiment
library(SingleCellExperiment)
library(scran)
library(scater)

#First, the import of the data from previously
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]

#For these downsteam analyses, we of course need all data to be present for all
#cells. For this reason, all cells not having a BCR are excluded here. We will
#keep the CD8 and the other B-cells for differing transcriptome analyses downstream
#however. 
csfSceBCR <- csfSce[,which(colnames(csfSce) %in% BCR_H$CELL)]
BCR_H <- BCR_H[which(BCR_H$CELL %in% colnames(csfSce)),]
BCR_H_ordered <- BCR_H[match(colnames(csfSceBCR), BCR_H$CELL),]

identical(BCR_H_ordered$CELL, colnames(csfSceBCR))
#TRUE

#So now, we can start integrating straight. 
colData(csfSceBCR) <- cbind(colData(csfSceBCR), BCR_H_ordered[,c("V_CALL","JUNCTION",
                                                          "JUNCTION_LENGTH",
                                                          "ISOTYPE","HamClone",
                                                          "light_type","Clonal",
                                                          "All_mutations",
                                                          "Non_silent_mutations",
                                                          "Specific", "LRR", "EPTP",
                                                          "Specific_UCA",
                                                          "SubClone")])

#Here, we evalueate the cell type 
table(csfSceBCR$clustersLouvain)
#ASC    B 
#269   82   

saveRDS(csfSceBCR, "Data/SingleCellExpFiles/csfSce_3_BCR_containers.rds")

#Now, we create a massive file, containing all cells, regardless of if they have 
#the above markers or not. 
csfSceNonBCR <- csfSce[,-which(colnames(csfSce) %in% BCR_H$CELL)]

colDatMat <- matrix(NA, ncol(csfSceNonBCR), 14)
colnames(colDatMat) <- c("V_CALL","JUNCTION",
                         "JUNCTION_LENGTH",
                         "ISOTYPE","HamClone",
                         "light_type","Clonal",
                         "All_mutations",
                         "Non_silent_mutations",
                         "Specific", "LRR", "EPTP",
                         "Specific_UCA",
                         "SubClone")
colData(csfSceNonBCR) <- cbind(colData(csfSceNonBCR), colDatMat)

totalCsfSce <- cbind(csfSceBCR, csfSceNonBCR)
#Here, we remove the two CD8T cells with BCRs: 
csfSceNoBCRCD8 <- csfSce[,which(colnames(csfSce) %in% colnames(totalCsfSce))]
totalCsfSce <- totalCsfSce[,match(colnames(csfSceNoBCRCD8), colnames(totalCsfSce))]
identical(colnames(totalCsfSce), colnames(csfSceNoBCRCD8))
#TRUE
saveRDS(totalCsfSce, "Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")

#Here, we also export a file with the B-lineage cell names, to be used for RNA velocity
write.csv(colnames(totalCsfSce)[-which(totalCsfSce$clustersLouvain == "CD8T")], 
          "../External/Data/All_full_cells.csv", row.names = FALSE)

#Sadly, as can be viewed with the figures in 6b, there are essentially no non-ASCs that are
#positive. 
table(csfSceBCR$clustersLouvain, csfSceBCR$Specific)

#    FALSE Not_tested TRUE
#ASC    18        132  119
#B      12         66    4

#These are distributed among donors 1166 (5 cells) and 1284 (2 cells).
#As the B cells are so clearly skewed towards non-specificity, we will exclude them
#from any further transcriptome analyses here, as we will otherwise get a 
#result downstream just showing that ASCs genes are overrepresented among the specifics
#which we already know now. 

csfSce <- csfSceBCR[,-which(csfSceBCR$clustersLouvain == "B")]

saveRDS(csfSce, "Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")
