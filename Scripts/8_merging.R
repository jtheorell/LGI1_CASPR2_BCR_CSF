#Here, we will add a lot of BCR and TCR-related information to the experiment
library(SingleCellExperiment)
library(scran)
library(scater)

#First, the import of the data from step 3: 
csfSce <- readRDS("SingleCellExpFiles/csfSce_4_with_prot_and_QUMI.rds")

#Now, we also import the files from step 4. 

#Which led to the loss of 11 cells, as two of the 13 were already excluded. 
#First, we import and add the BCR information. 
BCR_data <- read.csv("All_bracer_results/condensed_file_one_BCR_per_row.csv", 
                     row.names = 1)
#Here, we start by excluding all the cells that are not present in the csfSce
BCR_data_matched <- BCR_data[which(row.names(BCR_data) %in% colnames(csfSce)),]

#We sadly lose 72 cells here (20%)
#Very bad indeed. We will investigate this a bit. 
donorPre <- factor(substr(row.names(BCR_data), 1, 6))
donorPost <- substr(row.names(BCR_data_matched), 1, 6)
donorPost <- factor(donorPost, levels = levels(donorPre))
(table(donorPre)-table(donorPost))/table(donorPre)

#This shows a quite representative loss of sequences all over the board, with 
#between 10 and 20% osses, apart from with the very low numbers. The acute samples
#generally have slightly below average losses, between 10 and 15%.

#We will also reduce the data the other way around, by removing all B cells from the
#csfSce that do not have an associated full BCR
BcsfSce <- csfSce[,which(csfSce$proteinCluster == "B_lin")]
BcsfSceFull <- BcsfSce[,which(colnames(BcsfSce) %in% row.names(BCR_data_matched))]

#This makes us lose another 56 cells. After this, we are down to 477 completely
#complete B cells
csfSce <- cbind(csfSce[,-which(csfSce$proteinCluster == "B_lin")],
                BcsfSceFull)

#Now, we will make the BCR_data_matched file as long as the csfSce is wide, 
#and then we will merge them. 
bigDataFrame <- data.frame(matrix(nrow = ncol(csfSce), 
                                  ncol = ncol(BCR_data_matched),
                                  dimnames = list(colnames(csfSce),
                                                  colnames(BCR_data_matched))))
#And now we fill it with the data.
for(i in row.names(BCR_data_matched)){
    bigDataFrame[which(row.names(bigDataFrame) == i),] <- 
        BCR_data_matched[which(row.names(BCR_data_matched) == i),]
}

#And now, this is added to the big thing
colData(csfSce) <- cbind(colData(csfSce), bigDataFrame)

#########
#TCR ab information
TCR_AB_data <- read.csv("All_tracer_ab_results/complete_Tracer_AB_file_with_clonal_collapses.csv")
table(substr(TCR_AB_data$cell_name, 1,6))

#0022_1 0051_2 0147_1 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1 
#   270    343    157    219    148    230    202    200    410 


#Now, we remove all the rows from the TCR data that we do not have RNA information
#for
TCR_AB_data_matched <- TCR_AB_data[which(TCR_AB_data$cell_name %in% colnames(csfSce)),]
table(substr(TCR_AB_data_matched$cell_name, 1,6))

#0022_1 0051_2 0147_1 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1 
#262    259    145    216    147    288    254    229    426 

#This removed 271 cells (12%), but we do need to do this anyway.

#At this stage, we also include the gamma-deltas
TCR_GD_data <- read.csv("All_tracer_gd_results/complete_Tracer_GD_file_with_clonal_collapses.csv")
table(substr(TCR_GD_data$cell_name, 1,6))

#0022_1 0051_2 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1 147_1_ 
#     4      3      4      2      5      5      4      4      1 
#So these are very few now

#Now, we remove all the rows from the TCR data that we do not have RNA information
#for
TCR_GD_data_matched <- TCR_GD_data[which(TCR_GD_data$cell_name %in% colnames(csfSce)),]

table(substr(TCR_GD_data_matched$cell_name,1,6))
#0022_1 0051_2 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1 
#72     21     47     54     11     35     54     72

#This removed 6 cells. 

#Now, we combine these two files. 
uniqueTCells <- unique(unlist(TCR_AB_data_matched$cell_name, TCR_GD_data_matched$cell_name))
#There is no overlap between these vectors now, as we have focused on the complete
#TCR cells
TCR_data_full <- data.frame("cell_name" = colnames(csfSce), 
                            "A1" = NA, "A2" = NA, "B1" = NA,
                            "B2" = NA, "AB_clone" = NA, "D1" = NA, "D2"= NA,
                            "G1" = NA, "G2" = NA, "GD_clone" = NA)
colnames(TCR_data_full)[c(2:5,7:10)] <- paste0(substr(colnames(TCR_data_full)[c(2:5,7:10)], 1,1),
                                       rep(c("_productive", "_unproductive"), 
                                           times = 4))
row.names(TCR_data_full) <- TCR_data_full$cell_name
TCR_data_full <- TCR_data_full[,-1]

for(i in uniqueTCells){
    locRow <- which(row.names(TCR_data_full) == i)
    if(length(locRow) == 1){
        locFull <- TCR_data_full[locRow,]
        
        locABDat <- TCR_AB_data_matched[which(TCR_AB_data_matched$cell_name == i),]
        locFull[,c("A_productive", "A_unproductive",
                   "B_productive", "B_unproductive", "AB_clone")] <- 
            locABDat[,c("A_productive", "A_unproductive",
                        "B_productive", "B_unproductive", "collapsedClone")]
        locGDDat <- TCR_GD_data_matched[which(TCR_GD_data_matched$cell_name == i),]
        locFull[,c("D_productive", "D_unproductive",
                   "G_productive", "G_unproductive", "GD_clone")] <- 
            locGDDat[,c("D_productive", "D_unproductive",
                        "G_productive", "G_unproductive", "collapsedClone")]
        TCR_data_full[locRow,] <- locFull
    }
}

#Now, wwe exclude the T-cells that lack a TCR from the csfSce
TcsfSce <- csfSce[,which(csfSce$proteinCluster %in% c("CD8T", "CD4T", "DPT", "DNT"))]
TcsfSceFull <- TcsfSce[,which(colnames(TcsfSce) %in% uniqueTCells)]

#569 cells are lost here, so we are down to 1901 T cells. 
csfSce <- cbind(csfSce[,-which(csfSce$proteinCluster %in% c("CD8T", "CD4T", "DPT", "DNT"))],
                TcsfSceFull)

TCR_data_complete <- TCR_data_full[which(row.names(TCR_data_full) %in% colnames(csfSce)),]

#And now we combine these 
identical(row.names(TCR_data_complete), colnames(csfSce))
#FALSE

#So we reorder them both. 
TCR_data_complete <- TCR_data_complete[order(row.names(TCR_data_complete)),]
csfSce <- csfSce[,order(colnames(csfSce))]
identical(row.names(TCR_data_complete), colnames(csfSce))

#And with that, we can combine them: 
colData(csfSce) <- cbind(colData(csfSce), TCR_data_complete)

#Now, we check if there are cells with a TCR or a BCR that are not categorized 
#accordingly in the transcriptome/protein data or if there are doublets. 
noBCR <- is.na(csfSce$Junction_H)
table(noBCR, csfSce$proteinCluster)
#noBCR  B_lin CD4T CD8T  DNT  DPT   NK
# TRUE      0 1273  493   26  104   96
# FALSE   477    4    1    0    0    1

#So a few cells need to be excluded here. 
csfSce <- csfSce[,-which(anyBCR & csfSce$proteinCluster != "B_lin")]

#Now, over to the T-cell side
noTCR <- apply(colData(csfSce)[34:43], 1, function(x) all(is.na(x)))

table(noTCR, csfSce$proteinCluster)

#Once again, 7 cells are on the wrong side and will be killed

csfSce <- csfSce[,-which(noTCR == FALSE & csfSce$proteinCluster %in% c("B_lin", "NK"))]


noTCR <- apply(colData(csfSce)[34:43], 1, function(x) all(is.na(x)))
noBCR <- is.na(csfSce$Junction_H)
table(noTCR, csfSce$proteinCluster)
table(noBCR, csfSce$proteinCluster)
table(noTCR, noBCR)

#And here we are, with 1896 T cells with a TCR, 471 B cells with a BCR and 95 NK cells
#without any of them! Time to save. 

saveRDS(csfSce, "SingleCellExpFiles/csfSce_5_with_BCR_TCR.rds")
