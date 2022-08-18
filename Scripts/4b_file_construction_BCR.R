library(ggplot2)
library(shazam)
library(alakazam)

#And now over to the BCR data
BCR_data <- 
    read.csv("All_bracer_results/IMGT_gapped_db_complete_post_clonality.csv")


table(substr(BCR_data$CELL, 1, 6))
#0022_1 0051_2 0147_1 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1
#     2     14      2      8     35    458    202     32    454 

#Now, we add mutational information
BCR_data_full <- observedMutations(BCR_data, sequenceColumn = "SEQUENCE_IMGT",
                                   germlineColumn = "GERMLINE_IMGT",
                                   regionDefinition = IMGT_V)

uniqueCells <- unique(BCR_data_full$CELL)
uniqueCells <- uniqueCells[-which(is.na(uniqueCells))]

interestingCols <- c("Cell", "HamClone", "intraClonalDistance", "Isotype", "Junction_H",
                     "V_gene_H", "TPM_H", "Rep_CDR_mut_H", "Sil_CDR_mut_H", 
                     "Rep_FWR_mut_H", "Sil_FWR_mut_H", "All_mut_H", "Func_IgK", 
                     "V_gene_K", "TPM_K", "Rep_CDR_mut_K", "Sil_CDR_mut_K", 
                     "Rep_FWR_mut_K", "Sil_FWR_mut_K", "All_mut_K", "Func_IgL",  
                     "V_gene_L", "TPM_L", "Rep_CDR_mut_L", "Sil_CDR_mut_L", 
                     "Rep_FWR_mut_L", "Sil_FWR_mut_L", "All_mut_L")  
bcrDf <- data.frame(matrix(nrow = length(uniqueCells), 
                           ncol = length(interestingCols), 
                           dimnames = list(uniqueCells, interestingCols)))
#So now, we run this very dirty script below. 
BCR_H <- BCR_data_full[which(BCR_data_full$LOCUS == "H"),]
for(i in uniqueCells){
    if(any(BCR_H$CELL == i)){
        locCell <- BCR_H[which(BCR_H$CELL == i),]
        locDf <- bcrDf[which(row.names(bcrDf) == i),]
        locDf$Cell <- locCell$CELL
        locDf$HamClone <- locCell$HamClone
        locDf$intraClonalDistance <- locCell$intraClonalDistance
        locDf$Isotype <- locCell$ISOTYPE
        locDf$Junction_H <- locCell$JUNCTION
        locDf$V_gene_H <- locCell$V_CALL
        locDf$TPM_H <- locCell$TPM
        locDf$Rep_CDR_mut_H <- locCell$mu_count_cdr_r
        locDf$Sil_CDR_mut_H <- locCell$mu_count_cdr_s
        locDf$Rep_FWR_mut_H <- locCell$mu_count_fwr_r
        locDf$Sil_FWR_mut_H <- locCell$mu_count_fwr_s
        locDf$All_mut_H <- sum(locCell$mu_count_cdr_r, 
                               locCell$mu_count_cdr_s, 
                               locCell$mu_count_fwr_r, 
                               locCell$mu_count_fwr_s)
        bcrDf[which(row.names(bcrDf) == i),] <- locDf
    }
}

BCR_K <- BCR_data_full[which(BCR_data_full$LOCUS == "K"),]
for(i in uniqueCells){
    if(any(BCR_K$CELL == i)){
        locCell <- BCR_K[which(BCR_K$CELL == i),]
        locDf <- bcrDf[which(row.names(bcrDf) == i),]
        locDf$Func_IgK <- locCell$FUNCTIONAL
        locDf$V_gene_K <- locCell$V_CALL
        locDf$TPM_K <- locCell$TPM
        locDf$Rep_CDR_mut_K <- locCell$mu_count_cdr_r
        locDf$Sil_CDR_mut_K <- locCell$mu_count_cdr_s
        locDf$Rep_FWR_mut_K <- locCell$mu_count_fwr_r
        locDf$Sil_FWR_mut_K <- locCell$mu_count_fwr_s
        locDf$All_mut_K <- sum(locCell$mu_count_cdr_r, 
                               locCell$mu_count_cdr_s, 
                               locCell$mu_count_fwr_r, 
                               locCell$mu_count_fwr_s)
        bcrDf[which(row.names(bcrDf) == i),] <- locDf
    }
}

#This is ugly coding, but today we will accept that. 
BCR_L <- BCR_data_full[which(BCR_data_full$LOCUS == "L"),]
for(i in uniqueCells){
    if(any(BCR_L$CELL == i)){
        locCell <- BCR_L[which(BCR_L$CELL == i),]
        locDf <- bcrDf[which(row.names(bcrDf) == i),]
        locDf$Func_IgL <- locCell$FUNCTIONAL
        locDf$V_gene_L <- locCell$V_CALL
        locDf$TPM_L <- locCell$TPM
        locDf$Rep_CDR_mut_L <- locCell$mu_count_cdr_r
        locDf$Sil_CDR_mut_L <- locCell$mu_count_cdr_s
        locDf$Rep_FWR_mut_L <- locCell$mu_count_fwr_r
        locDf$Sil_FWR_mut_L <- locCell$mu_count_fwr_s
        locDf$All_mut_L <- sum(locCell$mu_count_cdr_r, 
                               locCell$mu_count_cdr_s, 
                               locCell$mu_count_fwr_r, 
                               locCell$mu_count_fwr_s)
        bcrDf[which(row.names(bcrDf) == i),] <- locDf
    }
}


#And with this, we are ready to export the file. 
write.csv(bcrDf, "All_bracer_results/condensed_file_one_BCR_per_row.csv")
