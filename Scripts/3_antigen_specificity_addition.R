#As this is such a crucial thing, and as the selection for specificity
#testing is a separate script between, it will be added here as an individual
#file, despite its brevity. 
BCR_all <- read.csv("Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv")

BCR_all$Specific <- "Not_tested"

#This file encompasses both the cells that were actually tested and the clonal
#relatives of these that have identical sequences, where the antigen specificity has been inferred. 
specInfo <- read.csv("Data/BCR_auxiliaries/AntigenSpecificInfo.csv")

BCR_all$Specific[which(BCR_all$CELL %in% specInfo$CELL)] <- FALSE

BCR_all$Specific[which(BCR_all$CELL %in% 
                           specInfo$CELL[which(specInfo$Specific == TRUE)])] <- 
    TRUE

table(BCR_all$Specific)
#FALSE Not_tested       TRUE 
#87        422        237

#This looks right. 
write.csv(BCR_all, "Data/BCR_database_versions/6_Specificity_included.csv", 
          row.names = FALSE)

