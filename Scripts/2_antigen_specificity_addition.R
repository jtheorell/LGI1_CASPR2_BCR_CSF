#As this is such a crucial thing, and as the selection for specificity
#testing is a separate script between, it will be added here as an individual
#file, despite its brevity. 
BCR_all <- read.csv("Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv")

BCR_all$Specific <- "Not_tested"

#This file encompasses both the cells that were actually tested and the clonal
#relatives of these that have identical sequences, where the antigen specificity 
#has been inferred. This might theoretically not be justified
#in all cases, but it has been in all that have been tested.

specInfo <- read.csv("Data/BCR_auxiliaries/AntigenSpecificInfo.csv")

BCR_all$Specific[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$Specific[which(BCR_all$CELL %in% 
                           specInfo$CELL[which(specInfo$Specific == TRUE)])] <- 
    TRUE

table(BCR_all$Specific)
#FALSE Not_tested       TRUE 
#93        459        231

#For some cells, the whole clone might not be present in the above file, so 
#a further attempt to include all cells in the clones are added here: 
posClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "TRUE")])
negClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "FALSE")])
BCR_all$Specific[which(BCR_all$HamClone %in% posClones)] <- TRUE
BCR_all$Specific[which(BCR_all$HamClone %in% negClones)] <- FALSE

table(BCR_all$Specific)
#FALSE Not_tested       TRUE 
#100        439        244

#This looks right. 
write.csv(BCR_all, "Data/BCR_database_versions/6_Specificity_included.csv", 
          row.names = FALSE)

