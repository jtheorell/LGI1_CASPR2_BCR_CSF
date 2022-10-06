#As this is such a crucial thing, and as the selection for specificity
#testing is a separate script between, it will be added here as an individual
#file, despite its brevity. 
BCR_all <- read.csv("Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv")

BCR_all$Specific <- BCR_all$Specific_UCA <- "Not_tested"

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
#92        459        232

#For some cells, the whole clone might not be present in the above file, so 
#a further attempt to include all cells in the clones are added here: 
posClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "TRUE")])
negClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "FALSE")])
BCR_all$Specific[which(BCR_all$HamClone %in% posClones)] <- TRUE
BCR_all$Specific[which(BCR_all$HamClone %in% negClones)] <- FALSE

table(BCR_all$Specific)
#FALSE Not_tested       TRUE 
#92        439        252

##############
#POTENTIALLY FOR FIGURE 1
##############
table(BCR_all$Sample[which(BCR_all$LOCUS == "H")], 
      BCR_all$Specific[which(BCR_all$LOCUS == "H")],
      BCR_all$Clonal[which(BCR_all$LOCUS == "H")])
#Not clonal
#        FALSE Not_tested TRUE
#1166_1     7        110   19
#1227_1     2         44    5
#1284_1     8         49    1

#Clonal
#        FALSE Not_tested TRUE
#1166_1    13          4   59
#1227_1     2          5    7
#1284_1    13          3   30
##############


#Now the exact same thing is performed for the unmutated common ancestors

BCR_all$Specific_UCA[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$Specific_UCA[which(BCR_all$CELL %in% 
                               specInfo$CELL[which(specInfo$Specific_UCA == TRUE)])] <- 
    TRUE

table(BCR_all$Specific_UCA)
#FALSE Not_tested       TRUE 
#278        459         46

#For some cells, the whole clone might not be present in the above file, so 
#a further attempt to include all cells in the clones are added here: 
posClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific_UCA == "TRUE")])
negClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific_UCA == "FALSE")])
#Futher, as the UCAs have been generated with more than one method for some
#cells, we identify a cell with one specific UCA as having a specific UCA, 
#as a negative test result is more likely to be wrong for technical reasons. 
negClones <- negClones[-which(negClones %in% posClones)]

BCR_all$Specific_UCA[which(BCR_all$HamClone %in% posClones)] <- TRUE
BCR_all$Specific_UCA[which(BCR_all$HamClone %in% negClones)] <- FALSE

table(BCR_all$Specific_UCA)
#FALSE Not_tested       TRUE 
#254        455         74

#This looks right. 
write.csv(BCR_all, "Data/BCR_database_versions/6_Specificity_included.csv", 
          row.names = FALSE)

