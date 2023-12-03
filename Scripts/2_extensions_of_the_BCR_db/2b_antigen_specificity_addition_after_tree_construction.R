#As this is such a crucial thing, and as the selection for specificity
#testing is a separate script between, it will be added here as an individual
#file, despite its brevity. 
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

BCR_all$Specific <- BCR_all$Specific_UCA <- BCR_all$EPTP <- BCR_all$LRR <- "Not_tested"

#This file encompasses both the cells that were actually tested and the clonal
#relatives of these that have identical sequences, where the antigen specificity 
#has been inferred. This might theoretically not be justified
#in all cases, but it has been in all that have been tested.

specInfo <- read.csv("Data/BCR_auxiliaries/AntigenSpecificInfo.csv")

BCR_all$Specific[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$Specific[which(BCR_all$CELL %in% 
                           specInfo$CELL[which(specInfo$Specific == TRUE)])] <- 
    TRUE

#A few cells have also been reclassified since this phase:
BCR_all$Specific[which(BCR_all$HamClone == 75)] <- TRUE
BCR_all$Specific[which(BCR_all$HamClone == 91)] <- FALSE


table(BCR_all$Specific)
#FALSE Not_tested       TRUE 
#64        430        268

##For some cells, the whole clone might not be present in the above file, so 
##a further attempt to include all cells in the clones are added here: 
#posClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "TRUE")])
#negClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific == "FALSE")])
#BCR_all$Specific[which(BCR_all$HamClone %in% posClones)] <- TRUE
#BCR_all$Specific[which(BCR_all$HamClone %in% negClones)] <- FALSE
#
#table(BCR_all$Specific)
##FALSE Not_tested       TRUE 
##80        439        264

#Now the exact same thing is performed for the unmutated common ancestors

BCR_all$Specific_UCA[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$Specific_UCA[which(BCR_all$CELL %in% 
                           specInfo$CELL[which(specInfo$Specific_UCA == TRUE)])] <- 
    TRUE

table(BCR_all$Specific_UCA)
#FALSE Not_tested       TRUE 
#262        430         70

##For some cells, the whole clone might not be present in the above file, so 
##a further attempt to include all cells in the clones are added here: 
#posClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific_UCA == "TRUE")])
#negClones <- unique(BCR_all$HamClone[which(BCR_all$Clonal & BCR_all$Specific_UCA == "FALSE")])
##Futher, as the UCAs have been generated with more than one method for some
##cells, we identify a cell with one specific UCA as having a specific UCA, 
##as a negative test result is more likely to be wrong for technical reasons. 
##This is actually deprecated, but we stick with it anyway. 
#negClones <- negClones[-which(negClones %in% posClones)]
#
#BCR_all$Specific_UCA[which(BCR_all$HamClone %in% posClones)] <- TRUE
#BCR_all$Specific_UCA[which(BCR_all$HamClone %in% negClones)] <- FALSE
#
#table(BCR_all$Specific_UCA)
##FALSE Not_tested       TRUE 
##270        439         74 

#And the LRR and EPTP of course. As can be noted above, there is no update on the numbers of cells
#currently, so we will only add the information. 

BCR_all$LRR[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$LRR[which(BCR_all$CELL %in% 
                           specInfo$CELL[which(specInfo$LRR == TRUE)])] <- 
    TRUE

table(BCR_all$LRR)
#FALSE Not_tested       TRUE 
#226        430        106
BCR_all$EPTP[which(BCR_all$CELL%in% specInfo$CELL)] <- FALSE

BCR_all$EPTP[which(BCR_all$CELL %in% 
                      specInfo$CELL[which(specInfo$EPTP == TRUE)])] <- 
    TRUE

table(BCR_all$EPTP)
#FALSE Not_tested       TRUE 
#271        439         73
table(BCR_all$EPTP, BCR_all$LRR)
#           FALSE Not_tested TRUE
#FALSE        154          0  106
#Not_tested     0        430    0
#TRUE          72          0    0
#Great, no overlaps. 



#Now, one thing that turned out to be important in the review process
#was the more clear division of the LGI1 and CASPR2 cells/BCRs, so here
#we introduce a new colData column
BCR_all$Specific_antigen <- BCR_all$Specific
BCR_all$Specific_antigen[which(BCR_all$Specific == "TRUE")] <- "LGI1"
BCR_all$Specific_antigen[which(BCR_all$Specific == "TRUE" &
                                   BCR_all$Sample == "1284_1")] <- "CASPR2"
BCR_all$Specific_antigen <- factor(BCR_all$Specific_antigen, 
                                   levels = c("FALSE", "Not_tested", 
                                              "CASPR2", "LGI1"))
#This looks right. 
write.csv(BCR_all, "Data/BCR_database_versions/7_Subclones_included.csv", 
          row.names = FALSE)