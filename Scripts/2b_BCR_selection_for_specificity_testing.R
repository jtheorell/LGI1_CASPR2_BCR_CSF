#Here, we take the file from 2 and first exclude all non-IgGs as the project
#will focus on them. Then we for each clone select one BCR with the shortest
#intraclonal distance. After this, we randomly select 25% of the non-clonal
#sequences from the three donors separately. 
BCR_data <- read.csv("Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv")

##############
#FOR FIGURE 1
#############
table(BCR_data$Sample[which(BCR_data$LOCUS == "H")])

#1166_1 1227_1 1284_1 
#211     63     91

#As a few chains are dubiously short, we will exclude the shortest 10% here. 

hist(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "H")]))
quantile(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "H")]), 0.1)
#339.4
hist(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "K")]))
quantile(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "K")]), 0.1)
#275.2
hist(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "L")]))
quantile(nchar(BCR_data$SEQUENCE_VDJ_UNGAPPED[which(BCR_data$LOCUS == "L")]), 0.1)
#304.6

#This shows that all three chain loci have significant lower tails. We will
#exclude cells that has the heavy and the light chain, or in the case of K/L cells
#both of them, in this shortest 10% fraction. It turns out that of the 16 cells
#that have both a K and an L chain, 7 have one that is among the shortest 10%, so
#that probably shows the most likely candidate in all those cases too. Those
#cells will not be excluded here. 

short_chains <- lapply(split(BCR_data, f = BCR_data$LOCUS), function(x){
                                               lowQuant <- 
                                                   quantile(nchar(x$SEQUENCE_VDJ_UNGAPPED), 0.1)
                                               x$CELL[which(nchar(x$SEQUENCE_VDJ_UNGAPPED) < lowQuant)]
                                           })
#Are there any cells that have both a too short K and a too shot L chain?
length(which(short_chains$K %in% short_chains$L))
#2
#So we will need to do this separately for the double cells. 
doubleShort <- short_chains[which(short_chains$K %in% short_chains$L)]

shortAll <- unique(unlist(short_chains))
shortNonDouble <- 
    shortAll[-which(shortAll %in% BCR_data$CELL[which(BCR_data$light_type == "double")])]

realShort <- c(doubleShort, shortNonDouble)

length(realShort)/length(unique(BCR_data$CELL))
#0.1780822
#This means that 18% of the cells have at least one chain that is short. 

#Here, we focus the analysis on the heavy chains. 
BCR_data_H <- BCR_data[which(BCR_data$LOCUS == "H"),]

#And now, on the cells with all chains full-length.
BCR_data_H_long <- BCR_data_H[-which(BCR_data_H$CELL %in% realShort),]

##############
#FOR FIGURE 1
#############
length(realShort)
#65

table(BCR_data_H_long$ISOTYPE, BCR_data_H_long$Clonal)
#      FALSE TRUE
#IGHG1    75   17
#IGHG2    24   10
#IGHG3     6    0
#IGHG4    95   75

#We save the IgG3s anyway
BCR_data_H_clonal <- BCR_data_H_long[which(BCR_data_H_long$Clonal),]

BCR_H_clonal_selected <- do.call("rbind",lapply(split(BCR_data_H_clonal, 
                                      f = BCR_data_H_clonal$HamClone),
                                function(x){
                                    x[which.min(x$intraClonalDistance)[1],]
                                }))

#And now over to the others: 
BCR_data_H_nonclonal <- BCR_data_H_long[-which(BCR_data_H_long$Clonal),]

BCR_H_nonclonal_selected <- do.call("rbind",lapply(split(BCR_data_H_nonclonal , 
                                                      f = BCR_data_H_nonclonal$Sample),
                                                function(x){
                                                    x[sample(1:nrow(x), 
                                                             round(nrow(x)/4)),]
                                                }))

BCR_H_selected <- rbind(BCR_H_clonal_selected, BCR_H_nonclonal_selected)

BCR_selected <- BCR_data[which(BCR_data$CELL %in% BCR_H_selected$CELL),]

write.csv(BCR_selected, "Data/BCR_database_versions/5b_BCR_selected_db.csv", row.names = FALSE)

#As the exclusion criteria thresholds have previously been somehat harsher (especially
#exclusion of cells not fitting the surface B cell phenotype, a process that previously
#was unduly crude), the current, above selection would have resulted in approximately
#8 more cells being included. Furthermore, as the original selection was done without
#setting the pseudo-random seed, it would have been impossible to reproduce the exact 
#same set even if the exclusion criteria had been static. Below is therefore a file
#with a list of the cells in the original selection. The numbers for this selection
#is what is reported in the paper. As can be seen above, the intention was to select
#25% of the sequences, but as it got closer to 20% in the end, we report this number.
original_selection <- read.csv("Data/BCR_auxiliaries/Originally_selected_BCRs.csv")[,1]

originalClonality <- BCR_data_H_long$Clonal[which(BCR_data_H_long$CELL 
                                                  %in% original_selection)]

originalSample <- BCR_data_H_long$Sample[which(BCR_data_H_long$CELL 
                                               %in% original_selection)]

table(originalSample, originalClonality)

originalClonality
#originalSample FALSE TRUE
#1166_1    27   26
#1227_1     7    4
#1284_1    10   11
