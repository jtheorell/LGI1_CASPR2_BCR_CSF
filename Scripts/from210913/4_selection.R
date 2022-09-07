#Here, we take the file from 2 and first exclude all non-IgGs as the project
#will focus on them. Then we for each clone select one BCR with the shortest
#intraclonal distance. After this, we randomly select 25% of the non-clonal
#sequences from the three donors separately. 
BCR_data_B_cells <- read.csv("Complete_BCR_db.csv")

#Here we exclude the non-acute samples
BCR_data_B_cells_acute <- BCR_data_B_cells[which(substr(BCR_data_B_cells$CELL, 1, 6) %in% 
                                                     c("1166_1", "1227_1", "1284_1")),]

#As a few chains are dubiously short, we will exclude them here. 

hist(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "H")]))
quantile(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "H")]), 0.1)
#341
hist(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "K")]))
quantile(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "K")]), 0.1)
#292
hist(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "L")]))
quantile(nchar(BCR_data_B_cells_acute$SEQUENCE_VDJ_UNGAPPED[which(BCR_data_B_cells_acute$LOCUS == "L")]), 0.1)
#318

#This shows that all three chain loci have significant lower tails. We will
#exclude cells that has one chain in the lowest 10% from any of the three loci.
short_chains <- unique(unlist(lapply(split(BCR_data_B_cells_acute, 
                                           f = BCR_data_B_cells_acute$LOCUS), function(x){
                                               lowQuant <- 
                                                   quantile(nchar(x$SEQUENCE_VDJ_UNGAPPED), 0.1)
                                               x$CELL[which(nchar(x$SEQUENCE_VDJ_UNGAPPED) < lowQuant)]
                                           })))

length(short_chains)/length(unique(BCR_data_B_cells_acute$CELL))
#0.1832359
#This means that 19% of the cells have at least one chain that is short. 

#Here, we focus the analysis on the heavy chains. 
BCR_data_H <- BCR_data_B_cells_acute[which(BCR_data_B_cells_acute$LOCUS == "H"),]

#And now, on the cells with all chains full-length.
BCR_data_H_long <- BCR_data_H[-which(BCR_data_H$CELL %in% short_chains),]

#Now, we exclude all non-IgG. 
BCR_data_H_IgG <- BCR_data_H_long[grep("IGHG", BCR_data_H_long$ISOTYPE),]

table(BCR_data_H_IgG$ISOTYPE, BCR_data_H_IgG$Clonal)
#FALSE TRUE
#IGHG1    67   16
#IGHG2    22    9
#IGHG3     6    0
#IGHG4    84   69

#We save the IgG3s anyway
BCR_data_H_clonal <- BCR_data_H_IgG[which(BCR_data_H_IgG$Clonal),]

BCR_H_clonal_selected <- do.call("rbind",lapply(split(BCR_data_H_clonal, 
                                      f = BCR_data_H_clonal$HamClone),
                                function(x){
                                    x[which.min(x$intraClonalDistance)[1],]
                                }))

#And now over to the others: 
BCR_data_H_nonclonal <- BCR_data_H_IgG[-which(BCR_data_H_IgG$Clonal),]

BCR_H_nonclonal_selected <- do.call("rbind",lapply(split(BCR_data_H_nonclonal , 
                                                      f = substr(BCR_data_H_nonclonal$CELL, 1, 4)),
                                                function(x){
                                                    x[sample(1:nrow(x), 
                                                             round(nrow(x)/4)),]
                                                }))

#Now, we go into the final phase. Here, it feels appropriate to exclude the
#non-acute samples
BCR_H_selected <- rbind(BCR_H_clonal_selected, BCR_H_nonclonal_selected)

#We also exclude the sequences where we already have the answer
BCR_H_selected <- BCR_H_selected[which(is.na(BCR_H_selected$LGI1Score)),]

BCR_selected <- BCR_data_B_cells[which(BCR_data_B_cells$CELL %in% BCR_H_selected$CELL),]

write.csv(BCR_selected, "BCR_selected_db.csv", row.names = FALSE)

