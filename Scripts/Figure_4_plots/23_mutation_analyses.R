library(ggplot2)

#For the same clones, we also identify how many mutations they have, and how they relate to another. 
allBCR <- read.csv("Data/BCR_database_versions/8_new_mutations_added.csv")

#Now, we are going to identify the full VDJ sequence distance between all
#clonal members in all clones, both L and H. 
#How large is the largest clone?
clonalBCR <- allBCR[which(allBCR$Clonal),]
BCRList <- split(clonalBCR, f = paste0(clonalBCR$HamClone, "_", clonalBCR$LOCUS))
#Now, what is the largest number of clonal members? We can of course remove one cell
#as we are going to compare intraclonal differences, and we do not need to make the
#calculations for the cell against itself. 
maxCloneSize <- max(unlist(lapply(BCRList, nrow)))-1

vdjDist_intraclonalIMGT <- do.call("rbind", lapply(clonalBCR$SEQUENCE_ID, function(x){
    print(x)
    locCell <- clonalBCR[which(clonalBCR$SEQUENCE_ID == x),]
    locCloneName <- paste0(locCell$HamClone, "_", locCell$LOCUS)
    locClone <- BCRList[[which(names(BCRList) == locCloneName)]][,c("CELL", "SEQUENCE_IMGT", "JUNCTION_LENGTH")]
    locClone <- locClone[-which(locClone$CELL == locCell$CELL),]
    locClone$MUT_COMP <- locCell$SEQUENCE_IMGT
    mutCols <- colnames(locCell)[grep("mu_count_", colnames(locCell))]
    totalCols <- c(colnames(locClone), mutCols)
    locRes <- tryCatch(observedMutations(locClone, 
                                         sequenceColumn = "SEQUENCE_IMGT",
                                         germlineColumn = "MUT_COMP",
                                         regionDefinition = IMGT_VDJ_BY_REGIONS,
                                         juncLengthColumn = "JUNCTION_LENGTH",
                                         nproc = 1),
                       error=function(e){
                           errDf <- as.data.frame(matrix(NA, nrow(locClone),
                                                         length(totalCols),
                                                         dimnames = list(1:nrow(locClone),
                                                                         totalCols)))
                           errDf$CELL <- locClone$CELL
                           errDf
                       })
    locResMut <- locRes[,grep("mu_count_", colnames(locRes))]
    locResMutRow <- do.call("cbind", lapply(1:maxCloneSize, function(y){
        if(y <= nrow(locResMut)){
            locCell <- locRes$CELL[y]
            locDat <- locResMut[which(locRes$CELL == locCell),]
            locDat$name <- locCell
            colnames(locDat) <- paste0("Cell_", y, "_", colnames(locDat))
            locDat
        } else {
            as.data.frame(matrix(NA, 1, 15, dimnames = list(NA, paste0("Cell_", y, "_", c(colnames(locResMut), "name")))))
        }
    }))
}))

dir.create("Results/Figure_4_plots/Mutation_analysis")

#Now, we will divide the data into the heavy and light gene parts, and then generate
#graphs for them separately. 
h_or_l <- clonalBCR$LOCUS
h_or_l[-which(clonalBCR$LOCUS == "H")] <- "L"
vdjDistSplit <- split(vdjDist_intraclonalIMGT, f = h_or_l)

regionMutNameVec <- c("cdr1_r", "cdr1_s",
                      "cdr2_r", "cdr2_s",
                      "cdr3_r", "cdr3_s",
                      "fwr1_r", "fwr1_s",
                      "fwr2_r", "fwr2_s",
                      "fwr3_r", "fwr3_s",
                      "fwr4_r", "fwr4_s")

right_region_order <- c("fwr1", "cdr1", 
                        "fwr2", "cdr2",
                        "fwr3", "cdr3", 
                        "fwr4")
regionNameVec <- unique(substr(regionMutNameVec, 1, 4))

#We also need the average length of each of the segments. They are howqever all standardised in length
#so we could as well just take the first number. 
regionLengths <- sapply(colnames(clonalBCR)[grep("CDR._IMGT|FWR._IMGT", colnames(clonalBCR))], function(x){
    median(nchar(clonalBCR[,which(colnames(clonalBCR) == x)]))
})

names(regionLengths) <- c("fwr1", "fwr2", "fwr3", "fwr4", "cdr1", "cdr2", "cdr3")
regionLengthsNorm <- regionLengths/max(regionLengths)

lapply(names(vdjDistSplit), function(z){
    locVdj <- vdjDistSplit[[which(names(vdjDistSplit) == z)]]
    #Now, we sum up the differences per region
    
    s_and_r_mu_counts <- do.call("cbind", lapply(regionMutNameVec, function(x){
        unlist(locVdj[,grep(x, colnames(locVdj))])
    }))
    colnames(s_and_r_mu_counts) <- regionMutNameVec
    #And we remove all rows with NAs in them. they are always complete NA rows, so we will not lose any
    #data this way. 
    s_and_r_mu_complete <- s_and_r_mu_counts[complete.cases(s_and_r_mu_counts),]
    
    r_mu_counts <- as.data.frame(s_and_r_mu_complete[,grep("_r", colnames(s_and_r_mu_complete))])
    s_mu_counts <- as.data.frame(s_and_r_mu_complete[,grep("_s", colnames(s_and_r_mu_complete))])
    
    all_mu_counts <- as.data.frame(do.call("cbind", lapply(regionNameVec, function(x){
        locR <- r_mu_counts[,which(grepl(x, colnames(r_mu_counts)))]
        locS <- s_mu_counts[,which(grepl(x, colnames(s_mu_counts)))]
        locRes <- locR+locS
    })))
    colnames(all_mu_counts) <- regionNameVec
    
    all_mu_ggdat <- reshape2::melt(all_mu_counts)
    colnames(all_mu_ggdat) <- c("BCR_region", "basepair_distance")
    all_mu_ggdat$BCR_region <- factor(all_mu_ggdat$BCR_region, levels = right_region_order)
    
    ggplot(all_mu_ggdat, aes(x = BCR_region, y = basepair_distance)) +
        geom_violin() +theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/All_mutations_intraclonal_", z, "_chain.pdf"))
    
    #Now for the replacements
    colnames(r_mu_counts) <- substr(colnames(r_mu_counts), 1,4)
    r_mu_ggdat <- reshape2::melt(r_mu_counts)
    colnames(r_mu_ggdat) <- c("BCR_region", "basepair_distance")
    r_mu_ggdat$BCR_region <- factor(r_mu_ggdat$BCR_region, levels = right_region_order)
    
    ggplot(r_mu_ggdat, aes(x = BCR_region, y = basepair_distance)) +
        geom_violin() + theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/Replacement_mutations_intraclonal_", z, "_chain.pdf"))
    
    #And now the same, but corrected for the number of bases in each segment: 
    all_mu_ggdat_corr <- do.call("rbind", lapply(unique(all_mu_ggdat$BCR_region), function(x){
        locDat <- all_mu_ggdat[which(all_mu_ggdat$BCR_region == x),]
        locDat$corr_dist <- locDat$basepair_distance/regionLengthsNorm[which(names(regionLengthsNorm) == x)]
        locDat
    }))
    ggplot(all_mu_ggdat_corr, aes(x = BCR_region, y = corr_dist)) +
        geom_violin() +theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/All_mutations_intraclonal_", z, "_chain_region_length_corrected.pdf"))
    
    r_mu_ggdat_corr <- do.call("rbind", lapply(unique(r_mu_ggdat$BCR_region), function(x){
        locDat <- all_mu_ggdat[which(all_mu_ggdat$BCR_region == x),]
        locDat$corr_dist <- locDat$basepair_distance/regionLengthsNorm[which(names(regionLengthsNorm) == x)]
        locDat
    }))
    ggplot(r_mu_ggdat_corr, aes(x = BCR_region, y = corr_dist)) +
        geom_violin() + theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/Replacement_mutations_intraclonal_", z, "_chain_region_length_corrected.pdf"))
    
})

#Now, we make similar graphs for the total number of mutations from germline. 
BCR_mut_split <- split(clonalBCR[,grep("mu_", colnames(clonalBCR))], f = h_or_l)

lapply(names(BCR_mut_split), function(z){
    locVdj <- BCR_mut_split[[which(names(BCR_mut_split) == z)]]
    
    r_mu_counts <- locVdj[,grep("_r", colnames(locVdj))]
    s_mu_counts <- locVdj[,grep("_s", colnames(locVdj))]
    
    all_mu_counts <- as.data.frame(do.call("cbind", lapply(regionNameVec, function(x){
        locR <- r_mu_counts[,which(grepl(x, colnames(r_mu_counts)))]
        locS <- s_mu_counts[,which(grepl(x, colnames(s_mu_counts)))]
        locRes <- locR+locS
    })))
    colnames(all_mu_counts) <- regionNameVec
    
    all_mu_ggdat <- reshape2::melt(all_mu_counts)
    colnames(all_mu_ggdat) <- c("BCR_region", "basepair_distance")
    all_mu_ggdat$BCR_region <- factor(all_mu_ggdat$BCR_region, levels = right_region_order)
    
    ggplot(all_mu_ggdat, aes(x = BCR_region, y = basepair_distance)) +
        geom_violin() +theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/All_mutations_to_germline_", z, "_chain.pdf"))
    
    #Now for the replacements
    colnames(r_mu_counts) <- substr(colnames(r_mu_counts), 10,13)
    r_mu_ggdat <- reshape2::melt(r_mu_counts)
    colnames(r_mu_ggdat) <- c("BCR_region", "basepair_distance")
    r_mu_ggdat$BCR_region <- factor(r_mu_ggdat$BCR_region, levels = right_region_order)
    
    ggplot(r_mu_ggdat, aes(x = BCR_region, y = basepair_distance)) +
        geom_violin() + theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/Replacement_mutations_to_germline_", z, "_chain.pdf"))
    
    #And now the same, but corrected for the number of bases in each segment: 
    all_mu_ggdat_corr <- do.call("rbind", lapply(unique(all_mu_ggdat$BCR_region), function(x){
        locDat <- all_mu_ggdat[which(all_mu_ggdat$BCR_region == x),]
        locDat$corr_dist <- locDat$basepair_distance/regionLengthsNorm[which(names(regionLengthsNorm) == x)]
        locDat
    }))
    ggplot(all_mu_ggdat_corr, aes(x = BCR_region, y = corr_dist)) +
        geom_violin() +theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/All_mutations_to_germline_", z, "_chain_region_length_corrected.pdf"))
    
    r_mu_ggdat_corr <- do.call("rbind", lapply(unique(r_mu_ggdat$BCR_region), function(x){
        locDat <- all_mu_ggdat[which(all_mu_ggdat$BCR_region == x),]
        locDat$corr_dist <- locDat$basepair_distance/regionLengthsNorm[which(names(regionLengthsNorm) == x)]
        locDat
    }))
    ggplot(r_mu_ggdat_corr, aes(x = BCR_region, y = corr_dist)) +
        geom_violin() + theme_bw()
    ggsave(paste0("Results/Figure_4_plots/Mutation_analysis/Replacement_mutations_to_germline_", z, "_chain_region_length_corrected.pdf"))
    
})


#Now, as it all seemed to make quite a lot of sense, we are now going on to summing up. 


complete_VDJ_df <- cbind(clonalBCR[,c("CELL", "LOCUS", "HamClone", 
                                      colnames(clonalBCR)[grep("mu_",colnames(clonalBCR))])],
                         vdjDist_intraclonalIMGT)

germMuCols <- which(grepl("mu_", colnames(complete_VDJ_df )) & grepl("Cell_", colnames(complete_VDJ_df)) == FALSE)

colnames(complete_VDJ_df)[germMuCols] <- paste0("Germline_", colnames(complete_VDJ_df)[germMuCols])

VDJ_df_only_mut <- complete_VDJ_df[,grep("mu_", colnames(complete_VDJ_df))]
    
groups <- unique(gsub("|_mu_count_.+", "", colnames(VDJ_df_only_mut)))

groupedMutStats <- do.call("cbind", lapply(groups, function(x){
    locDat <- VDJ_df_only_mut[,grep(x, colnames(VDJ_df_only_mut))]
    locAllMut <- apply(locDat, 1, sum)
    locRepMut <- apply(locDat[,grep("_r", colnames(locDat))], 1, sum)
    locDf <- data.frame(locAllMut,locRepMut)
    colnames(locDf) <- paste0(x, c("_all_mut", "_replacement_mut"))
    locDf
}))

#Now, we sum this up on a per-cell basis, 

mutStatsPerCell <- do.call("rbind", lapply(unique(clonalBCR$CELL), function(x){
    locClone <- clonalBCR$HamClone[which(clonalBCR$CELL == x)][1]
    locSubClone <- clonalBCR$SubClone[which(clonalBCR$CELL == x)][1]
    data.frame("CELL" = x, "HamClone" = locClone, 
                        "SubClone" = locSubClone, 
                        t(colSums(groupedMutStats[which(clonalBCR$CELL == x),])))
}))

#And now, we calculate the furthest intraclonal distance for each chain. 
mutStatsPerCell$maxAllMutDist <- apply(mutStatsPerCell[,-grep("CELL|Clone|Germline",colnames(mutStatsPerCell))],1,max, na.rm = TRUE)
mutStatsPerCell$maxRepMutDist <- apply(mutStatsPerCell[,-grep("CELL|Clone|Germline|_all",colnames(mutStatsPerCell))],1,max, na.rm = TRUE)


write.csv(mutStatsPerCell, "Results/Figure_4_plots/Mutation_analysis/mutStatsPerCell.csv", row.names = FALSE)

#Here, we now do the same, but only with the CDR3 mutations
CDR3MutStats <- VDJ_df_only_mut[,grep("cdr3", colnames(VDJ_df_only_mut))]
CDR3MutStatsPerCell <- do.call("rbind", lapply(unique(clonalBCR$CELL), function(x){
    locClone <- clonalBCR$HamClone[which(clonalBCR$CELL == x)][1]
    locSubClone <- clonalBCR$SubClone[which(clonalBCR$CELL == x)][1]
    data.frame("CELL" = x, "HamClone" = locClone, 
               "SubClone" = locSubClone, 
               t(colSums(CDR3MutStats[which(clonalBCR$CELL == x),])))
}))
write.csv(CDR3MutStatsPerCell, "Results/Figure_4_plots/Mutation_analysis/CDR3mutStatsPerCell.csv", row.names = FALSE)

#And now, we sum up all until the CDR3. 
nonCDR3MutStats <- VDJ_df_only_mut[,-grep("cdr3|fwr4", colnames(VDJ_df_only_mut))]
nonCDR3GroupedMutStats <- do.call("cbind", lapply(groups, function(x){
    locDat <- nonCDR3MutStats[,grep(x, colnames(nonCDR3MutStats))]
    locAllMut <- apply(locDat, 1, sum)
    locRepMut <- apply(locDat[,grep("_r", colnames(locDat))], 1, sum)
    locDf <- data.frame(locAllMut,locRepMut)
    colnames(locDf) <- paste0(x, c("_all_mut", "_replacement_mut"))
    locDf
}))

nonCDR3MutStatsPerCell <- do.call("rbind", lapply(unique(clonalBCR$CELL), function(x){
    locClone <- clonalBCR$HamClone[which(clonalBCR$CELL == x)][1]
    locSubClone <- clonalBCR$SubClone[which(clonalBCR$CELL == x)][1]
    data.frame("CELL" = x, "HamClone" = locClone, 
               "SubClone" = locSubClone, 
               t(colSums(nonCDR3GroupedMutStats[which(clonalBCR$CELL == x),])))
}))
write.csv(nonCDR3MutStatsPerCell, "Results/Figure_4_plots/Mutation_analysis/Non-CDR3mutStatsPerCell.csv", row.names = FALSE)


