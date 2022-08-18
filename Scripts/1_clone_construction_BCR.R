library(ggplot2)
library(shazam)
library(alakazam)
#First, we import all,and combine them into one big frame
bracerDirs <- paste0(list.files("Data/Bracer_output", full.names = TRUE, pattern = "filtered"),
                     "/IMGT_gapped_db.tab")

bracerFileList <- lapply(bracerDirs, read.delim)
bracerFile <- do.call("rbind",bracerFileList)

#One thing we want to do early is to add the TPM information to all the sequences
#as it has bearings at multiple stages. 

TPM_data_dirs <- paste0(list.files("Data/Bracer_output", 
                                   pattern = "filtered_BCR_summary", 
                                   full.names = TRUE),
                        "/changeodb.tab")

TPMFileList <- lapply(TPM_data_dirs, read.delim)
TPMFile <- do.call("rbind",TPMFileList)

#Now, we check if the two files have identical rows. They did not at first: 
#two very incomplete sequences for JR1166 - take 1 (between position 300 and 310
#as well as 500 and 510 in the file) had been cleared from the  gapped file. 
#They are now removed also from the TPM file. In addition, one cell had moved 
#two cells down in donor JR1227 take 2, plate 4, positions 1030:1035. I moved it
#to the matching place. After this, I checked that the sequence IDs were right. 
#All were, apart from for JR1227_2 plate 4, where the BCR_data sequence IDs for
#some reason had been converted from 4 to 43. I changed that in the file, and
#suddently it all made sense. 
identical(bracerFile$SEQUENCE_ID, TPMFile$SEQUENCE_ID)
identical(bracerFile$CELL, TPMFile$CELL)

#Now, we have two matching files, and can thus add the TPM information and the non-gapped VDJ
bracerFile$TPM <- TPMFile$TPM
bracerFile$SEQUENCE_VDJ_UNGAPPED <- TPMFile$SEQUENCE_VDJ

#In a piece of important work, presented in the file 
#All_bracer_results/Clone_veracity_check, a number of cells as contaminated with
#BCRs from surrounding cells are defined. Two extra (1284_1_5_N01 and 1284_1_5_H02)
#have been found in iterations of clustering and close consideration in 4a1 (close wells having 
#being "clonal" as the cell in question expresses very, tiny levels of a highly expressed
#transcript from the other well)
BCR_contaminations <- c("1284_1_5_N01", "1284_1_5_N02", 
                        "1284_1_5_N05", "1284_1_5_H02",
                        "1284_1_1_F07", "1284_1_1_I07",
                        "1284_1_5_O04", "1284_1_5_O05", 
                        "1284_1_5_O06", "1284_1_5_H03",
                        "1284_1_5_P02", "1284_1_5_P03", 
                        "1284_1_5_P04", "1284_1_5_P05",
                        "1284_1_5_P06")

#These are of course removed
bracerFile <- bracerFile[-which(bracerFile$CELL %in% BCR_contaminations),]

#Before going in to the clustering steps, we will remove incomplete BCRs and likely doublet cells

#We remove cells with incomplete BCRs or non-functional, as they will not be possible to
#work with downstream anyway with the current primary focus on the BCR and TCR and
#only secondary interest in the transcriptomes. 
uniqueCells <- unique(bracerFile$CELL)
locFunc <- paste0(bracerFile$LOCUS, bracerFile$FUNCTIONAL)
fullCells <- unlist(lapply(uniqueCells, function(x){
    locChains <- locFunc[which(bracerFile$CELL == x)]
    if(length(locChains) > 1){
        if(any(locChains == "HTRUE") && any(locChains %in% c("KTRUE", "LTRUE"))){
            TRUE
        } else {
            FALSE
        }
    } else {
        FALSE
    }
}))

table(fullCells)
#This makes us lose 94 incomplete cells
uniqueCells <- uniqueCells[which(fullCells)]

write.csv(bracerFile[which(bracerFile$CELL %in% uniqueCells),], "Data/IMGT_gapped_db_TPM_added.csv", 
          row.names =FALSE)

#Here cells with more than one functional H and K or L that have an edit distance
#above 10 are removed, as such cells are highly likely to be doublets, and thus
#should be excluded from all downstream analyses
funcChains <- bracerFile[which(bracerFile$FUNCTIONAL),]
doubletCells <- unlist(lapply(uniqueCells, function(x){
    locCell <- funcChains[which(funcChains$CELL == x),]
    doubleChains <- unlist(lapply(c("H", "K", "L"), function(y){
        locChain <- locCell[which(locCell$LOCUS == y),]
        if(nrow(locChain) > 1){
            distMat <- adist(locChain$JUNCTION)
            if(max(distMat) > 5){
                TRUE
            } else {
                FALSE
            }
        } else {
            FALSE
        }
    }))
    if(length(which(doubleChains)) > 1){
        TRUE
    } else {
        FALSE
    }
}))

table(doubletCells)

#This removes 4 cells that have one or more unrelated extra functional chain, for at least 
#two of H, K and L. 
uniqueCells <- uniqueCells[-which(doubletCells)]

#We of course also remove all empty wells, denoted by at least one X in the file
#name
uniqueCells <- uniqueCells[-grep("X", uniqueCells)]

#Making us lose another 9 cells, but after this, we are down to 568 nice cells

bracerFile <- bracerFile[which(bracerFile$CELL %in% uniqueCells),]

#In a previous test, we established that with the current set of cells, that are likely to be all
#true, there are only four that have double heavy chains and that are clonal. Among 
#these four, selecting the dominant chain will always lead to inclusion into the largest
#clone, so this is done here, for all chains, as it makes the downstream analyses
#so much cleaner. After this, all cells will have three chains, but in most
#instances, one of K and L will not be present. From analyses in 4a2c, it is clear
#that a significant umber of cells do express both K and L, so they are both kept, 
#instead of just keeping a "light chain". 
bracerFile3Chains <- do.call("rbind", lapply(c("H", "K", "L"), function(i){
    bracerFile_focus <- bracerFile[which(bracerFile$LOCUS == i),]
    bracerFile_focus_condensed <- 
        do.call("rbind", lapply(uniqueCells,  function(x){
            locRows <- which(bracerFile_focus$CELL == x)
            if(length(locRows) > 0){
                locDat <- bracerFile_focus[locRows,]
                #In the unlikely scenario that there are 
                #two functional chains, or only 
                #two non-functional chains for that 
                #matter, we here select the chain
                #with the highest expression
                if(nrow(locDat) > 1){
                    locDat[which.max(locDat$TPM),]
                } else {
                    locDat
                }
            } else {
                NA
            }
            
        }))
}))


#Now, we are down to 1515 chains from the remaining, complete cells. 
#And this is saved
write.csv(bracerFile3Chains, "Data/IMGT_gapped_db_TPM_added_3_chains.csv", 
          row.names =FALSE)

#Now, we use the standard Hamming-based definition of clones from Chanteo/shazam.
#We leave the one provided by BraCeR, as it is hard to trace what has actually happened there. 
#Step one here is to extract the relevant sequences, 
bracerFileH <- bracerFile3Chains[which(bracerFile3Chains$LOCUS == "H"),]

#Now, we cluster the data with a conservative, Hamming-based method, which 
#has been compared to other methods and been identified as giving idetical results.
dist_ham <- distToNearest(bracerFileH, sequenceColumn="JUNCTION", 
                          vCallColumn="V_CALL", jCallColumn="J_CALL",
                          model="ham", normalize="len", nproc=6)

#Among the threshold selection methods, it is the simpler density method that
#works best for this data, considering both dist methods. 
output <- findThreshold(dist_ham$dist_nearest, method="density")
plot(output, binwidth=0.02, title="Density Method")

#Now, we go on to ChangeO and we get the clustering results, but we start by
changeoClusteringFile <- bracerFileH[,c("CELL", "SEQUENCE_ID", "V_CALL", "J_CALL",
                                        "JUNCTION", "TPM")]
write.table(changeoClusteringFile, "Data/changeoClusteringFile.tsv", 
            sep = "\t",row.names = FALSE)

system(paste0("DefineClones.py -d Data/changeoClusteringFile.tsv -o Data/changeoClusteringFile_Ham_result.tsv --act set --model ham --norm len --format changeo --sf JUNCTION --vf V_CALL --jf J_CALL --dist ", output@threshold))

#Now, we re-import these to see how much they overlap
hamClones <- read.table("Data/changeoClusteringFile_Ham_result.tsv", header = TRUE)

uniqueCells <- unique(hamClones$CELL)

bracerFileHam <- do.call("rbind", lapply(uniqueCells, function(x){
    locCell <- bracerFile3Chains[which(bracerFile3Chains$CELL == x),]
    locCell$HamClone <- hamClones$CLONE[which(hamClones$CELL == x)]
    locCell
}))

#Here, we make one final addition before saving: we add the heavy junction 
#edit distance within each clone.

bracerFileHam$intraClonalDistance <- NA
cloneNames <- unique(bracerFileHam$HamClone)
for(i in cloneNames){
    locRows <- which(bracerFileHam$HamClone == i & bracerFileHam$LOCUS == "H")
    uniqueLocNames <- unique(bracerFileHam$CELL[locRows])
    if(length(uniqueLocNames) > 1){
        locNames <- bracerFileHam$CELL[locRows]
        locCloneDist <- adist(bracerFileHam[locRows,"JUNCTION"])
        #To ensure that we do not get zero everywhere, we remove the diagonal
        diag(locCloneDist) <- 1000
        colnames(locCloneDist) <- locNames
        for(j in uniqueLocNames){
            if(length(which(locNames == j)) > 1){
                innerCloneDist <- locCloneDist[,which(locNames == j)]
                minVal <- min(innerCloneDist[,which.min(apply(innerCloneDist,2,min))])
            } else {
                minVal <- min(locCloneDist[,j])
            }
            bracerFileHam$intraClonalDistance[which(bracerFileHam$CELL == j)] <- 
                minVal
        }
    }
}

table(bracerFileHam$intraClonalDistance)
#  0   1   2   3   4   5   6 
#262  20   6  17   6   5   4

#We also un-name all clones that consist of one cell only
bracerFileHam$HamClone[which(is.na(bracerFileHam$intraClonalDistance))] <- NA


#And with that, we have a very clean dataset, with 538 cells, where the largest
#intra-clonal junction edit distance is 6

write.csv(bracerFileHam, "Data/IMGT_gapped_db_complete_post_clonality.csv", 
          row.names = FALSE)

