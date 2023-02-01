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

#And here, we save the fully complete dataset. 
dir.create("Data/BCR_database_versions")

##############
#Removals
#We start by removing all non-cells, that have slipped through as well as the positive
#controls in A1. 
preAnalyticRemoval <- grep("X|A01", bracerFile$CELL)
bracerFile <- bracerFile[-preAnalyticRemoval,]

bracerFile$Sample <- substr(bracerFile$CELL, 1,6)

##############
#FOR TABLE 1
length(unique(bracerFile$CELL))
#This gives us 603 in total
lapply(split(bracerFile$CELL, bracerFile$Sample), function(x) length(unique(x)))
#$`1166_1`
#[1] 234
#
#$`1227_1`
#[1] 104
#
#$`1284_1`
#[1] 265
##############

#Before going in to the clustering steps, we will remove incomplete BCRs and likely doublet cells

#First, we remove all non-functional chains. 
table(bracerFile$FUNCTIONAL)
#FALSE  TRUE 
#245  1256

bracerFile <- bracerFile[which(bracerFile$FUNCTIONAL),]
dim(bracerFile)
#1256   59
length(unique(bracerFile$CELL))
#603 
#So no cells have only non-functional chains. 

#We remove cells with incomplete BCRs here. 
uniqueCells <- unique(bracerFile$CELL)
fullCells <- unlist(lapply(uniqueCells, function(x){
    locChains <- bracerFile$LOCUS[which(bracerFile$CELL == x)]
    if(length(locChains) > 1){
        if(any(locChains == "H") && any(locChains %in% c("K", "L"))){
            TRUE
        } else {
            FALSE
        }
    } else {
        FALSE
    }
}))

##############
#FOR TABLE 1
table(fullCells)
#fullCells
#FALSE  TRUE 
#67     536
##############

#This makes us lose 67 incomplete cells
uniqueCells <- uniqueCells[which(fullCells)]

bracerFile <- bracerFile[which(bracerFile$CELL %in% uniqueCells),]
dim(bracerFile)
#1184   59

length(unique(bracerFile$CELL))
#536

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
dim(bracerFile)
#1184   59
bracerFile <- bracerFile[-which(bracerFile$CELL %in% BCR_contaminations),]
dim(bracerFile)
#1131   59

length(unique(bracerFile$CELL))
#521
#So this made us lose another 15 cells. 

#Here cells with more than one functional H and K or L that have an edit distance
#above 10 are removed, as such cells are highly likely to be doublets, and thus
#should be excluded from all downstream analyses

doubletCells <- unlist(lapply(unique(bracerFile$CELL), function(x){
    locCell <- bracerFile[which(bracerFile$CELL == x),]
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
#FALSE  TRUE 
#520     1 

#This removes 1 cell that has one or more unrelated extra functional chain, for at least 
#two of H, K and L. 
singletCells <- unique(bracerFile$CELL)[-which(doubletCells)]
bracerFile <- bracerFile[which(bracerFile$CELL %in% singletCells),]

#And so, we are ready to go on to the next step. 
write.csv(bracerFile, "Data/BCR_database_versions/1_IMGT_gapped_db_TPM_added.csv", 
          row.names =FALSE)

#In a previous test, we established that with the current set of cells, that are likely to be all
#true, there are only four that have double heavy chains and that are clonal. Among 
#these four, selecting the dominant chain will always lead to inclusion into the largest
#clone, so this is done here, for all chains, as it makes the downstream analyses
#so much cleaner.
bracerFile1H <- do.call("rbind", lapply(c("H", "K", "L"), function(i){
    bracerFile_focus <- bracerFile[which(bracerFile$LOCUS == i),]
    bracerFile_focus_condensed <- 
        do.call("rbind", lapply(unique(bracerFile$CELL),  function(x){
            locRows <- which(bracerFile_focus$CELL == x)
            if(length(locRows) > 0){
                locDat <- bracerFile_focus[locRows,]
                #Here the chain with the highest expression is selected.
                if(nrow(locDat) > 1){
                    locDat[which.max(locDat$TPM),]
                } else {
                    locDat
                }
            }
        }))
}))

dim(bracerFile1H)
#1074   59
#Furthermore, in the cases where there are two light chains (K+L), it is always the more
#highly expressed one that is specific, after a number of tests. For this reason, 
#we can exclude all lowly expressed chains in cases where two are present. 
bracerFile2Chain <- do.call("rbind", lapply(unique(bracerFile1H$CELL), function(x){
    locDat <- bracerFile1H[which(bracerFile1H$CELL == x),]
    if(nrow(locDat) > 2){
        locLightDat <- locDat[which(locDat$LOCUS != "H"),]
        locDat <- rbind(locDat[which(locDat$LOCUS == "H"),],
                        locLightDat[which.max(locLightDat$TPM[]),])
    }
    locDat
}))
dim(bracerFile2Chain)
#1040   59
#Now, we are down to 1040 chains from the remaining, complete cells, i.e. 2 chains per cell. 
#And this is saved
write.csv(bracerFile2Chain, "Data/BCR_database_versions/3_IMGT_gapped_db_TPM_added_1H_1L.csv", 
          row.names =FALSE)

#Now, we use the standard Hamming-based definition of clones from Chanteo/shazam.
#We leave the one provided by BraCeR, as it is hard to trace what has actually happened there. 
#Step one here is to extract the relevant sequences, 
bracerFileH <- bracerFile2Chain[which(bracerFile2Chain$LOCUS == "H"),]

#Now, we cluster the data with a conservative, Hamming-based method, which 
#has been compared to other methods and been identified as giving idetical results.
dist_ham <- distToNearest(bracerFileH, sequenceColumn="JUNCTION", 
                          vCallColumn="V_CALL", jCallColumn="J_CALL",
                          model="ham", normalize="len", nproc=8)

#Among the threshold selection methods, it is the simpler density method that
#works best for this data, considering both dist methods. 
output <- findThreshold(dist_ham$dist_nearest, method="density")
plot(output, binwidth=0.02, title="Density Method")

#Now, we go on to ChangeO and we get the clustering results, but we start by
changeoClusteringFile <- bracerFileH[,c("CELL", "SEQUENCE_ID", "V_CALL", "J_CALL",
                                        "JUNCTION", "TPM")]
write.table(changeoClusteringFile, "Data/BCR_database_versions/changeoClusteringFile.tsv", 
            sep = "\t",row.names = FALSE)

system(paste0("DefineClones.py -d Data/BCR_database_versions/changeoClusteringFile.tsv -o Data/BCR_database_versions/changeoClusteringFile_Ham_result.tsv --act set --model ham --norm len --format changeo --sf JUNCTION --vf V_CALL --jf J_CALL --dist ", output@threshold))

#Now, we re-import these. 
hamClones <- read.table("Data/BCR_database_versions/changeoClusteringFile_Ham_result.tsv", header = TRUE)

#At this stage, we actually exchange the clone numbers for a standardised set used
#in previous rounds of analysis. They are in all respects equal
#to the ones generated above, apart from the actual numbers. 
originalCloneNames <- read.csv("Data/BCR_auxiliaries/Original_clone_names.csv")

bracerFileHam <- bracerFile2Chain
bracerFileHam$HamClone <- NA
for(i in bracerFileH$CELL){
    bracerFileHam$HamClone[which(bracerFileHam$CELL == i)] <- 
        originalCloneNames$HamClone[which(originalCloneNames$CELL == i)]
}
#There are warnings here, as we are using a single row as input for many rows,
#but that is correct. 

#The original code used to introduce the hamClone names was this: 
#uniqueCells <- unique(hamClones$CELL)
#
#bracerFileHam <- do.call("rbind", lapply(uniqueCells, function(x){
#    locCell <- bracerFile1H[which(bracerFile1H$CELL == x),]
#    locCell$HamClone <- hamClones$CLONE[which(hamClones$CELL == x)]
#    locCell
#}))

#Here we add the heavy and light junction 
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
#0   1   2   3   4   5   6 
#240  20   6  12   6   4   4

bracerFileHam$intraClonalDistanceL <- NA
cloneNames <- unique(bracerFileHam$HamClone)
for(i in cloneNames){
    locRows <- which(bracerFileHam$HamClone == i & bracerFileHam$LOCUS != "H")
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
            bracerFileHam$intraClonalDistanceL[which(bracerFileHam$CELL == j)] <- 
                minVal
        }
    }
}

table(bracerFileHam$intraClonalDistanceL)
#0   1   2   3   4  12 
#230  40   8   6   4   4 

#We also un-name all clones that consist of one cell only
bracerFileHam$HamClone[which(is.na(bracerFileHam$intraClonalDistance))] <- NA

#And with that, we have a very clean dataset, with 520 cells, where the largest
#intra-clonal junction edit distance for the heavy chain is 6

#Here, we add a column for the light type, which will be useful downstream
bracerFileHam$light_type <- NA
for(i in unique(bracerFileHam$CELL)){
    locRows <- which(bracerFileHam$CELL == i)
    locLoc <- bracerFileHam$LOCUS[locRows]
    lightLoc <- locLoc[which(locLoc != "H")]
    if(length(lightLoc) == 2){
        bracerFileHam$light_type[locRows] <-  "double"
    } else {
        bracerFileHam$light_type[locRows] <-  lightLoc
    }
}

table(bracerFileHam$light_type)
#K   L 
#618 422
#Which is great, as it means that no cells have more than one light chain now. 

#Now, we add mutational information
bracerFileHam <- observedMutations(bracerFileHam, sequenceColumn = "SEQUENCE_IMGT",
                                   germlineColumn = "GERMLINE_IMGT",
                                   regionDefinition = IMGT_V)

#We also add a column with the total number of mutations
bracerFileHam$All_mutations <- rowSums(bracerFileHam[,c("mu_count_cdr_r",
                                                    "mu_count_cdr_s",
                                                    "mu_count_fwr_r",
                                                    "mu_count_fwr_s"),])

bracerFileHam$Non_silent_mutations <- rowSums(bracerFileHam[,c("mu_count_cdr_r",
                                                        "mu_count_fwr_r"),])

#And a simple column denoting if the cell is clonal or not. 
bracerFileHam$Clonal <- FALSE
bracerFileHam$Clonal[-which(is.na(bracerFileHam$intraClonalDistance))] <- TRUE
write.csv(bracerFileHam, "Data/BCR_database_versions/4_IMGT_gapped_db_complete_post_clonality.csv", 
          row.names = FALSE)

#As we are now through all steps that involve thresholding etc, we will now
#exclude all the non-IgGs: 
IgG <- unique(bracerFileHam$CELL[grep("IGHG", bracerFileHam$ISOTYPE)])
bracerFileIgG <- bracerFileHam[which(bracerFileHam$CELL %in% IgG),]

##############
#FOR TABLE 1
uniqueHamCells <- unique(bracerFileHam$CELL)
length(unique(uniqueHamCells))-length(unique(bracerFileIgG$CELL))
#139
length(unique(bracerFileIgG$CELL))
#381
table(bracerFileIgG$Sample[which(bracerFileIgG$LOCUS == "H")])
#1166_1 1227_1 1284_1 
#212     65    104
#############
#And this is saved. 
write.csv(bracerFileIgG, "Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv", 
          row.names = FALSE)