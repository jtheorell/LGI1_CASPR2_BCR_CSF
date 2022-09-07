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
#1504   59
bracerFile <- bracerFile[-which(bracerFile$CELL %in% BCR_contaminations),]
dim(bracerFile)
#1450   59

#Here, we remove all non-functional chains
table(bracerFile$FUNCTIONAL)
#FALSE  TRUE 
#244  1203

bracerFile <- bracerFile[which(bracerFile$FUNCTIONAL),]
dim(bracerFile)
#1203   59

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
#This makes us lose 67 incomplete cells
uniqueCells <- uniqueCells[which(fullCells)]

write.csv(bracerFile[which(bracerFile$CELL %in% uniqueCells),], "Data/BCR_database_versions/1_IMGT_gapped_db_TPM_added.csv", 
          row.names =FALSE)

#Here cells with more than one functional H and K or L that have an edit distance
#above 10 are removed, as such cells are highly likely to be doublets, and thus
#should be excluded from all downstream analyses

doubletCells <- unlist(lapply(uniqueCells, function(x){
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
uniqueCells <- uniqueCells[-which(doubletCells)]

bracerFile <- bracerFile[which(bracerFile$CELL %in% uniqueCells),]
dim(bracerFile)
#1127   59

##############
#FOR TABLE 1
length(unique(bracerFile$CELL))
#This gives us 520 in total, and in other words, we have lost 83 cells due to 
#incomplete chains. 
##############

###############
#Exclusions of cells of dubious character. 

#Exclusion of cells with a TCR. TraCeR has been used to identify these: 
tcrContainers <- read.csv("Data/BCR_auxiliaries/tcr_bcr_doublets.csv")[,2]

bracerFile <- bracerFile[-which(bracerFile$CELL %in% tcrContainers),]
#In total, this makes us lose 10 cells. 

#Now, we import information about the cell type, and exclude the discordant cells
#i.e., cells that carry a BCR but that does not have a B-cell phenotype. 
cellTypeData<- read.csv("Data/flowDataPlusIndexAndcellType.csv", row.names = 1)

bracerFile$Cell_type <- NA
for(i in unique(bracerFile$CELL)){
    bracerFile$Cell_type[which(bracerFile$CELL == i)] <- 
        cellTypeData$Cell_type[which(cellTypeData$Cell == i)]
}

table(bracerFile$Cell_type[which(bracerFile$LOCUS == "H")])
#ASC    B CD4T CD8T 
#288  229    4    3

#So these 7 non-B are excluded here. 
bracerFile <- bracerFile[-grep("T", bracerFile$Cell_type),]

write.csv(bracerFile, "Data/BCR_database_versions/2_IMGT_gapped_db_all_dubious_excluded.csv", 
          row.names =FALSE)

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

#Now, we are down to 1560 chains from the remaining, complete cells. 
#And this is saved
write.csv(bracerFile3Chains, "Data/BCR_database_versions/3_IMGT_gapped_db_TPM_added_3_chains.csv", 
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
write.table(changeoClusteringFile, "Data/BCR_database_versions/changeoClusteringFile.tsv", 
            sep = "\t",row.names = FALSE)

system(paste0("DefineClones.py -d Data/BCR_database_versions/changeoClusteringFile.tsv -o Data/changeoClusteringFile_Ham_result.tsv --act set --model ham --norm len --format changeo --sf JUNCTION --vf V_CALL --jf J_CALL --dist ", output@threshold))

#Now, we re-import these. 
hamClones <- read.table("Data/BCR_database_versions/changeoClusteringFile_Ham_result.tsv", header = TRUE)

#At this stage, we actually exchange the clone numbers for a standardised set used
#in previous rounds of analysis. They are in all respects equal
#to the ones generated above, apart from the actual numbers. 
originalCloneNames <- read.csv("Data/BCR_auxiliaries/Original_clone_names.csv")

bracerFileHam <- bracerFile3Chains
bracerFileHam$HamClone <- NA

for(i in unique(bracerFileHam$CELL)){
    bracerFileHam$HamClone[which(bracerFileHam$CELL == i)] <- 
        originalCloneNames$HamClone[which(originalCloneNames$CELL == i)]
}
#There are warnings here, as we are using a single row as input for many rows,
#but that is correct. 

#The original code used to introduce the hamClone names was this: 
#uniqueCells <- unique(hamClones$CELL)
#
#bracerFileHam <- do.call("rbind", lapply(uniqueCells, function(x){
#    locCell <- bracerFile3Chains[which(bracerFile3Chains$CELL == x),]
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
#212  24   6  13   6   5   4

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
#    0   1   2   3   4   12
#   209  37  11   7   2  4

#We also un-name all clones that consist of one cell only
bracerFileHam$HamClone[which(is.na(bracerFileHam$intraClonalDistance))] <- NA

#And with that, we have a very clean dataset, with 503 cells, where the largest
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
#double      K      L 
#   87    586    362 

#Now, we add mutational information
bracerFileHam <- observedMutations(bracerFileHam, sequenceColumn = "SEQUENCE_IMGT",
                                   germlineColumn = "GERMLINE_IMGT",
                                   regionDefinition = IMGT_V)



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
length(unique(bracerFileHam$CELL[-which(bracerFileHam$CELL %in% 
                                            bracerFileIgG$CELL)]))
#139
length(unique(bracerFileIgG$CELL))
#365
table(bracerFileIgG$Sample[which(bracerFileIgG$LOCUS == "H")])
#1166_1 1227_1 1284_1 
#211     63     91
#############
#And this is saved. 
write.csv(bracerFileIgG, "Data/BCR_database_versions/5_IMGT_gapped_db_IgG.csv", 
          row.names = FALSE)

uniqueCells <- unique(bracerFileIgG$CELL)

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
BCR_H <- bracerFileIgG[which(bracerFileIgG$LOCUS == "H"),]
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

BCR_K <- bracerFileIgG[which(bracerFileIgG$LOCUS == "K"),]
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
BCR_L <- bracerFileIgG[which(bracerFileIgG$LOCUS == "L"),]
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

#And with this, we are ready to export the file and to export a file used for
#integration with the flow cytometry data
write.csv(bcrDf, "Data/BCR_database_versions/condensed_file_one_BCR_per_row.csv")
