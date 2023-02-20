#I am here heavily reliant on 
#https://www.antibodysociety.org/wordpress/wp-content/uploads/2021/11/Immcantation-webinar-slides.pdf
#What we however do here is that we, for the first part, move out to terminal
#and run immcantation through Docker. 
#########
#TERMINAL
#########

#docker run -it -v ~/Labbet/2022/220818_full_LGI1_B-cell_analysis/For_github/:/data:z immcantation/suite:4.3.0 bash
R
setwd("/data")
library(dowser)

BCR_clonal <- read.csv("Data/BCR_database_versions/6_Specificity_included.csv")
BCR_clonal$seq_id <- paste0(BCR_clonal$CELL, BCR_clonal$LOCUS)
clones = formatClones(BCR_clonal,id = "seq_id",
                      seq = "SEQUENCE_IMGT", 
                      germ = "GERMLINE_IMGT",
                      v_call = "V_CALL", j_call = "J_CALL", 
                      junc_len = "JUNCTION_LENGTH",
                      traits = c("Sample", "ISOTYPE", "Specific"),
                      clone = "HamClone", heavy = "H", cell = "CELL", 
                      locus = "LOCUS", add_count = TRUE, collapse = FALSE)

#The reason we are using Docker is to be able to access this igphyml, which
#is a tree-construction algorithm that takes all the quirks of B-cell mutational
#information into account when constructing the trees. 
trees = getTrees(clones, build="igphyml",
                 exec="/usr/local/share/igphyml/src/igphyml", nproc=7)

#And with that, we are done in the immcantation world and can go back here. 
saveRDS(trees, "Data/BCR_auxiliaries/IgphyML_trees.rds")

library(dowser)
library(ggtree)
BCR_all <- read.csv("Data/BCR_database_versions/6_Specificity_included.csv")

trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")
    
plots = plotTrees(trees, tips="ISOTYPE", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("IGHG4" ="#238A8D",
                                                 "IGHG1" = "#B1CFFF",
                                                 "IGHG2" = "#9E4532",
                                                 "Germline" = "black"))
dir.create("Results/Trees")
treesToPDF(plots, file="Results/Trees/BCR_trees.pdf",
           nrow=2, ncol=2)

#We also create trees with named leaves
testPlots <- lapply(plots, function(x) x + geom_tiplab())
treesToPDF(testPlots, file="Results/Trees/BCR_trees_w_labels.pdf",
           nrow=2, ncol=2)

#These data structures are completely crazy. Luckly there are others who understand
#so here is a link for an explainer: 
#https://mohaksharda.medium.com/understanding-the-structure-of-an-object-of-class-phylo-92102cf357a6
#Sadly, this explanation above, which does make perfect sense, does for
#some reason not apply here. Instead the order is very messed up. However, 
#there are definable relations between the leaves, i.e. the BCRs, and their
#edge lengths, so that a BCR that is part of a network of edges with a length
#of 0, are identical. This can be harnessed for downstream plotting.
#It will inevitably be a bit nested though. 
#We will now create a new column called "subclone" to be used for cells that are
#Identical. The thing is that in each clone there might be multiple such, so
#the old standard of just creating a intraClonalDistance thing, and merging
#all zeros is not complex enough.
BCR_all$SubClone <- NA
hamClones <- unique(BCR_all$HamClone)
hamClones <- hamClones[-which(is.na(hamClones))]
for(i in hamClones){
    locTree <- trees[which(trees$clone_id == i),]$trees[[1]]
    #This way of correctly identifying which cell is associated to which
    #number is taken from this conversation: 
    #https://stackoverflow.com/questions/34364660/how-to-get-correct-order-of-tip-labels-in-ape-after-calling-ladderize-function
    is_tip <- locTree$edge[,2] <= length(locTree$tip.label)
    ordered_tips <- locTree$edge[is_tip, 2]
    locTree$tip.label[ordered_tips]
    nodeCellDf <- data.frame("Tip_label" = ordered_tips,
                                    "Cell" = locTree$tip.label[ordered_tips])
    nodeCellDf$SubClone <- 0
    edgeDf <- data.frame(locTree$edge, locTree$edge.length)
    colnames(edgeDf) <- c("Other_edge_ends", "Tips_and_nodes", "Edge_length")
    #Now, we identify cells that have only zero-length edges to other cells.
    #With this particular model, however, there are no absolute zeros. 
    #for this reason, we introduce a cutoff here. It is arbitrary, but judging from
    #all the clonal representations, it gives identical results to the graphical
    #output. 
    edgeDf$Edge_length[which(edgeDf$Edge_length < 0.0001)] <- 0
    subCloneNum <- 1
    while(any(nodeCellDf$SubClone == 0)){
        locPos <- which(nodeCellDf$SubClone == 0)[1]
        locTip <- nodeCellDf$Tip_label[locPos]
        
        #And now to the hard part. 
        #FIist, we identify the row with the cell of interest.
        locEdgeRow <- which(edgeDf$Tips_and_nodes == locTip)
        #Then, we check if this cells edge to the next node has no length, 
        #and if so all other cells in connection with this node are identified
        #and clumped together in the same subclone
        
        if(edgeDf$Edge_length[locEdgeRow] == 0){
            stillSameSubclone <- TRUE
            zeroNodes <- edgeDf$Other_edge_ends[locEdgeRow]
            #Here we iteratively add more and more rows with zero edges, until
            #there are none more to be found. 
            while(stillSameSubclone){
                allZeroNodeRows <- edgeDf[which((edgeDf$Other_edge_ends %in% zeroNodes |
                                                     edgeDf$Tips_and_nodes %in% zeroNodes) &
                                                    edgeDf$Edge_length == 0),]
                allZeroNodes <- unique(unlist(allZeroNodeRows[,1:2]))
                if(identical(zeroNodes,allZeroNodes)){
                    
                    stillSameSubclone <- FALSE
                    
                } else {
                    zeroNodes <- allZeroNodes
                }
            }
            nodeCellDf$SubClone[which(nodeCellDf$Tip_label %in% zeroNodes)] <- 
                subCloneNum
            
        } else {
            nodeCellDf$SubClone[locPos] <- subCloneNum
        }
        
        subCloneNum <- subCloneNum +1
    }
    #And this information is returned to the outer world. We need to shave off
    #the added H to the cell names though. 
    for(j in 1:nrow(nodeCellDf)){
        locCell <- substr(nodeCellDf$Cell[j], 1, 12)
        BCR_all$SubClone[which(BCR_all$CELL == locCell)] <- 
            nodeCellDf$SubClone[j]
    }
}

#And there we are done and ready for the next step!
write.csv(BCR_all, "Data/BCR_database_versions/7_Subclones_included.csv",
          row.names = FALSE)