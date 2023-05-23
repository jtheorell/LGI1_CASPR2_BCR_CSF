#Here, we are going to plot a number of interesting factors on the clonal trees:
#CD20prot, CD20gene, CD27prot, CD27gene, CD38prot, CD38gene, HLA-DR genes, end-point dilution and 
#cell type
library(viridis)
library(dowser)
library(ggtree)
#Here, we are going to overlay a number of factors on the trees. 
trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")

#Before doing anything else, we need to rename these clones so that they also
#inform us about the donor and the specificity
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

trees$data <- lapply(trees$data, function(x){
    locClone <- x@clone
    locSample <- x@data$Sample[1]
    locSpecs <- x@data$Specific
    if(any(locSpecs == "TRUE")){
        locSpec <- "TRUE"
    } else if(any(locSpecs == "FALSE")){
        locSpec <- "FALSE"
    } else {
        locSpec <- "Not_tested"
    }
    
    x@clone <- paste0(locSample, "_", locSpec, "_", locClone)
    x
})

for(x in seq_along(trees$trees)){
    trees$trees[[x]]$name <- trees$data[[x]]@clone
}

names(trees$trees) <- unlist(lapply(trees$data, function(x) x@clone))


###########
#EPD
###########
epdDat <- read.csv("Data/BCR_auxiliaries/End-point_dilutions.csv")

#These values are now changed to something more meaningful, as a low value means
#a high concentration/avidity
epdDat$minusLog10EPD <- round(-log10(epdDat$EPD), 1)

#Here, we bin them: 
epdDat$minusLog10EPDBinned <- as.character(cut(epdDat$minusLog10EPD, 10, 
                                                 include.lowest = T))

trees$data <- lapply(trees$data, function(x){
    x@data$locDat <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- substr(y, 1,12)
        if(locName %in% epdDat$CELL){
            epdDat$minusLog10EPDBinned[which(epdDat$CELL == locName)]
        } else {
            "Not_present"
        }
    }))
    x
})

plots = plotTrees(trees, tips="locDat", tipsize = 5, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("(2.81,3.2]" = "#67001F",
                                                 "(2.42,2.81]" = "#B2182B",
                                                 "(2.03,2.42]" = "#D6604D",
                                                 "(1.64,2.03]" = "#F4A582",
                                                 "(1.25,1.64]" = "#FDDBC7",
                                                 "(0.86,1.25]" = "white",
                                                 "(0.47,0.86]" = "#D1E5F0",
                                                 "(0.08,0.47]" = "#92C5DE",
                                                 "(-0.31,0.08]" = "#4393C3",
                                                 "[-0.704,-0.31]" = "#2166AC",
                                                 "Germline" = "#000000FF",
                                                 "Not_present" = "grey"))

treesToPDF(plots, file="Results/Figure_4_plots/Trees/BCR_trees_with_EPD.pdf",
           nrow=2, ncol=2)

#As the data format is so extremely convoluted, I will add the colors to the 
#germline sequences myself. 

#We also add another factor here, namely the CV within and betwen all clones
epdDat$HamClone <- unlist(sapply(epdDat$CELL, function(x){
    if(x %in% colnames(aeSce)){
        aeSce$HamClone[which(colnames(aeSce) == x)]
    } else {
        NA
    }
    
}))
#Now, we calculate the coefficient of variation for the full shebang
fullCV <- sd(epdDat$EPD)/mean(epdDat$EPD)*100

#Zooming in on the originally tested clones
original_selection <- read.csv("Data/BCR_auxiliaries/Originally_selected_BCRs.csv")[,1]

epdDatTested <- epdDat[which(epdDat$CELL %in% original_selection ),]

#And here for the individual clones
isCloneLargeTested <- unlist(lapply(split(epdDatTested$HamClone, epdDatTested$HamClone), 
                                    function(x){
                                        if(length(x) > 1){
                                            TRUE
                                        } else {
                                            FALSE
                                        }
                                    }))

largeTested <- names(isCloneLargeTested)[which(isCloneLargeTested)]

clonalEPDDat <- epdDatTested[which(epdDatTested$HamClone %in% largeTested),]
#For some reason, 158 is not included here. 
clonalEPDDat <- rbind(clonalEPDDat, epdDat[which(epdDat$HamClone == 158),])
cvVec <- unlist(lapply(split(clonalEPDDat$EPD, clonalEPDDat$HamClone), function(x){
    sd(x) / mean(x) * 100
}))

names(cvVec) <- paste0("Clone_", names(cvVec))

cvVec$allCellCV <- fullCV
round(unlist(cvVec))

#Clone_5  Clone_83  Clone_84 Clone_158 Clone_202 Clone_212 Clone_214 Clone_291 allCellCV 
#0        58         0        85        46        87        75         0       438 
#Now, we are going to calculate the intraclonal distances to the germline and the 
#longest distance between any clonal members. 

########################
#ALT 1: germ to founder and max intraclonal

treeDists <- as.data.frame(do.call("rbind", lapply(seq_along(trees$trees), function(i){
    locDat <- trees$trees[[i]]
    #This includes the germline which is always the last number of these. 
    numTips <- length(locDat$tip.label)
    #First: distance from each tip to the central point of the tree, where all three
    #branches (the first branching point, and the stem from the germline) meet. 
    #This point is always two points away from the germline. 
    centralPoint <- numTips+2
    edgeDat <- data.frame(locDat$edge)
    edgeDat$length <- locDat$edge.length
    lengthToCent <- lapply(1:(numTips-1), function(x){
        centPointNotFound <- TRUE
        locNum <- x
        listVal <- 1
        locDist <- 0
        locResList <- list()
        while(centPointNotFound){
            locRow <- edgeDat[which(edgeDat$X2 == locNum),]
            locDist <- locDist+locRow$length
            locRes <- cbind(locRow[1:2],locDist)
            colnames(locRes) <- c("To_node", "From_node", "Acc_dist")
            locResList[[listVal]] <- locRes
            locNum <- as.numeric(locRow[1])
            listVal <- listVal+1
            if(locNum == centralPoint){
                centPointNotFound <- FALSE
            }
        }
        do.call("rbind", locResList)
        
    })
    #The length from the central point to the germline is essentially
    germToCent <- edgeDat$length[which(edgeDat$X2 == numTips)] + 
        edgeDat$length[which(edgeDat$X2 == centralPoint)]
    
    #Now, what is the lowest sum of germToCent plus one of the others?
    minLengthToCent <- min(unlist(lapply(lengthToCent, function(x){
        max(x$Acc_dist)
    })))
    #So the shortest distance is
    germToFounder <- germToCent+minLengthToCent
    
    #Now, what is the longest leaf-to-leaf distance? 
    maxIntraClonal <- max(unlist(lapply(lengthToCent, function(x){
        max(unlist(lapply(lengthToCent, function(y){
            locToNodeX <- x$To_node
            locToNodeY <- y$To_node
            xPos <- which(locToNodeX %in% locToNodeY)[1]
            yPos <- which(locToNodeY %in% locToNodeX)[1]
            #This provision is brought in to make sure that the longest distance
            #is not ever going to be from one cell to itself, due to a computational
            #oddity
            if(identical(x$Acc_dist, y$Acc_dist)){
                0
            }else{
                x$Acc_dist[xPos]+y$Acc_dist[yPos] 
            }
        })))
    })))
    outerRes <- c(germToFounder, maxIntraClonal)
    names(outerRes) <- c("germToFounder", "maxIntraClonal")
    outerRes
})))

treeDists$cloneNum <- sapply(trees$trees, function(x) x$name)
treeDists$size<- unlist(trees$seqs)
plot(treeDists$size, treeDists$maxIntraClonal)
plot(treeDists$size, treeDists$germToFounder)

cor.test(treeDists$size, treeDists$maxIntraClonal, method = "spearman")
#Not very unexpectedly, there is a very strong correlation here. 
#p-value = 4.638e-07
cor.test(treeDists$size, treeDists$maxIntraClonal, method = "pearson")
#The Pearson correlation coefficient is 0.658497

#And this is exported for the UCA delta comparison analyses. 
write.csv(treeDists, "Results/Figure_4_plots/Trees/treeDists.csv", row.names = FALSE)

########################
#ALT 1: germ to founder and max intraclonal

treeDists2 <- as.data.frame(do.call("rbind", lapply(seq_along(trees$trees), function(i){
    locDat <- trees$trees[[i]]
    #This includes the germline which is always the last number of these. 
    numTips <- length(locDat$tip.label)
    #First: distance from each tip to the central point of the tree, where all three
    #branches (the first branching point, and the stem from the germline) meet. 
    #This point is always two points away from the germline. 
    centralPoint <- numTips+2
    edgeDat <- data.frame(locDat$edge)
    edgeDat$length <- locDat$edge.length
    lengthToCent <- lapply(1:(numTips-1), function(x){
        centPointNotFound <- TRUE
        locNum <- x
        listVal <- 1
        locDist <- 0
        locResList <- list()
        while(centPointNotFound){
            locRow <- edgeDat[which(edgeDat$X2 == locNum),]
            locDist <- locDist+locRow$length
            locRes <- cbind(locRow[1:2],locDist)
            colnames(locRes) <- c("To_node", "From_node", "Acc_dist")
            locResList[[listVal]] <- locRes
            locNum <- as.numeric(locRow[1])
            listVal <- listVal+1
            if(locNum == centralPoint){
                centPointNotFound <- FALSE
            }
        }
        do.call("rbind", locResList)
        
    })
    #The length from the central point to the germline is essentially
    germToCent <- edgeDat$length[which(edgeDat$X2 == numTips)] + 
        edgeDat$length[which(edgeDat$X2 == centralPoint)]
    
    #Now, what is the lowest sum of germToCent plus one of the others?
    founder <- which.min(unlist(lapply(lengthToCent, function(x){
        max(x$Acc_dist)
    })))
    minLengthToCent <- unlist(lapply(lengthToCent, function(x){
        max(x$Acc_dist)
    }))[founder]
    #So the shortest distance is
    germToFounder <- germToCent+minLengthToCent
    
    #Now, what is the longest distance from the founder?
    furthest <- which.max(unlist(lapply(lengthToCent, function(x){
        max(x$Acc_dist)
    })))
    
    
    maxToFounder <- if(founder == furthest){
        0
    } else {
        foundDat <- lengthToCent[[founder]]
        furthDat <- lengthToCent[[furthest]]
        foundToNode <- foundDat$To_node
        furthToNode <- furthDat$To_node
        foundPos <- which(foundToNode %in% furthToNode)[1]
        furthPos <- which(furthToNode %in% foundToNode)[1]
        foundDat$Acc_dist[foundPos]+furthDat$Acc_dist[furthPos] 
    }
    outerRes <- c(germToFounder, maxToFounder)
    names(outerRes) <- c("germToFounder", "maxToFounder")
    outerRes
})))

treeDists2$cloneNum <- sapply(trees$trees, function(x) x$name)
treeDists2$size<- unlist(trees$seqs)
plot(treeDists2$size, treeDists2$maxToFounder)
plot(treeDists2$size, treeDists2$germToFounder)

cor.test(treeDists2$size, treeDists2$maxToFounder, method = "spearman")
#Not very unexpectedly, there is a very strong correlation here. 
#p-value = 3.39e-06
cor.test(treeDists2$size, treeDists2$maxToFounder, method = "pearson")
#The Pearson correlation coefficient is 0.5386276

#And this is exported for the UCA delta comparison analyses. 
write.csv(treeDists2, "Results/Figure_4_plots/Trees/treeDists_founder_to_furthest.csv", row.names = FALSE)

#Here, we are also going to print the clonal members for four shown clones and their
#respective cell subtype according to the developmental model. The reason for this
#is that this information is going to be manually added to the figure. 
aeScePCDat <- readRDS("Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")
focClones <- c(5, 83, 84, 158, 202, 212, 214, 291)
focPCDat <- as.data.frame(colData(aeScePCDat)[which(aeScePCDat$HamClone %in% focClones),c("llpcSingler", "HamClone")])
focPCDatOrd <- focPCDat[order(focPCDat$HamClone, focPCDat$llpcSingler),]

#             llpcSingler HamClone
#1166_1_2_H14          PB        5
#1166_1_2_J02          PB        5
#1166_1_2_N14          PC        5
#1166_1_1_H01          PC       83
#1166_1_2_E06          PC       83
#1166_1_2_J04          PC       83
#1166_1_2_F13          PB       84
#1166_1_2_G13          PB       84
#1166_1_2_P12          PB       84
#1166_1_2_P15          PB       84
#1166_1_2_I08          PB      158
#1166_1_2_M08          PB      158
#1166_1_2_N02          PB      158
#1166_1_2_D03          PB      202
#1166_1_2_E05          PB      202
#1166_1_2_F12          PB      202
#1166_1_2_I11          PB      202
#1166_1_2_L14          PB      202
#1166_1_2_N05          PB      202
#1166_1_2_J03          PC      202
#1166_1_2_D09          PB      212
#1166_1_2_E15          PB      212
#1166_1_2_J12          PB      212
#1166_1_2_N10          PC      212
#1166_1_2_G08       prePB      214
#1166_1_2_L16       prePB      214
#1166_1_2_H09          PB      214
#1166_1_2_I16          PB      214
#1166_1_2_M05          PB      214
#1166_1_2_F14         MBC      291
#1166_1_2_J07       prePB      291
#1166_1_2_J16       prePB      291

#we also need information on the UCA binding of clone 5, to be able to manually
#change that color. 
EPD_UCA <- read.csv("Data/BCR_auxiliaries/UCA_EPD_and_mut.csv")
round(-log10(EPD_UCA$EPD_.ug.mL._UCA_BS[which(EPD_UCA$CELL == "1166_1_2_H14")])[1],1)
#1.7

#Now after looking at the trees, it is inevitable to want to look at cell type vs
#EPD and UCA binding. 
aeScePCDatEPD <- aeScePCDat[,which(colnames(aeScePCDat) %in% epdDat$CELL)]
epdDatCELL <- epdDat[which(epdDat$CELL %in% colnames(aeScePCDatEPD)),]
epdDatCELLOrd <- epdDatCELL[match(colnames(aeScePCDatEPD), epdDatCELL$CELL),]
identical(epdDatCELLOrd$CELL, colnames(aeScePCDatEPD)) #TRUE
epdDatCELLOrd$cellType <- aeScePCDatEPD$llpcSingler
epdDatCELLOrd$cloneSize <- aeScePCDatEPD$ClonalSize
epdDatCELLOrd$Clonal <- aeScePCDatEPD$Clonal
epdDatCELLOrd <- epdDatCELLOrd[-which(is.na(epdDatCELLOrd$cellType)),]

library(ggforce)
epdDatCELLOrdSmall <- epdDatCELLOrd[which(epdDatCELLOrd$cloneSize != "≥4"),]
ggplot(epdDatCELLOrdSmall, aes(x = cellType, y = minusLog10EPD)) +
     geom_violin() + geom_sina(jitter_y = FALSE) + theme_bw()

epdDatCELLOrdLarge <- epdDatCELLOrd[which(epdDatCELLOrd$cloneSize == "≥4"),]
ggplot(epdDatCELLOrdLarge, aes(x = cellType, y = minusLog10EPD)) +
    geom_violin() + geom_sina(jitter_y = FALSE) + theme_bw()

#This is for all. Now we split by clonal size. 


#Further, what about UCA binding strength?
epdDatCELLOrd$UCA_EPC <- sapply(epdDatCELLOrd$CELL, function(x){
    if(x %in% EPD_UCA$CELL){
        EPD_UCA$EPD_.ug.mL._UCA_BS[which(EPD_UCA$CELL == x)][1]
    } else {
        100
    }
})
epdDatCELLOrd$UCAPos <- FALSE
epdDatCELLOrd$UCAPos[which(epdDatCELLOrd$UCA_EPC < 100)] <- TRUE

#        MBC prePB PB PC
#FALSE   1     6 53 17
#TRUE    0     0 13  5
#Of course there is no difference here. 

