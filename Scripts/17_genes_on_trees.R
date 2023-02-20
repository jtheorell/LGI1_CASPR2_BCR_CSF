#Here, we are going to plot a number of interesting factors on the clonal trees:
#CD20prot, CD20gene, CD27prot, CD27gene, CD38prot, CD38gene, HLA-DR genes, end-point dilution and 
#cell type
library(viridis)
library(dowser)
#Here, we are going to overlay a number of factors on the trees. 
trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")

#Before doing anything else, we need to rename these clones so that they also
#inform us about the donor and the specificity
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")

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

#############
#PROTEIN
#############
#Now, we create a neat dataset with this information. 
csfIndexDf <- read.csv("Data/Cytometry/flowDataPlusIndex.csv")

#To be able to plot this, we need to bin the data
treePlotBinList <- lapply(c("CD20", "CD27", "CD38"), function(x){
    locDat <- csfIndexDf[,which(colnames(csfIndexDf) == x)]
    as.factor(cut(locDat, 10, include.lowest = T, labels = 0:9))
})

names(treePlotBinList) <- c("CD20", "CD27", "CD38")

#And here it comes!
for(i in names(treePlotBinList)){
    locColDat <- treePlotBinList[[i]]
    trees$data <- lapply(trees$data, function(x){
        x@data$locDat <- factor(unlist(sapply(x@data$sequence_id, function(y){
            locName <- substr(y, 1,12)
            if(locName %in% csfIndexDf$Cell){
                locColDat[which(csfIndexDf$Cell == locName)]
            } else {
                "Not_present"
            }
        })), levels = c("Not_present", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
        x
    })
    plots = plotTrees(trees, tips="locDat", tipsize = 3, common_scale = TRUE,
                      labelsize = 8, tip_palette = c("11" ="#FFFFFF",
                                                     "10" ="#FFFFFF",
                                                     "9" ="#FCFFA4FF",
                                                     "8" = "#F6D645FF",
                                                     "7" = "#FCA50AFF",
                                                     "6" = "#F3771AFF",
                                                     "5" = "#DD513AFF",
                                                     "4" ="#BB3754FF",
                                                     "3" = "#932667FF",
                                                     "2" = "#6B186EFF",
                                                     "1" = "#420A68FF",
                                                     "0" = "#170C3AFF",
                                                     "Germline" = "#000000FF",
                                                     "Not_present" = "grey"))
    
    treesToPDF(plots, file=paste0("Results/Trees/BCR_trees_with_", i, ".pdf"),
               nrow=2, ncol=2)
}

#############
#RNA
#############
#Now, we create a neat dataset with this information. 
aeSceB <- aeSce[,-which(aeSce$clustersLouvain == "CD8T")]
rnaFocRows <- which(rowData(aeSceB)$hgnc_symbol %in% c("CD38", "CD27", "MS4A1", "HLA-DRB1", "HLA-DRA", "HLA-DRB5"))
rnaFocDat <- t(logcounts(aeSceB[rnaFocRows,]))
colnames(rnaFocDat) <- paste0("RNA_", rowData(aeSceB)$hgnc_symbol[rnaFocRows])

#To be able to plot this, we need to bin the data
treePlotBinList <- lapply(colnames(rnaFocDat), function(x){
    locDat <- rnaFocDat[,which(colnames(rnaFocDat) == x)]
    as.factor(cut(locDat, 10, include.lowest = T, labels = 0:9))
})

names(treePlotBinList) <- colnames(rnaFocDat)

#And here it comes!
for(i in names(treePlotBinList)){
    locColDat <- treePlotBinList[[i]]
    trees$data <- lapply(trees$data, function(x){
        x@data$locDat <- factor(unlist(sapply(x@data$sequence_id, function(y){
            locName <- substr(y, 1,12)
            if(locName %in% rownames(rnaFocDat)){
                locColDat[which(rownames(rnaFocDat) == locName)]
            } else {
                "Not_present"
            }
        })), levels = c("Not_present", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9"))
        x
    })
    plots = plotTrees(trees, tips="locDat", tipsize = 3, common_scale = TRUE,
                      labelsize = 8, tip_palette = c("11" ="#FFFFFF",
                                                     "10" ="#FFFFFF",
                                                     "9" ="#FCFFA4FF",
                                                     "8" = "#F6D645FF",
                                                     "7" = "#FCA50AFF",
                                                     "6" = "#F3771AFF",
                                                     "5" = "#DD513AFF",
                                                     "4" ="#BB3754FF",
                                                     "3" = "#932667FF",
                                                     "2" = "#6B186EFF",
                                                     "1" = "#420A68FF",
                                                     "0" = "#170C3AFF",
                                                     "Germline" = "#000000FF",
                                                     "Not_present" = "grey"))
    
    treesToPDF(plots, file=paste0("Results/Trees/BCR_trees_with_", i, ".pdf"),
               nrow=2, ncol=2)
}

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

plots = plotTrees(trees, tips="locDat", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("(2.81,3.2]" ="#FCFFA4FF",
                                                 "(2.42,2.81]" = "#F6D645FF",
                                                 "(2.03,2.42]" = "#FCA50AFF",
                                                 "(1.64,2.03]" = "#F3771AFF",
                                                 "(1.25,1.64]" = "#DD513AFF",
                                                 "(0.86,1.25]" ="#BB3754FF",
                                                 "(0.47,0.86]" = "#932667FF",
                                                 "(0.08,0.47]" = "#6B186EFF",
                                                 "(-0.31,0.08]" = "#420A68FF",
                                                 "[-0.704,-0.31]" = "#170C3AFF",
                                                 "Germline" = "#000000FF",
                                                 "Not_present" = "grey"))

treesToPDF(plots, file="Results/Trees/BCR_trees_with_EPD.pdf",
           nrow=2, ncol=2)

#As the data format is so extremely convoluted, I will add the colors to the 
#germline sequences myself. 


###########
#Cell type
###########

trees$data <- lapply(trees$data, function(x){
    x@data$locDat <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- substr(y, 1,12)
        if(locName %in% colnames(aeSceB)){
            aeSceB$clustersLouvain[which(colnames(aeSceB) == locName)]
        } else {
            "Not_present"
        }
    }))
    x
})

plots = plotTrees(trees, tips="locDat", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c(
                                                 "B" = "#FCA50AFF",
                                                 
                                                 "ASC" = "#BB3754FF",
                                                 "Germline" = "#000000FF",
                                                 "Not_present" = "grey"))

treesToPDF(plots, file="Results/Trees/BCR_trees_with_cellType.pdf",
           nrow=2, ncol=2)



#We also import the new deepcycle data

deepDat <- read.csv("Data/DeepCycle/DeepCycle_csv_results/obs.csv")

deepDat$Cell <- sapply(deepDat$CellID, function(x) gsub("onefilepercell_JR1166_1_1_B02_and_others_N8PVQ:|.out.bam", "", x))


trees$data <- lapply(trees$data, function(x){
    x@data$locDat <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- substr(y, 1,12)
        if(locName %in% deepDat$Cell){
            deepDat$cell_cycle_theta[which(deepDat$Cell == locName)]
        } else {
            "Not_present"
        }
    }))
    x
})

















cloneSce <- csfSceB[,which(csfSceB$Clonal)]

perCloneSlingDivergence <- unlist(lapply(split(cloneSce$slingPseudotime_1, cloneSce$HamClone), function(x){
    if(any(is.na(x))){
        locRange <- range(x[-which(is.na(x))])
    } else {
        locRange <- range(x)
    }
    locRange[2]-locRange[1]
}))

densDat <- data.frame("Specificity" = unlist(lapply(split(cloneSce$Specific, cloneSce$HamClone), "[[",1)),
                      "Pseudotime_divergence" = perCloneSlingDivergence)
ggplot(densDat, aes(y=Pseudotime_divergence, x = Specificity, fill=Specificity)) +
    geom_dotplot(binaxis='y', stackdir='center') + theme_bw()
ggsave("Results/Slingshot/Intraclonal_pseudotime_divergence.pdf", width = 8, height = 5)                            


