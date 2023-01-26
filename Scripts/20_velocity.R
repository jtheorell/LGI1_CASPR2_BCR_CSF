#This analyis is a small endpoint of a vast dive into the velocities. See separate .sh and .py docs. 
#In this case, we are going to include the B-cells. 
library(SingleCellExperiment)
library(ggplot2)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")

veloRes <- read.csv("Data/Velocity/velocyto_out/All_cells/obs.csv")

nameVec <- gsub(".+_others_.....:JR|\\.out\\.bam", "",
                veloRes$CellID)

veloPseudo <- data.frame("Cell" = nameVec, "pseudotime" = veloRes$velocity_pseudotime)

csfSce$pseudotime <- sapply(colnames(csfSce), function(x){
  if(x %in% veloPseudo$Cell){
      veloPseudo$pseudotime[which(veloPseudo$Cell == x)]
  } else {
      NA
  }
})

summary(csfSce$pseudotime[which(csfSce$Specific == "TRUE")])
summary(csfSce$pseudotime[which(csfSce$Specific == "FALSE")])

wilcox.test(csfSce$pseudotime[which(csfSce$Specific == "TRUE")], csfSce$pseudotime[which(csfSce$Specific == "FALSE")])
#p-value = 0.001061
#The difference is of course highly significant, as the specific cluster is almost
#exclusively ASC and the non-specific are all over the place. 
nonNASce <- csfSce[,-which(is.na(csfSce$pseudotime))]
nonNASce$Specific[which(is.na(nonNASce$Specific))] <- "Not_tested"

densDat <- data.frame("Specificity" = nonNASce$Specific,
                      "Pseudotime" = nonNASce$pseudotime)
ggplot(densDat, aes(x=Pseudotime, fill=Specificity)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_specificity_with_B-cells.pdf")
densDat <- data.frame("Specificity" = nonNASce$Specific[which(nonNASce$clustersLouvain == "ASC")],
                      "Pseudotime" = nonNASce$pseudotime[which(nonNASce$clustersLouvain == "ASC")])
ggplot(densDat, aes(x=Pseudotime, fill=Specificity)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_specificity_ASC.pdf")


#THis is a slight violation, as a few of them theoretically couldbe clonal. 
nonNASce$Clonal[which(is.na(nonNASce$Clonal))] <- FALSE
densDat <- data.frame("Clonality" = nonNASce$Clonal,
                      "Pseudotime" = nonNASce$pseudotime)
ggplot(densDat, aes(x=Pseudotime, fill=Clonality)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_clonality_with_B-cells.pdf")

densDat <- data.frame("Clonality" = nonNASce$Clonal[which(nonNASce$clustersLouvain == "ASC")],
                      "Pseudotime" = nonNASce$pseudotime[which(nonNASce$clustersLouvain == "ASC")])
ggplot(densDat, aes(x=Pseudotime, fill=Clonality)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_clonality_ASC.pdf")


#We of course also need to overlay this on the hierarchical trees.
trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")

#Now, we bin the pseudotime. 
floatPseudo <- veloPseudo$pseudotime
definedPseudo <- .bincode(floatPseudo, seq(0,1, by = 0.1))
veloPseudo$binnedPseudo <- as.factor(definedPseudo)
library(viridis)

trees$data <- lapply(trees$data, function(x){
    x@data$pseudotime <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- gsub("|H", "", y)
        if(locName %in% veloPseudo$Cell){
            veloPseudo$binnedPseudo[which(veloPseudo$Cell == locName)]
        } else {
            "Germline"
        }
    }))
    x
})

plots = plotTrees(trees, tips="pseudotime", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("10" ="#FCFFA4FF",
                                                 "9" = "#F6D645FF",
                                                 "8" = "#FCA50AFF",
                                                 "7" = "#F3771AFF",
                                                 "6" = "#DD513AFF",
                                                 "5" ="#BB3754FF",
                                                 "4" = "#932667FF",
                                                 "3" = "#6B186EFF",
                                                 "2" = "#420A68FF",
                                                 "1" = "#170C3AFF",
                                                 "Germline" = "#000004FF"))

treesToPDF(plots, file="Results/Trees/BCR_trees_with_velo_pseudotime.pdf",
           nrow=2, ncol=2)

#Now, do trees with high intraclonal celocuty diversity differ in specificity from the others?
cloneSce <- csfSce[,which(csfSce$Clonal)]

perCloneVeloDivergence <- unlist(lapply(split(cloneSce$pseudotime, cloneSce$HamClone), function(x){
    if(any(is.na(x))){
        locRange <- range(x[-which(is.na(x))])
    } else {
        locRange <- range(x)
    }
    locRange[2]-locRange[1]
}))

densDat <- data.frame("Specificity" = unlist(lapply(split(cloneSce$Specific, cloneSce$HamClone), "[[",1)),
                      "Pseudotime_divergence" = perCloneVeloDivergence)
densDat <- densDat[-which(densDat$Pseudotime_divergence < 0),]
ggplot(densDat, aes(y=Pseudotime_divergence, x = Specificity, fill=Specificity)) +
    geom_dotplot(binaxis='y', stackdir='center') + theme_bw()
ggsave("Results/Velocity/Intraclonal_pseudotime_divergence.pdf", width = 8, height = 5)                            

