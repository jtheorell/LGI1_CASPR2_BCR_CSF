#This analyis is a small endpoint of a vast dive into the velocities. See separate .sh and .py docs. 
library(SingleCellExperiment)
library(ggplot2)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

veloRes <- read.csv("Data/Velocity/All_cells/obs.csv")

nameVec <- gsub("onefilepercell_JR...._1_._..._and_others_.....:JR|\\.out\\.bam", "",
                veloRes$CellID)

veloPseudo <- data.frame("Cell" = nameVec, "pseudotime" = veloRes1166$velocity_pseudotime)

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
#p-value = 0.1628
#So the difference is not significant. 

densDat <- data.frame("Specificity" = csfSce$Specific,
                      "Pseudotime" = csfSce$pseudotime)
ggplot(densDat, aes(x=Pseudotime, fill=Specificity)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_specificity.pdf")

densDat <- data.frame("Specificity" = csfSce$Specific[which(csfSce$Cell_type == "ASC")],
                      "Pseudotime" = csfSce$pseudotime[which(csfSce$Cell_type == "ASC")])
ggplot(densDat, aes(x=Pseudotime, fill=Specificity)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_specificity_ASC.pdf")

#We of course also need to overlay this on the hierarchical trees.
trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")

#Now, we bin the pseudotime. 
floatPseudo <- veloPseudo$pseudotime
definedPseudo <- .bincode(floatPseudo, seq(0,1, by = 0.1))
veloPseudo$binnedPeudo <- as.factor(definedPseudo)
library(viridis)

trees$data <- lapply(trees$data, function(x){
    x@data$pseudotime <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- gsub("|H", "", y)
        if(locName %in% veloPseudo$Cell){
            veloPseudo$binnedPeudo[which(veloPseudo$Cell == locName)]
        } else {
            "Germline"
        }
    }))
    x
})

plots = plotTrees(trees, tips="pseudotime", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("10" ="#FFEA46FF",
                                                 "9" = "#E4CF5BFF",
                                                 "8" = "#C4B56CFF",
                                                 "7" = "#A69D75FF",
                                                 "6" = "#8A8779FF",
                                                 "5" ="#707173FF",
                                                 "4" = "#575C6DFF",
                                                 "3" = "#39486BFF",
                                                 "2" = "#00336FFF",
                                                 "1" = "#00204DFF",
                                                 "Germline" = "black"))

treesToPDF(plots, file="Results/Trees/BCR_trees_with_velo_pseudotime.pdf",
           nrow=2, ncol=2)

#My overarching impression is that these cells are very early in the progression. 
#So how do the non-clonal look?

densDat <- data.frame("Clonality" = csfSce$Clonal,
                      "Pseudotime" = csfSce$pseudotime)
ggplot(densDat, aes(x=Pseudotime, fill=Clonality)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Velocity/Velocity_vs_clonality.pdf")

#No, that was not the case. Most are early in this progression. 
#However, with only 8% of the trancripts being unspliced, it is likely that the
#cells generally have 


