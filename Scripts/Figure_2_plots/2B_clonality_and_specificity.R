#Here, we remake the plot for 2B, as Ruby clearly shows something else. 
library(ggplot2)
BCR_all <- read.csv("Data/BCR_database_versions/8_new_mutations_added.csv")
BCR_h <- BCR_all[which(BCR_all$LOCUS == "H" & BCR_all$Specific != "Not_tested"),]
BCR_h$ClonalSize <- 1
BCR_h$ClonalSize[which(BCR_h$Clonal)] <- sapply(BCR_h$CELL[which(BCR_h$Clonal)], function(x){
    locCloneSize <- length(which(BCR_h$HamClone == BCR_h$HamClone[which(BCR_h$CELL == x)]))
    if(locCloneSize >= 4){
        "≥4"
    } else if(locCloneSize >= 2){
        "2-3"
    } else {
        as.character(locCloneSize)   
        }
})

clonSpecMat <- as.matrix(table(BCR_h$ClonalSize, BCR_h$Specific))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "Specific", "Freq")
clonSpecDf$CloneSize <- factor(clonSpecDf$CloneSize, levels = c("1", "2-3", "≥4"))
ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = Specific)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("black", "orange")) + 
    scale_y_continuous(limits = c(0, 80), expand = c(0, 0, 0, 0))
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity.pdf", width = 4.5, height = 4)

#And the percentages: 
round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))
#  nnnFALSE TRUE
#1      38   62
#2-3    21   79
#≥4      0  100

#And we also create one for the supplement with the percentages of other neuronal 
#tissue-specific cells in it. 
BCR_h$SpecTiss <- BCR_h$Specific
BCR_h$SpecTiss[which(BCR_h$CELL %in% c("1166_1_2_H12",
                                             "1166_1_2_A15",
                                             "1227_1_2_J08", 
                                             "1227_1_2_O05",
                                             "1166_1_2_J05", 
                                             "1166_1_2_L15",
                                             "1166_1_2_C11", 
                                             "1166_1_2_H08", 
                                             "1166_1_2_N07"))] <- "Tissue"

clonSpecMat <- as.matrix(table(BCR_h$ClonalSize, BCR_h$SpecTiss))
clonSpecDf <- as.data.frame(clonSpecMat)
colnames(clonSpecDf) <- c("CloneSize", "Specific_incl_tissue", "Freq")
clonSpecDf$CloneSize <- factor(clonSpecDf$CloneSize, levels = c("1", "2-3", "≥4"))
ggplot(clonSpecDf, aes(x = CloneSize, y = Freq, fill = Specific_incl_tissue)) + 
    geom_bar(stat="identity") + theme_bw() + 
    scale_fill_manual(values = c("black", "yellow", "orange")) + 
    scale_y_continuous(limits = c(0, 80), expand = c(0, 0, 0.05, 0))
ggsave("Results/Figure_2_plots/2B_clonality_vs_specificity_with_tissue.pdf", width = 4.5, height = 4)

#And the percentages: 
round(100*t(apply(clonSpecMat, 1, function(x) x/sum(x))))
#      FALSE Tissue TRUE
#1      33      5   62
#2-3    12      9   79
#≥4      0      0  100


