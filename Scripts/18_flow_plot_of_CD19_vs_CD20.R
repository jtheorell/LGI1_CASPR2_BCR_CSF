library(SingleCellExperiment)
library(ggplot2)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")

aeSceB <- aeSce[,-which(aeSce$clustersLouvain == "CD8T")]

protPlotDat <- as.data.frame(t(normcounts(altExp(aeSceB, "flowData"))))
protPlotDat <- cbind(protPlotDat, colData(aeSceB))
protPlotDat$Specific[which(is.na(protPlotDat$Specific))] <- "Not_tested"
p <- ggplot(protPlotDat, aes(x=CD20, y=CD19, shape = donor)) + 
    geom_point() + xlim(-2, 6) + ylim(-2, 6) + theme_bw()
p
ggsave("Results/Flow_specific/CD20_vs_CD19.pdf")
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave("Results/Flow_specific/CD20_vs_CD19_stripped.pdf", width = 5, height = 5)


#And the same is done for a few more. 
#FOr CD3, we need to investigate what the high levels are, to adjust the scale correctly
which.max(normcounts(altExp(aeSce, "flowData"))["CD3",])

p <- ggplot(protPlotDat, aes(x=CD3, y=CD19, shape = donor)) + 
    geom_point() + xlim(-3, 8) + ylim(-2, 6) + theme_bw()
p
ggsave("Results/Flow_specific/CD3_vs_CD19.pdf")
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave("Results/Flow_specific/CD3_vs_CD19_stripped.pdf", width = 5, height = 5)


p <- ggplot(protPlotDat, aes(x=CD27, y=CD38, shape = donor)) + 
    geom_point() + theme_bw()
p
ggsave("Results/Flow_specific/CD27_vs_CD38.pdf")
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave("Results/Flow_specific/CD27_vs_CD38_stripped.pdf", width = 5, height = 5)

p <- ggplot(protPlotDat, aes(x=CD20, y=CD38, shape = donor)) + 
    geom_point() + xlim(-2, 6) + theme_bw()
p
ggsave("Results/Flow_specific/CD20_vs_CD38.pdf")
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave("Results/Flow_specific/CD20_vs_CD38_stripped.pdf", width = 5, height = 5)





