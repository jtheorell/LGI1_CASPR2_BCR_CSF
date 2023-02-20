#We start by importing the cell cycle data
deepCycleRes <- read.csv("Data/DeepCycle/DeepCycle_csv_results/obs.csv")
deepCycleRes$Cell <- gsub("onefilepercell_JR...._._._..._and_others_.....:JR|.out.bam",
                          "", deepCycleRes$CellID)

aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

aeSce$cellCycle <- sapply(colnames(aeSce), function(x){
    if(x %in% deepCycleRes$Cell){
        deepCycleRes$cell_cycle_theta[which(deepCycleRes$Cell == x)]
    } else {
        NA
    }
})

library(ggplot2)
plotDat <- as.data.frame(colData(aeSce))

dir.create("Results/Cell_cycle")
ggplot(plotDat[which(plotDat$cellType == "ASC" & plotDat$Specific != "Not_tested"),], 
       aes(x = Specific, y = cellCycle, fill = Specific)) +
    geom_violin(alpha = 1) + geom_jitter(shape=16, position=position_jitter(0.1)) +
    scale_fill_manual(values = c("#555555", "orange")) + theme_bw()
ggsave("Results/Cell_cycle/Violin_only_tested.pdf")