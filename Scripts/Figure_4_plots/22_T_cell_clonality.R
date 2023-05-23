#Here, we are going to clarify what rates of clonality we are seeing among the CD4 T cells
TCRdb <- read.csv("Data/TCR_data/cell_data_all_together.csv")
TCRdb$Sample <- substr(TCRdb$cell_name, 3,8)
TCRdb$Cell <- substr(TCRdb$cell_name, 3,16)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")


#First, we exclude all that are A1 wells, as these are positive control mini-bulks, 
#or that are present among the B-cells.
nrow(TCRdb) #1066
TCRdbRed <- TCRdb[-which(TCRdb$cell_name %in% colnames(aeSce) | grepl("A01|X",TCRdb$cell_name)),]
nrow(TCRdbRed) #1054
#Here, we do not exclude possibly incomplete cells, as we are interested in getting an
#as rich count as possible. 

TCRdbRed$newClone <- paste0(TCRdbRed$Sample, TCRdbRed$clonal_group)
cloneTab <- as.data.frame(table(TCRdbRed$newClone))
#And the non-real clones are excluded
cloneTab <- cloneTab[-which(cloneTab$Freq == 1 | grepl("NA", cloneTab$Var1)),]
TCRdbRed$Clonal <- FALSE
TCRdbRed$Clonal[which(TCRdbRed$newClone %in% cloneTab$Var1)] <- TRUE

#So now, we are going to introduce the CD4/CD8 division here. 
flowData <- read.csv("Data/Cytometry/flowDataPlusIndex.csv")

flowDataTCR <- flowData[which(flowData$Cell %in% TCRdbRed$Cell),]

plot(flowDataTCR$CD4, flowDataTCR$CD8)

flowDataTCR$T_cell_subset <- "DNT"
flowDataTCR$T_cell_subset[which(flowDataTCR$CD4 > 2.2 & flowDataTCR$CD8 > 2.5)] <- "DPT"
flowDataTCR$T_cell_subset[which(flowDataTCR$CD4 > 2.2 & flowDataTCR$CD8 <= 2.5)] <- "CD4T"
flowDataTCR$T_cell_subset[which(flowDataTCR$CD4 <= 2.2 & flowDataTCR$CD8 > 2.5)] <- "CD8T"

library(ggplot2)

ggplot(flowDataTCR, aes(x = CD4, y = CD8, color = T_cell_subset)) +
    geom_point() + theme_bw()

#This makes sense and can therefore be added
flowDataTCROrdered <- flowDataTCR[match(TCRdbRed$Cell, flowDataTCR$Cell),]
identical(flowDataTCROrdered$Cell, TCRdbRed$Cell)
#TRUE
TCRdbRed$Subset <- flowDataTCROrdered$T_cell_subset

#We are focusing on the CD4 compartment here. 
TCRdbCD4 <- TCRdbRed[which(TCRdbRed$Subset == "CD4T"),]

#Now, we inflate the clones
cloneDbAll <- TCRdbCD4[,c("Sample", "newClone", "Clonal")]
cloneDbAll$Size <- 1
nonClonalDb <- cloneDbAll[which(cloneDbAll$Clonal == FALSE),]
clonalDb <- cloneDbAll[which(cloneDbAll$Clonal),]
cloneDbCollapsed <- do.call("rbind", lapply(unique(clonalDb$newClone), function(x){
    locClone <- clonalDb[which(clonalDb$newClone == x),]
    resClone <-  locClone[1,]
    resClone$Size <- nrow(locClone)
    resClone
}))

fullCollapsedDb <- rbind(nonClonalDb, cloneDbCollapsed)

donList <- split(fullCollapsedDb, f = fullCollapsedDb$Sample)

#Here, we create a scaling factor to be used below
sizeVals <- unlist(lapply(donList, nrow))
fracVals <- sizeVals/max(sizeVals)

#And now for the clonality plots. 
library(ggnetwork)

dir.create("Results/Figure_4_plots/T_cell_networks")
lapply(seq_along(donList), function(x){
    locDon <- donList[[x]]
    adjAll <- matrix(0, nrow(locDon), nrow(locDon))
    for(i in 1:nrow(adjAll)){
        adjAll[i,which(locDon$HamClone == locDon$HamClone[i])] <- 1
    }
    set.seed(1000)
    netAll <- ggnetwork(adjAll)
    netAll$x <- netAll$x*fracVals[x]
    netAll$y <- netAll$y*fracVals[x]
    netAll$Size <- factor(paste0(locDon$Size, "_cells"), levels = paste0(c(1,2,3), "_cells"))
    width1 <- height1 <- 5*sqrt(fracVals[x])
    
    p <- ggplot(netAll, aes(x = x, y = y, size = Size)) + 
        geom_point()  + theme_blank() +coord_equal() + scale_size_discrete(drop = FALSE)
    ggsave(paste0("Results/Figure_4_plots/T_cell_networks/", names(donList)[x], ".pdf"), plot = p, width = width1, height = height1) 
    p <- p + theme(legend.position = "none")
    ggsave(paste0("Results/Figure_4_plots/T_cell_networks/", names(donList)[x], "_no_legend.pdf"), 
           plot = p, width = width1, height = height1) 
})

#A few of the dots fall a bit outside the visual circle of the plots, and these
#have been manually moved in hte publication, as these networks merely are created for visual purposes. 
#Now, we are going to make a simple plot showing the fraction of cell htat are clonal

plotDat <- as.data.frame(table(TCRdbCD4$Sample, TCRdbCD4$Clonal))
colnames(plotDat) <- c("Pat", "Specific", "Freq")
ggplot(plotDat,                         
       aes(x = Pat,
           y = Freq,
           fill = Specific)) + 
    geom_bar(stat = "identity",
             position = "fill", color = "black") + theme_bw() +
    scale_fill_manual(values = c("white", "black"))
ggsave("Results/Figure_4_plots/T_cell_networks/Fraction_clonal.pdf", width = 5, height = 5)

p <- ggplot(locDat,                         
            aes(x = Specific,
                y = value,
                fill = locDat[,1])) + 
    geom_bar(stat = "identity",
             position = "fill") + theme_bw() +
    facet_grid(~ Sample) + scale_fill_viridis(limits = scaleLims, 
                                              discrete = FALSE, 
                                              option = colorOption)

t(table(TCRdbCD4$Sample, TCRdbCD4$Clonal))

apply(as.matrix(table(TCRdbCD4$Sample, TCRdbCD4$Clonal)), 1, function(x) x/sum(x))
#            1166_1     1227_1    1284_1
#FALSE 0.97520661 0.97927461 0.8921283
#TRUE  0.02479339 0.02072539 0.1078717
