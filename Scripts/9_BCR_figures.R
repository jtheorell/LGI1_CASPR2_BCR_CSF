library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Data/BCR_database_versions/8_surface_pheno_included.csv")

#For the visualisations, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]

#In this version, we will collapse the BCRs within the same subclone, as defined
#by the hierarchical tree analysis. 
BCR_reduced <- BCR_H[,c("Sample","HamClone", "Mutations", "light_type", "ISOTYPE",
                        "Specific","SubClone", "Cell_type")]
BCR_reduced$n_identical <- 1

#Now, we separate information about the cell type, the isotype and the light type
#into separate columns
for(i in c("light_type", "ISOTYPE", "Cell_type")){
    newCols <- unique(BCR_reduced[,i])
    newDf <- as.data.frame(matrix(data = 0, nrow = nrow(BCR_reduced), 
                                  ncol = length(newCols),
                                  dimnames= list(seq(1, nrow(BCR_reduced)),newCols)))
    for(j in newCols){
        newDf[,j][which(BCR_reduced[,i] == j)] <- 1
    }
    BCR_reduced <- cbind(BCR_reduced[,-which(colnames(BCR_reduced) == i)], newDf)
    
}

#Now, we collapse the identical cells
BCR_collapsedClonal <- do.call("rbind",lapply(unique(BCR_reduced$HamClone), function(x){
    locBCR <- BCR_reduced[which(BCR_reduced$HamClone == x),]
    if(length(unique(locBCR$SubClone)) < nrow(locBCR)){
        do.call("rbind",lapply(unique(locBCR$SubClone), function(y){
            identicalRows <- which(locBCR$SubClone == y)
            commonRow <- data.frame("Sample" = locBCR$Sample[1], 
                                    "HamClone" = locBCR$HamClone[1], 
                                    "Mutations" = mean(locBCR$Mutations[identicalRows]), 
                                    "Specific" = locBCR$Specific[1],
                                    "SubClone" = y,
                                    "n_identical" = length(identicalRows)
            )
            commonRow <- cbind(commonRow, t(colSums(as.matrix(locBCR[,(1+which(colnames(locBCR) == "n_identical")):ncol(locBCR)]))))
        }))
        
    } else {
        locBCR
    }
}))

#And here we add back all the non-clonal cells
BCR_collapsed <- rbind(BCR_collapsedClonal, BCR_reduced[which(is.na(BCR_reduced$HamClone)),])

#Now over to the plotting. 

dir.create("Results/Network")

#Here we bring in a size constant to multiply with, which is dependent on 
#the number of individual BCRs for each donor. 
donorSplit <- split(BCR_collapsed, f = BCR_collapsed$Sample)

sizeVals <- unlist(lapply(donorSplit, nrow))
fracVals <- sizeVals/max(sizeVals)
donorSplit <- lapply(seq_along(donorSplit), function(x){
  locDon <- donorSplit[[x]]
  locDon$fracVal <- fracVals[x]
  locDon
})


#And here they are: all the plots!
lapply(donorSplit, function(x){
    adjAll <- matrix(0, nrow(x), nrow(x))
    for(i in 1:nrow(adjAll)){
        adjAll[i,which(x$HamClone == x$HamClone[i])] <- 1
    }
    set.seed(5)
    netAll <- ggnetwork(adjAll)
    #Here, we change the size of the network: 
    netCode <- unique(netAll$vertex.names)[order(unique(netAll$vertex.names))]
    netAllExtras <- matrix(NA, nrow(netAll), 9)
    
    netAllExtras <- as.data.frame(matrix(NA, nrow(netAll), ncol(x)-2))
    netAllExtrasData <- x[,-which(colnames(x) %in% c("Sample", "HamClone"))]
    for(i in seq_len(nrow(netAllExtras))){
        locCodes <- netCode[which(netCode == netAll$vertex.names[i])]
        locCodNum <- as.numeric(locCodes[1])
        netAllExtras[i,] <- netAllExtrasData[locCodNum,]
    }
    colnames(netAllExtras) <- colnames(netAllExtrasData)
    
    netAll2 <- cbind(netAll, netAllExtras)
    
    #Here, we introduce a cap for the n_identical, to avoid making certain points so
    #very small. 
    
    #To be able to plot the data on a meaningful scale in the pie chart plot,
    #we will have to
    #rescale the x and y variables, so we do this for all. Here, we also introduce
    #a scaling factor that takes the size of the total BCR pool into account. 
    netAll2$x <- sqrt(x$fracVal[1])*100*(netAll2$x-min(netAll2$x))/(max(netAll2$x)-min(netAll2$x))
    netAll2$y <- sqrt(x$fracVal[1])*100*(netAll2$y-min(netAll2$y))/(max(netAll2$y)-min(netAll2$y))
    netAll2$xend <- sqrt(x$fracVal[1])*100*(netAll2$xend-min(netAll2$xend))/(max(netAll2$xend)-min(netAll2$xend))
    netAll2$yend <- sqrt(x$fracVal[1])*100*(netAll2$yend-min(netAll2$yend))/(max(netAll2$yend)-min(netAll2$yend))
    
    #Here, we introduce a cap for the n_identical, to avoid making certain points so
    #very small. 
    netAll2$n_identical[which(netAll2$n_identical > 3)] <- 3
    
    commonPlotName <- "Results/Network/Network_"
    width1 <- height1 <- 5*sqrt(x$fracVal[1])
    #First, we create the non-scatterpie plots
    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = Mutations, size = n_identical^2.5))  + scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_viridis(option = "viridis", direction =-1) + theme_blank()+coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Sample[1], "_mutations.pdf"))
    
    p+ theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_mutations_no_leg.pdf"), height = height1, width = width1)
    
    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = Specific, size = n_identical^2.5))  + scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_manual(values = c("black", "grey", "orange")) + theme_blank() +coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Sample[1], "_specific.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_specific_no_leg.pdf"), height = height1, width = width1)
    
    #ANd now for the scatterpies, we need to create color scales, so that the colors
    #are the same, regardless of if all categories are present or not. 
    #And now we do a rather confusing manouvre: we turn the format long again!
    CellTypeData <- reshape2::melt(netAll2[,c("x","y","xend","yend","n_identical", "ASC","B")],
                                   measure.vars = c("ASC","B"))
    colnames(CellTypeData)[6] <- "CellType"
    
    p <- ggplot(CellTypeData, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_scatterpie(aes(x = x, y = y, r = n_identical+0.13), data=CellTypeData,
                        cols= "CellType", long_format = TRUE, color = NA) + geom_edges()+
        coord_equal()+ theme_blank() + 
        scale_fill_manual(values=c("#A65C85", "#DE7065"), drop = FALSE)
    p  + geom_scatterpie_legend(netAll2$n_identical, x=1, y=1)
    ggsave(paste0(commonPlotName, x$Sample[1], "_cell_type.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_cell_type_no_leg.pdf"), 
           height = height1, width = width1)
    
    #######
    LightTypeData <- reshape2::melt(netAll2[,c("x","y","xend","yend","n_identical", "K","L","double")],
                                    measure.vars = c("K","L","double"))
    colnames(LightTypeData)[6] <- "LightType"
    LightTypeData$LightType <- factor(LightTypeData$LightType, levels = c("K","L","double"))
    
    p <- ggplot(LightTypeData, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_scatterpie(aes(x = x, y = y, r = n_identical+0.13), data=LightTypeData,
                        cols= "LightType", long_format = TRUE, color = NA) + geom_edges()+
        coord_equal()+ theme_blank() + 
        scale_fill_manual(values=c("#D40041",  "#2E00B7", "#8F00B9"), drop = FALSE)
    p  + geom_scatterpie_legend(netAll2$n_identical, x=1, y=1)
    ggsave(paste0(commonPlotName, x$Sample[1], "_light_type.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_light_type_no_leg.pdf"), height = height1, width = width1)
    
    #######
    IsoTypeData <- reshape2::melt(netAll2[,c("x","y","xend","yend","n_identical", 
                                             "IGHG1", "IGHG2", "IGHG3", 
                                             "IGHG4")],
                                  measure.vars = c("IGHG1", "IGHG2", "IGHG3", 
                                                   "IGHG4"))
    colnames(IsoTypeData)[6] <- "IsoType"
    IsoTypeData$IsoType <- factor(IsoTypeData$IsoType, 
                                  levels = c("IGHG1", "IGHG2", "IGHG3", 
                                             "IGHG4"))
    
    p <- ggplot(IsoTypeData, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_scatterpie(aes(x = x, y = y, r = n_identical+0.13), data=IsoTypeData,
                        cols= "IsoType", long_format = TRUE, color = NA) + geom_edges()+
        coord_equal()+ theme_blank() + 
        scale_fill_manual(values=c("#B1CFFF", "#9E4532", "#541352", "#238A8D"), drop = FALSE)
    p  + geom_scatterpie_legend(netAll2$n_identical, x=1, y=1)
    ggsave(paste0(commonPlotName, x$Sample[1], "_Isotype.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_Isotype_no_leg.pdf"), height = height1, width = width1)
    
})

#Now, we will create simple stacked bar graphs for the above parameters. 
barGraphDat <- BCR_H[,c("Sample","Clonal", "Mutations", "light_type", "ISOTYPE",
                        "Specific", "Cell_type")]
barGraphDat$light_type <- factor(barGraphDat$light_type, levels = c("K", "double", "L"))
barGraphDat$ISOTYPE <- factor(barGraphDat$ISOTYPE)

barGraphDat$MutBinned <- "0"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 1:5)] <- "1-5"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 6:10)] <- "6-10"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 11:15)] <- "11-15"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 16:20)] <- "16-20"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 21:25)] <- "21-25"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 26:30)] <- "26-30"
barGraphDat$MutBinned[which(barGraphDat$Mutations %in% 31:35)] <- "31-35"
barGraphDat$MutBinned[which(barGraphDat$Mutations > 35)] <- ">35"
barGraphDat$MutBinned <- factor(barGraphDat$MutBinned, levels = c(">35", "31-35",
                                                                  "26-30", "21-25",
                                                                  "16-20", "11-15",
                                                                  "6-10", "1-5", "0"))

colsAndNames <- list(list("Clonal", c("grey","black")),
                     list("MutBinned", rev(viridis(9))),
                     list("light_type", c("#D40041", "#8F00B9",  "#2E00B7")),
                     list("ISOTYPE", c("#B1CFFF", "#9E4532", "#541352", "#238A8D")),
                     list("Cell_type", c("#A65C85", "#DE7065")))
dir.create("Results/CloneBarGraphs")
for(i in colsAndNames){
    locDat <- as.data.frame(table(barGraphDat[,i[[1]]], barGraphDat$Sample, barGraphDat$Specific))
    colnames(locDat) <- c(i[[1]], "Sample","Specific", "value")
    p <- ggplot(locDat,                         
           aes(x = Specific,
               y = value,
               fill = locDat[,1])) + 
        geom_bar(stat = "identity",
                 position = "fill") + theme_bw() +
        facet_grid(~ Sample) + scale_fill_manual(values=i[[2]], drop = FALSE)
    p
    ggsave(paste0("Results/CloneBarGraphs/", i[[1]], ".pdf"))
    p + theme_blank() + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 5, height = 6)
    
}

#And now, we display the data in another way, namely by stacked pie charts.

#After many trials and tribulations I have realized that hte best way forward
#is to produce three independent graphs and stack them manually in illustrator. 
#the reason for this is that the colors otherwise become very hard to unerstand. 
pieDat <- barGraphDat[-which(barGraphDat$Specific == "Not_tested"),]
pieDatSplit <- split(pieDat, paste0(pieDat$Sample, "_", pieDat$Specific))
#We also add all donors as a fourth category
pieDatSplit$All_TRUE <- pieDat[which(pieDat$Specific == "TRUE"),]
pieDatSplit$All_FALSE <- pieDat[which(pieDat$Specific == "FALSE"),]

dir.create("Results/Raw_pies")
lapply(names(pieDatSplit), function(x){
    locDf <- pieDatSplit[[x]]
    locDfSplit <- split(locDf, locDf$Clonal)
    lapply(names(locDfSplit), function(y){
        locInnerDf <- locDfSplit[[y]]
        ####
        firstTable <- as.data.frame(table(locInnerDf$ISOTYPE))
        p <- ggplot(firstTable, aes(y=Freq)) +
            geom_bar(aes(fill=factor(Var1), x=0), width=.5, 
                     stat='identity', color = "black") +
            coord_polar(theta='y') + theme_minimal() +
            scale_fill_manual(values = c("#B1CFFF", "#9E4532", "#541352", "#238A8D"))
        p
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_first.pdf"))
        p + theme_void() + theme(legend.position = "none")
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_first_no_legend.pdf"))
        ####
        secondTable <- as.data.frame(table(locInnerDf$light_type,
                                           locInnerDf$ISOTYPE))
        secondTable$CombVar <- paste0(secondTable$Var1, "_",
                                      secondTable$Var2)
        secondTable$CombVar <- factor(secondTable$CombVar,
                                      levels = secondTable$CombVar)
        p <- ggplot(secondTable, aes(y=Freq)) +
            geom_bar(aes(fill=factor(CombVar), x=0), width=.5, 
                     stat='identity', color = "black") +
            coord_polar(theta='y') + theme_minimal() +
            scale_fill_manual(values = rep(c("#D40041", "#8F00B9",  "#2E00B7"),
                              length.out= 12))
        p
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_second.pdf"))
        p + theme_void() + theme(legend.position = "none")
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_second_no_legend.pdf"))
        ####
        thirdTable <- as.data.frame(table(locInnerDf$MutBinned,
                                          locInnerDf$light_type,
                                          locInnerDf$ISOTYPE))
        thirdTable$CombVar <- paste0(thirdTable$Var1, "_",
                                     thirdTable$Var2, "_",
                                     thirdTable$Var3)
        thirdTable$CombVar <- factor(thirdTable$CombVar,
                                     levels = thirdTable$CombVar)
        p <- ggplot(thirdTable, aes(y=Freq)) +
            geom_bar(aes(fill=factor(CombVar), x=0), width=.5, 
                     stat='identity', color = "black") +
            coord_polar(theta='y') + theme_minimal() +
            scale_fill_manual(values = rep(rev(viridis(9)), 
                                           length.out= 108))
        p
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_third.pdf"))
        p + theme_void() + theme(legend.position = "none")
        ggsave(paste0("Results/Raw_pies/", x, "_", y, "_third_no_legend.pdf"))
    })
})
