library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

#For the visualisations, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]

#In this version, we will collapse the BCRs within the same subclone, as defined
#by the hierarchical tree analysis. 
BCR_reduced <- BCR_H[,c("Sample","HamClone", "Mutations", "light_type", "ISOTYPE",
                        "Specific","SubClone")]
BCR_reduced$n_identical <- 1

#Now, we separate information about the cell type, the isotype and the light type
#into separate columns
for(i in c("light_type", "ISOTYPE")){
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

