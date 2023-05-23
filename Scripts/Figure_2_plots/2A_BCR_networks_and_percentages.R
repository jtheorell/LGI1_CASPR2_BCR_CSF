library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
dir.create("Results/Figure_2_plots")
#For the visualisations, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_all[which(BCR_all$LOCUS == "H"),]

#In this version, we will collapse the BCRs within the same subclone, as defined
#by the hierarchical tree analysis. 
BCR_reduced <- BCR_H[,c("Sample","HamClone", "Specific","SubClone")]
BCR_reduced$n_identical <- 1

#We also add a separate column, where the few cells that do not bind LGI1 or CASPR2
#but that do bind tissue are included separately. 
BCR_reduced$SpecTiss <- BCR_reduced$Specific
BCR_reduced$SpecTiss[which(BCR_H$CELL %in% c("1166_1_2_H12",
                                             "1166_1_2_A15",
                                             "1227_1_2_J08", 
                                             "1227_1_2_O05",
                                             "1166_1_2_J05", 
                                             "1166_1_2_L15",
                                             "1166_1_2_C11", 
                                             "1166_1_2_H08", 
                                             "1166_1_2_N07"))] <- "Tissue"

#Now, we collapse the identical cells
BCR_collapsedClonal <- do.call("rbind",lapply(unique(BCR_reduced$HamClone), function(x){
    locBCR <- BCR_reduced[which(BCR_reduced$HamClone == x),]
    if(length(unique(locBCR$SubClone)) < nrow(locBCR)){
        do.call("rbind",lapply(unique(locBCR$SubClone), function(y){
            identicalRows <- which(locBCR$SubClone == y)
            commonRow <- data.frame("Sample" = locBCR$Sample[1], 
                                    "HamClone" = locBCR$HamClone[1], 
                                    "Specific" = locBCR$Specific[1],
                                    "SpecTiss" = locBCR$SpecTiss[1],
                                    "SubClone" = y,
                                    "n_identical" = length(identicalRows))
            #commonRow <- cbind(commonRow, t(colSums(as.matrix(locBCR[,(1+which(colnames(locBCR) == "n_identical")):ncol(locBCR)]))))
        }))
        
    } else {
        locBCR
    }
}))

#And here we add back all the non-clonal cells
BCR_collapsed <- rbind(BCR_collapsedClonal, BCR_reduced[which(is.na(BCR_reduced$HamClone)),])
#Now over to the plotting. 

dir.create("Results/Figure_2_plots/Network")

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
    netAll2$SpecTiss <- factor(netAll2$SpecTiss, levels = c("TRUE", "FALSE", "Not_tested", "Tissue"))
    
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
    
    commonPlotName <- "Results/Figure_2_plots/Network/Network_"
    width1 <- height1 <- 5*sqrt(x$fracVal[1])

    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = Specific, size = n_identical^2.5))  + scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_manual(values = c("black", "grey", "orange")) + theme_blank() +coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Sample[1], "_specific.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_specific_no_leg.pdf"), height = height1, width = width1)
    
    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = SpecTiss, size = n_identical^2.5))  + 
        scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_manual(values = c("orange", "black", "grey", "yellow")) + theme_blank() +coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Sample[1], "_specTiss.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Sample[1], "_specTiss_no_leg.pdf"), height = height1, width = width1)
    
})

#Here, we create one plot for the percentage of positives for the three individuals, 
#This then being the percentage of all IgG BCRs. 

BCR_IgG_test <- BCR_all[which(grepl("IGHG", BCR_all$ISOTYPE) & BCR_all$Specific != "Not_tested"),]

pdf("Results/Figure_2_plots/Percentage_specific.pdf", width = 5, height = 5)
plot(table(substr(BCR_IgG_test$CELL, 1, 4), BCR_IgG_test$Specific), color = c("black", "orange"),
     main = "")
dev.off()

#Now, for the percentages
locTab <- table(substr(BCR_IgG_test$CELL, 1, 4), BCR_IgG_test$Specific)
locTab/rowSums(locTab)
#     FALSE      TRUE
#1166 0.1836735 0.8163265
#1227 0.2500000 0.7500000
#1284 0.2115385 0.7884615

#Number of tested BCRs: 
rowSums(locTab)
#1166 1227 1284 
#98   16   52 

#We also create a variant of the percentage specific plot, showing the number of cells that are specific to other brain antigens

BCR_spec <- BCR_reduced[-which(BCR_reduced$Specific == "Not_tested"),]
BCR_spec$SpecTiss <- factor(BCR_spec$SpecTiss, levels = c("Tissue", "FALSE", "TRUE"))

pdf("Results/Figure_2_plots/Percentage_specific_plus_tissue.pdf", width = 5, height = 5)
plot(table(BCR_spec$Sample, BCR_spec$SpecTiss), color = c("black", "yellow", "orange"),
     main = "")
dev.off()
#And the percentages are
tissTab <- as.matrix(table(BCR_spec$Sample, BCR_spec$SpecTiss))
apply(tissTab, 1, function(x) x/sum(x))

#           1166_1 1227_1    1284_1
#FALSE  0.11224490  0.125 0.1923077
#Tissue 0.07142857  0.125 0.0000000
#TRUE   0.81632653  0.750 0.8076923
