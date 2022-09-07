#Now, we will in this case zoom in on the IgGs
#Here, we take the file, that has been modified in excel and where antibodies 
#have been selected, back into R and generate figures with it. 
#Strategy-wise, there are two different options. 

#The first aims for maximal
#diversity among the singlets. Here, all subclasses and light chain categories
#are represented if they exist, but not by more than 2 BCRs. This means that
#dominant isotypes, such as IgG1 and IgG4 are vastly underrepresented here, 
#which is a drawback, as they dominate the clones and thus are much more represented
#there, as we on the clonal side, in both strategies aim at finding if each clone 
#binds LGI1/CASPR2 or not and thus get a much more representative distribution, by
#looking at one BCR from each clone.

#Strategy two instead aimsa t representativity, with 25% of all from all groups
#being represented. This means that only BCRs belonging to a group larger than 2
#antibodies will be represented. Here, it is this strategy that has been chosen. 
library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Complete_BCR_db.csv")

BCR_specificity <- read.csv("~/Labbet/2022/220504_LGI1_pos_neg_GE_comparisons/Data/20220426 BCR_selected_db_with exclusions_light chains.csv")
#Step one is to introduce the selection column
BCR_all$Specific <- "Not_tested"
specific <- which(BCR_all$CELL %in% BCR_specificity$CELL[which(BCR_specificity$Antigen_positive == 1)])
BCR_all$Specific[specific] <- "yes"
nonSpecific <- which(BCR_all$CELL %in% BCR_specificity$CELL[which(BCR_specificity$Antigen_positive == 0)])
BCR_all$Specific[nonSpecific] <- "no"
write.csv(BCR_all, "Complete_BCR_db_with_specific.csv")

#In this case, we will zoom in on the IgG BCRs
BCR_IgG <- BCR_all[grep("IGHG", BCR_all$ISOTYPE),]

#And now over to the visualisation. For this, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_IgG[which(BCR_IgG$LOCUS == "H"),]

#We also remove the non-interesting donors. 
BCR_H <- BCR_H[which(substr(BCR_H$CELL, 1, 6) %in% 
                                           c("1166_1", "1227_1", "1284_1")),]

BCR_H$Donor <- substr(BCR_H$CELL, 1, 4)


#In this version, we will collapse the BCRs with no Hamming distance between 
#either the heavy or light chains. This also means that all unique categories are 
#lost, such as cell type. That might be brought back in later. 
BCR_reduced <- BCR_H[,c("Donor","HamClone", "Mutations", "light_type", "ISOTYPE",
                        "Specific","intraClonalDistance", "intraClonalDistanceL",
                        "Cell_type")]
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
BCR_collapsed <- do.call("rbind",sapply(unique(BCR_reduced$HamClone), function(x){
    locBCR <- BCR_reduced[which(BCR_reduced$HamClone == x),]
    if(nrow(locBCR) > 1){
        locBCR$HamClone
        if(any(locBCR$intraClonalDistance == 0 & locBCR$intraClonalDistanceL == 0)){
            identicalRows <- which(locBCR$intraClonalDistance == 0 & locBCR$intraClonalDistanceL == 0)
            commonRow <- data.frame("Donor" = locBCR$Donor[1], "HamClone" = locBCR$HamClone[1], 
                           "Mutations" = mean(locBCR$Mutations), 
                           "Specific" = if(any(locBCR$Specific == "yes")){
                               "yes"
                               } else if(any(locBCR$Specific == "no")){
                                   "no"
                               } else {"Not_tested"},
                           "intraClonalDistance" = 0, 
                           "intraClonalDistanceL" = 0,
                           "n_identical" = length(identicalRows)
                           )
            commonRow <- cbind(commonRow, t(colSums(as.matrix(locBCR[,(1+which(colnames(locBCR) == "n_identical")):ncol(locBCR)]))))
            newData <- rbind(locBCR[-identicalRows,],commonRow)
            newData
        }
        
    } else {
        locBCR
    }
}))

#Now over to the plotting. 

dir.create("Plots_IgG_BCRs")

#Here we bring in a size constant to multiply with, which is dependent on 
#the number of individual BCRs for each donor. 
donorSplit <- split(BCR_collapsed, f = BCR_collapsed$Donor)

sizeVals <- unlist(lapply(donorSplit, nrow))
fracVals <- sizeVals/max(sizeVals)
donorSplit <- lapply(seq_along(donorSplit), function(x){
  locDon <- donorSplit[[x]]
  locDon$fracVal <- fracVals[x]
  locDon
})


#And here they are: all the plots!
dir.create("Results/Graphics/Plots_IgG")
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
    
    netAllExtras <- as.data.frame(matrix(NA, nrow(netAll), ncol(x)-4))
    netAllExtrasData <- x[,-which(colnames(x) %in% c("Donor", "HamClone", 
                                                     "intraClonalDistance",
                                                     "intraClonalDistanceL"))]
    for(i in seq_len(nrow(netAllExtras))){
        locCodes <- netCode[which(netCode == netAll$vertex.names[i])]
        locCodNum <- as.numeric(locCodes[1])
        netAllExtras[i,] <- netAllExtrasData[locCodNum,]
    }
    colnames(netAllExtras) <- colnames(netAllExtrasData)
    
    netAll2 <- cbind(netAll, netAllExtras)
    
    #Here, we introduce a cap for the n_identical, to avoid making certain points so
    #very small. 
    netAll2$n_identical[which(netAll2$n_identical > 3)] <- 3
    
    #To be able to plot the data on a meaningful scale in the pie chart plot,
    #we will have to
    #rescale the x and y variables, so we do this for all. Here, we also introduce
    #a scaling factor that takes the size of the total BCR pool into account. 
    netAll2$x <- sqrt(x$fracVal[1])*100*(netAll2$x-min(netAll2$x))/(max(netAll2$x)-min(netAll2$x))
    netAll2$y <- sqrt(x$fracVal[1])*100*(netAll2$y-min(netAll2$y))/(max(netAll2$y)-min(netAll2$y))
    netAll2$xend <- sqrt(x$fracVal[1])*100*(netAll2$xend-min(netAll2$xend))/(max(netAll2$xend)-min(netAll2$xend))
    netAll2$yend <- sqrt(x$fracVal[1])*100*(netAll2$yend-min(netAll2$yend))/(max(netAll2$yend)-min(netAll2$yend))
    
    
    commonPlotName <- "Results/Graphics/Plots_IgG/Network_"
    width1 <- height1 <- 5*sqrt(x$fracVal[1])
    #First, we create the non-scatterpie plots
    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = Mutations, size = n_identical^2.5))  + scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_viridis(option = "viridis", direction =-1) + theme_blank()+coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Donor[1], "_mutations.pdf"))
    
    p+ theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Donor[1], "_mutations_no_leg.pdf"), height = height1, width = width1)
    
    p <- ggplot(netAll2, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = Specific, size = n_identical^2.5))  + scale_size_continuous(range = c(2.5, 8)) +
        geom_edges() + scale_color_manual(values = c("black", "grey", "orange")) + theme_blank() +coord_equal()
    p
    ggsave(paste0(commonPlotName, x$Donor[1], "_specific.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Donor[1], "_specific_no_leg.pdf"), height = height1, width = width1)
    
    #ANd now for the scatterpies, we need to create color scales, so that the colors
    #are the same, regardless of if all categories are present or not. 
    #And now we do a rather confusing manouvre: we turn the format long again!
    CellTypeData <- reshape2::melt(netAll2[,c("x","y","xend","yend","n_identical", "ASC","Memory","IgDnegCD27neg")],
                                   measure.vars = c("ASC","Memory","IgDnegCD27neg"))
    colnames(CellTypeData)[6] <- "CellType"
    
    p <- ggplot(CellTypeData, aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_scatterpie(aes(x = x, y = y, r = n_identical+0.13), data=CellTypeData,
                        cols= "CellType", long_format = TRUE, color = NA) + geom_edges()+
        coord_equal()+ theme_blank() + 
        scale_fill_manual(values=c("#A65C85", "#DE7065", "#FFCF20"), drop = FALSE)
    p  + geom_scatterpie_legend(netAll2$n_identical, x=1, y=1)
    ggsave(paste0(commonPlotName, x$Donor[1], "_cell_type.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Donor[1], "_cell_type_no_leg.pdf"), 
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
        scale_fill_manual(values=c("red", "blue", "purple"), drop = FALSE)
    p  + geom_scatterpie_legend(netAll2$n_identical, x=1, y=1)
    ggsave(paste0(commonPlotName, x$Donor[1], "_light_type.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Donor[1], "_light_type_no_leg.pdf"), height = height1, width = width1)
    
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
    ggsave(paste0(commonPlotName, x$Donor[1], "_Isotype.pdf"))
    
    p + theme(legend.position = "none") 
    ggsave(paste0(commonPlotName, x$Donor[1], "_Isotype_no_leg.pdf"), height = height1, width = width1)
    
})

#, this data is shown as 

#Now, we are going to construct a file that can be used for visualizations with 
#SPICE. 
#This is a beautiful example of what we want. It will be very easy to get there
#actually, when the dataframe has been designed. 
BCR_reduced2 <-  BCR_H[,c("Donor","HamClone", "Mutations", "light_type", "ISOTYPE",
                          "Specific","intraClonalDistance", "intraClonalDistanceL",
                          "Cell_type")]

#We will now reduce the number of options, by conflating a few of the isotuype options
#(IgE, IgD, unspecified, none)
BCR_reduced2$ISOTYPE[which(BCR_reduced2$ISOTYPE %in% 
                               c("IGHE", "IGHD", "Unknown", "None"))] <- "Other"
    
#Now, to get all levels with all the donors: 

BCR_reduced2$ISOTYPE <- factor(BCR_reduced2$ISOTYPE)
BCR_reduced2$light_type <- factor(BCR_reduced2$light_type)
BCR_reduced2$Clonal <- FALSE
BCR_reduced2$Clonal[-which(is.na(BCR_reduced2$intraClonalDistance))] <- TRUE
BCR_reduced2$MutBinned <- "0"
BCR_reduced2$MutBinned[which(BCR_reduced2$Mutations %in% 1:10)] <- "1-10"
BCR_reduced2$MutBinned[which(BCR_reduced2$Mutations %in% 11:20)] <- "11-20"
BCR_reduced2$MutBinned[which(BCR_reduced2$Mutations > 20)] <- ">20"
BCR_reduced2$MutBinned <- factor(BCR_reduced2$MutBinned, levels = c(">20","11-20", "1-10", "0"))

BCR_reducedSplit <- split(BCR_reduced2, f = BCR_reduced2$Donor)

#Here we new introduce alluvial maps. 
library(ggalluvial)
#T be able to include the mutations in this plot, we need to restructure the data
#into bins. 

dir.create("Results/Graphics/Alluvial_IgG")
commonPlotName <- "Results/Graphics/Alluvial_IgG/Alluvial_"
lapply(names(BCR_reducedSplit), function(x){
    locXDat <- BCR_reducedSplit[[x]]
    locXSize <- nrow(locXDat)
    splitLocXDat <- split(locXDat, locXDat$Specific)
    lapply(names(splitLocXDat), function(y){
        locYDat <- splitLocXDat[[y]]
        locYTab <- as.data.frame(table(
            locYDat$ISOTYPE, locYDat$Clonal, 
            locYDat$MutBinned, locYDat$light_type))
        colnames(locYTab) <- c("Isotype", "Clonal", "Mutations","Light_chain", 
                                   "freq")
        locY_lodes <- to_lodes_form(locYTab,
                                     axes = 1:4,
                                     id = "Aggregate")
        
        #When it comes to the height, we adjust it so that one gets a sense that
        #some groups dominate the others, without this effect taking over completely
        height2 <- 2+(6*nrow(locYDat)/locXSize)

        p <- ggplot(locY_lodes,
               aes(x = x, stratum = stratum, alluvium = Aggregate,
                   y = freq,
                   fill = stratum, label = stratum)) +
            scale_x_discrete(expand = c(.1, .1)) +
            geom_flow() +
            geom_stratum(alpha = .5) +
            theme_bw() + 
            scale_fill_manual(values=c("#B1CFFF", "#9E4532", "#541352", "#238A8D",
                                       "#F58216", "#006400", 
                                       "#E8FA5B", "#F7CB44", "#A65C85", "#042333",
                                       "#A020F0", "#FF0000", "#0000FF")) + 
            scale_y_continuous(expand = c(0, 0))
        p
        ggsave(paste0(commonPlotName, x, "_", y, ".pdf"))
        p + theme(legend.position = "none",
                  axis.text.x=element_blank(),
                  axis.text.y=element_blank(),) +xlab("") + ylab("")
        ggsave(paste0(commonPlotName, x, "_", y, "_no_legend.pdf"), 
               height = height2, width = width1)
    })
})

#Now we run a logistic regression model, more specifically a binominal such. 
binModDat <- BCR_reduced2[-which(BCR_reduced2$Specific == "Not_tested"),]
binModDat$Specific <- sapply(binModDat$Specific, switch, "no" = 0, "yes" = 1)
binModDat$Clonal <- sapply(as.character(binModDat$Clonal), switch, "FALSE" = 0, "TRUE" = 1)
binModDat$light_type <- sapply(as.character(binModDat$light_type), switch, 
                               "double" = 0.5, 
                               "K" = 0,
                               "L" = 1)
binModDat$ISOTYPE <- sapply(as.character(binModDat$ISOTYPE), 
                            switch, "IGHG1" = 0, "IGHG2" = 0.33, 
                            "IGHG3" = 0.66, "IGHG4" = 1)

binMod <- glm(Specific~Clonal+light_type+Mutations+ISOTYPE,
               family=binomial,data=binModDat) 
summary(binMod) #Gives: 
#Call:
#    glm(formula = Specific ~ Clonal + light_type + Mutations + ISOTYPE, 
#        family = binomial, data = binModDat)
#
#Deviance Residuals: 
#    Min       1Q   Median       3Q      Max  
#-2.1897  -1.0244   0.5219   0.8633   1.8232  
#
#Coefficients:
#    Estimate Std. Error z value Pr(>|z|)   
#(Intercept)  1.88535    0.93058   2.026  0.04277 * 
#    Clonal      -0.16397    0.54235  -0.302  0.76240   
#light_type  -1.00427    0.65687  -1.529  0.12629   
#Mutations   -0.10765    0.03936  -2.735  0.00624 **
#    ISOTYPE      1.54958    0.59062   2.624  0.00870 **
#    ---
#    Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#
#(Dispersion parameter for binomial family taken to be 1)
#
#Null deviance: 108.267  on 83  degrees of freedom
#Residual deviance:  91.427  on 79  degrees of freedom
#AIC: 101.43
#
#Number of Fisher Scoring iterations: 4
#