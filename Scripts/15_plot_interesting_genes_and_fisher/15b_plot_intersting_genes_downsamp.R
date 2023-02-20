library(SingleCellExperiment)
library(scater)
#Here, we are going to produce individual gene plots for a few marker genes as 
#well as the three top genes in the differential expression analysis
#After much going back and forth, we will here show DPEP1 and TCN2
focusGenes <- c("DPEP1", "TCN2", "CD19", "MS4A1", "XBP1", "SDC1", "IRF4", "PRDM1", "CD79A")
# PRDM1 = BLIMP1
# MS4A1 = CD20
sceList <- readRDS("Data/Comp_to_others/Downsampled_B_and_LGI1-C2_neg_included.rds")

sceListRed <- lapply(sceList, function(x) x[which(row.names(x) %in% focusGenes),])

plotDatList <- lapply(names(sceListRed), function(x){
    locDat <- sceListRed[[which(names(sceListRed) == x)]]
    list(logcounts(locDat), "colDat" = DataFrame("dataset" = rep(x, ncol(locDat)),
                                                 "cellType" = locDat$cellType))
})

plotSce <- SingleCellExperiment(assays = list("logcounts" = do.call("cbind", lapply(plotDatList, "[[", 1))),
                                              colData = do.call("rbind", lapply(plotDatList, "[[", 2)))

plotSce$dataCellType <- paste0(plotSce$dataset, ".", plotSce$cellType)

plotSce$dataCellType <- factor(plotSce$dataCellType, levels = c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                                                "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                                                "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                                                "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B",
                                                                "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                                                "Covid_PBMC.ASC", "Covid_PBMC.B",
                                                                "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                                                "Influensa_PBMC.ASC", "Influensa_PBMC.B"))

dir.create("Results/Gene_expression_plots/Downsamp")

for(i in focusGenes){
    p <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + ylim(0, 8)
    ggsave(paste0("Results/Gene_expression_plots/Downsamp/", i, ".pdf"), plot = p, height = 7, width = 8)
    
    q <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + 
        theme_void() +  ylim(0, 8) + 
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              strip.text.x = element_blank(),
              legend.position="none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
    ggsave(paste0("Results/Gene_expression_plots/Downsamp/", i, ".png"), plot = q, height = 7, width = 8)
}

#We also create a second one for DPEP1, as the scale looks so silly. 
i <- "DPEP1"
p <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave(paste0("Results/Gene_expression_plots/Downsamp/", i, ".pdf"), plot = p, height = 7, width = 8)

q <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + 
    theme_void() + 
    theme(axis.text.x = element_blank(), 
          axis.ticks.x = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(), 
          strip.text.x = element_blank(),
          legend.position="none",
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave(paste0("Results/Gene_expression_plots/Downsamp/", i, ".png"), plot = q, height = 7, width = 8)

#Now, looking at these plots, it is worth asking the question: is there even a statistical
#difference between our specific cells and all the other datasets combined, when it comes
#to the likelihood of expression? We check this in a crude, simple way, namely by Fisher. 
#Of course, this shows the down-sampled data, but in the raw data, we only have 15 cells positive. 

aeSpecDat <- plotSce[,which(plotSce$dataset == "AE_pos")]
specAll <- ncol(aeSpecDat)
specExpr <- length(which(logcounts(aeSpecDat)[which(row.names(aeSpecDat) == "DPEP1"),] > 0))
specNonExpr <- specAll-specExpr
otherDat <- plotSce[,-which(plotSce$dataset == "AE_pos")]
nonSpecAll <- ncol(otherDat)
nonSpecExpr <- length(which(logcounts(otherDat)[which(row.names(otherDat) == "DPEP1"),] > 0))
nonSpecNonExpr <- nonSpecAll-nonSpecExpr

fisherDf <- data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
                       "Spec_ASC" = c(specNonExpr, specExpr))
row.names(fisherDf) <- c("Not_expressed", "Expressed")

fisherDf
#              Non_spec Spec_ASC
#Not_expressed    41499      117
#Expressed           35        6

fisher.test(fisherDf, alternative = "greater")
#p-value = 2.421e-09, so highly significant. When only the raw data for our dataset
#is compared, we get a p-value of 0.04. What about the other one? 

aeSpecDat <- plotSce[,which(plotSce$dataset == "AE_pos")]
specAll <- ncol(aeSpecDat)
specExpr <- length(which(logcounts(aeSpecDat)[which(row.names(aeSpecDat) == "TCN2"),] > 0))
specNonExpr <- specAll-specExpr
otherDat <- plotSce[,-which(plotSce$dataset == "AE_pos")]
nonSpecAll <- ncol(otherDat)
nonSpecExpr <- length(which(logcounts(otherDat)[which(row.names(otherDat) == "TCN2"),] > 0))
nonSpecNonExpr <- nonSpecAll-nonSpecExpr

fisherDf <- data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
                       "Spec_ASC" = c(specNonExpr, specExpr))
row.names(fisherDf) <- c("Not_expressed", "Expressed")

fisherDf
#              Non_spec Spec_ASC
#Not_expressed    41486      119
#Expressed           48        4

fisher.test(fisherDf, alternative = "greater")
#p-value = 1.756e-05, so also here it is most significant. Worth checking too of course
#if it holds true even if the B-cells are removed. 
plotSceAsc <- plotSce[,which(plotSce$cellType == "ASC")]

aeSpecDat <- plotSceAsc[,which(plotSceAsc$dataset == "AE_pos")]
specAll <- ncol(aeSpecDat)
specExpr <- length(which(logcounts(aeSpecDat)[which(row.names(aeSpecDat) == "DPEP1"),] > 0))
specNonExpr <- specAll-specExpr
otherDat <- plotSceAsc[,-which(plotSceAsc$dataset == "AE_pos")]
nonSpecAll <- ncol(otherDat)
nonSpecExpr <- length(which(logcounts(otherDat)[which(row.names(otherDat) == "DPEP1"),] > 0))
nonSpecNonExpr <- nonSpecAll-nonSpecExpr

fisherDf <- data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
                       "Spec_ASC" = c(specNonExpr, specExpr))
row.names(fisherDf) <- c("Not_expressed", "Expressed")

fisherDf
#              Non_spec Spec_ASC
#Not_expressed     1474      113
#Expressed           34        6

fisher.test(fisherDf, alternative = "greater")
#p-value = 0.06659, so a border land. 

plotSceAsc <- plotSce[,which(plotSce$cellType == "ASC")]

aeSpecDat <- plotSceAsc[,which(plotSceAsc$dataset == "AE_pos")]
specAll <- ncol(aeSpecDat)
specExpr <- length(which(logcounts(aeSpecDat)[which(row.names(aeSpecDat) == "TCN2"),] > 0))
specNonExpr <- specAll-specExpr
otherDat <- plotSceAsc[,-which(plotSceAsc$dataset == "AE_pos")]
nonSpecAll <- ncol(otherDat)
nonSpecExpr <- length(which(logcounts(otherDat)[which(row.names(otherDat) == "TCN2"),] > 0))
nonSpecNonExpr <- nonSpecAll-nonSpecExpr

fisherDf <- data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
                       "Spec_ASC" = c(specNonExpr, specExpr))
row.names(fisherDf) <- c("Not_expressed", "Expressed")

fisherDf
#              Non_spec Spec_ASC
#Not_expressed     1486      115
#Expressed           22        4
fisher.test(fisherDf, alternative = "greater")
#p-value = 0.1166
