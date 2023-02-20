library(DepecheR)
library(scater)
#Here, we are going to make a plot showing the B-cell compartment in all its
#complexity and overlay the specificity information. 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

#We start by creating one plot for the percentage of positives for the three individuals, 
#This then being the percentage of all IgG BCRs. 

BCR_IgG_test <- BCR_all[which(grepl("IGHG", BCR_all$ISOTYPE) & BCR_all$Specific != "Not_tested"),]

pdf("Results/Percentage_specific.pdf", width = 5, height = 5)
plot(table(substr(BCR_IgG_test$CELL, 1, 4), BCR_IgG_test$Specific), color = c("black", "orange"),
     main = "")
dev.off()

#Now, for the percentages
locTab <- table(substr(BCR_IgG_test$CELL, 1, 4), BCR_IgG_test$Specific)
locTab/rowSums(locTab)
#     FALSE      TRUE
#1166 0.1836735 0.8163265
#1227 0.1250000 0.8750000
#1284 0.2115385 0.7884615

#Number of tested BCRs: 
rowSums(locTab)
#1166 1227 1284 
#98   16   52 

#And here the number of tested, unrelated BCRs
specIndiv <- BCR_IgG_test[which(BCR_IgG_test$Clonal == FALSE | 
                                    duplicated(BCR_IgG_test$HamClone) == FALSE),]
table(substr(specIndiv$CELL,1,4))
#1166 1227 1284 
#53   11   21 

#Now, over to the ASC question. 

#Now, we will generate an UMAP specifically to separate the ASC from the 
#B cells 


dir.create("Results/Flow_specific")
flowDatLoc <- as.data.frame(t(normcounts(altExp(csfSce, "flowData")))[,c("CD38", "CD138", "CD27", "IgD", "CD20")])
#As CD138 has an important role in separation of these cell types we will give it
#double weight
flowDatLoc$CD138 <- flowDatLoc$CD138*2
set.seed(111)
bUmap <- uwot::umap(flowDatLoc)
dColorPlot(flowDatLoc, xYData = bUmap, plotDir = "Results/Flow_specific/Markers")

#To make this plot look nice, I will reorder the events to show up with the 
#specifics on top. Importantly, a large group of cells have not been considered for this
#at all, as they are non-IgG. They are included as "not_tested" here. 
csfSce$Specific[which(is.na(csfSce$Specific))] <- "Not_tested"
flowDatLoc$Specific <- factor(csfSce$Specific, levels = c("TRUE", "FALSE", "Not_tested"))
bFlowReordered <- flowDatLoc[order(flowDatLoc$Specific, decreasing = TRUE),]
umapReordered <- bUmap[order(flowDatLoc$Specific, decreasing = TRUE),]

dColorPlot(as.character(bFlowReordered$Specific), colorScale = c("black", "grey", "orange"),
           xYData = umapReordered, plotName = "Results/Flow_specific/Specific", 
           densContour = FALSE)

#Now, we will make one further plot set, namely showing the five subsets as colors
#and shapes on the UMAP. 
plot(flowDatLoc$IgD, flowDatLoc$CD27)
plot(flowDatLoc$IgD[which(csfSce$cellType == "B")],
     flowDatLoc$CD27[which(csfSce$cellType == "B")])

#When investigating these two plots, it is clear that the bulk of the CD27 positive
#cells are ASC. 
flowDatLoc$fineGrainedCellType <- "ASC"
flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                       flowDatLoc$CD27 > 2 &
                                       flowDatLoc$IgD > 1)] <- "Unswitched_mem"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                       flowDatLoc$CD27 > 2 &
                                       flowDatLoc$IgD <= 1)] <- "Mem"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                       flowDatLoc$CD27 <= 2 &
                                       flowDatLoc$IgD > 1)] <- "Naïve"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                       flowDatLoc$CD27 <= 2 &
                                       flowDatLoc$IgD <= 1)] <- "Double_neg"

dColorPlot(flowDatLoc$fineGrainedCellType, 
           xYData = bUmap, plotName = "Results/Flow_specific/Cell_type", 
           colorScale = "dark_rainbow", densContour = FALSE)

#NOw, we are going to make bar graphs with this informaiton. 

igGDat <- flowDatLoc[which(colnames(csfSce) %in% BCR_IgG_test$CELL),]
pdf("Results/Specificity_vs_cell_type.pdf", width = 4, height = 6)
plot(table(as.character(igGDat$Specific), igGDat$fineGrainedCellType),
     color = dColorVector(1:5, colorScale = "dark_rainbow"), main = "")
dev.off()
length(which(igGDat$Specific == "TRUE" & igGDat$fineGrainedCellType == "ASC"))
#117
length(which(igGDat$Specific == "TRUE" & igGDat$fineGrainedCellType != "ASC"))
#5
#So a tiny fraction in the specific compatment are non-ASC. 

length(which(igGDat$Specific == "FALSE" & igGDat$fineGrainedCellType == "ASC"))
#18
length(which(igGDat$Specific == "FALSE" & igGDat$fineGrainedCellType != "ASC"))
#12


#Now, we also plot a few ASC genes
#Here, we are going to plot a few selected genes to show that the genes we 
#have selected are indeed ASC. 
csfSce$b_cell_type <- flowDatLoc$fineGrainedCellType
plotSce <- csfSce[which(rowData(csfSce)$hgnc_symbol %in% c("XBP1", "MS4A1")),]

#PRDM1 is the new name for BLIMP1, MS4A1 is CD20
rownames(plotSce) <- rowData(plotSce)$hgnc_symbol

g <- plotExpression(plotSce, features=rownames(plotSce), 
               x="b_cell_type", colour_by="b_cell_type") + theme_bw() + 
    scale_color_manual(values = dColorVector(1:5, colorScale = "dark_rainbow"))
ggsave("Results/ASC_vs_B.pdf", width = 6, height = 4)
g + theme_bw() +theme(axis.text.x=element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none"
) +xlab("") + ylab("")
ggsave("Results/ASC_vs_B_no_leg.pdf", width = 4, height = 4)
