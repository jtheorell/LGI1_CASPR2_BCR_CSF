library(DepecheR)
library(scater)
#Here, we are going to make a plot showing the B-cell compartment in all its
#complexity and overlay the specificity information. 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")


#Now, we will generate an UMAP specifically to separate the ASC from the 
#B cells 

dir.create("Results/Figure_2_plots/Flow_specific", recursive = TRUE)
flowDatLoc <- as.data.frame(t(normcounts(altExp(csfSce, "flowData")))[,c("CD38", "CD138", "CD27", "IgD", "CD20")])
#As CD138 has an important role in separation of these cell types we will give it
#extra weight
flowDatLoc$CD138 <- flowDatLoc$CD138*1.5
set.seed(11)
bUmap <- uwot::umap(flowDatLoc)
dColorPlot(flowDatLoc, xYData = bUmap, plotDir = "Results/Figure_2_plots/Flow_specific/Markers")

#To make this plot look nice, I will reorder the events to show up with the 
#specifics on top. Importantly, a large group of cells have not been considered for this
#at all, as they are non-IgG. They are included as "not_tested" here. 
csfSce$Specific[which(is.na(csfSce$Specific))] <- "Not_tested"
flowDatLoc$Specific <- factor(csfSce$Specific, levels = c("FALSE", "TRUE", "Not_tested"))
bFlowReordered <- flowDatLoc[order(flowDatLoc$Specific, decreasing = TRUE),]
umapReordered <- bUmap[order(flowDatLoc$Specific, decreasing = TRUE),]

dColorPlot(as.character(bFlowReordered$Specific), colorScale = c("black", "grey", "orange"),
           xYData = umapReordered, plotName = "Results/Figure_2_plots/Flow_specific/Specific", 
           densContour = FALSE)
dColorPlot(as.character(bFlowReordered$Specific)[-which(bFlowReordered$Specific == "Not_tested")], 
           colorScale = c("black", "grey", "orange"),
           xYData = umapReordered[-which(bFlowReordered$Specific == "Not_tested"),], 
           dotSize = 20,
           plotName = "Results/Figure_2_plots/Flow_specific/Specific_spec_only", 
           densContour = FALSE)


#Now, we will make one further plot set, namely showing the five subsets as colors
#and shapes on the UMAP. 
plot(flowDatLoc$IgD, flowDatLoc$CD27)
plot(flowDatLoc$IgD[which(csfSce$cellType == "B")],
     flowDatLoc$CD27[which(csfSce$cellType == "B")])

#Sadly, the protein CD27 expression is much poorer than the RNAseq
plot(logcounts(csfSce)[which(rowData(csfSce)$hgnc_symbol == "CD27"),
                       which(csfSce$cellType == "B")],
     flowDatLoc$CD27[which(csfSce$cellType == "B")])
abline(h = 2.3, v = 1)

CD27RNA <- logcounts(csfSce)[which(rowData(csfSce)$hgnc_symbol == "CD27"),]
#When investigating these two plots, it is clear that the bulk of the CD27 positive
#cells are ASC. 
flowDatLoc$fineGrainedCellType <- "ASC"
flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                       (flowDatLoc$CD27 > 2.3 | CD27RNA > 0) &
                                       flowDatLoc$IgD > 1)] <- "Unswitched_mem"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                         (flowDatLoc$CD27 > 2.3 | CD27RNA > 0) &
                                       flowDatLoc$IgD <= 1)] <- "Mem"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                         (flowDatLoc$CD27 <= 2.3 & CD27RNA <= 0) &
                                       flowDatLoc$IgD > 1)] <- "Naïve"

flowDatLoc$fineGrainedCellType[which(csfSce$cellType == "B" & 
                                         (flowDatLoc$CD27 <= 2.3 & CD27RNA <= 0) &
                                       flowDatLoc$IgD <= 1)] <- "Double_neg"

dColorPlot(flowDatLoc$fineGrainedCellType, 
           xYData = bUmap, plotName = "Results/Figure_2_plots/Flow_specific/Fine_grained_cell_type", 
           colorScale = c("#00204DFF","#D7C463FF", "#838079FF","#FFEA46FF", "#ABA074FF"), 
           densContour = FALSE)

#And we make a version of this with the non-fine grained cell type
flowDatLoc$cellType <- flowDatLoc$fineGrainedCellType
flowDatLoc$cellType[-which(flowDatLoc$cellType =="ASC")] <- "B"
dColorPlot(flowDatLoc$cellType, 
           xYData = bUmap, plotName = "Results/Figure_2_plots/Flow_specific/Cell_type", 
           colorScale = c("#00204DFF","red"), 
           densContour = FALSE)

dColorPlot(flowDatLoc$cellType[-which(flowDatLoc$Specific == "Not_tested")], 
           xYData = bUmap[-which(flowDatLoc$Specific == "Not_tested"),], plotName = "Results/Figure_2_plots/Flow_specific/Cell_type_spec", 
           colorScale = c("#00204DFF","red"), 
           dotSize = 20,
           densContour = FALSE)

#NOw, we are going to make bar graphs with this informaiton. 

igGDat <- flowDatLoc[-which(flowDatLoc$Specific == "Not_tested"),]
#pdf("Results/Figure_2_plots/Specificity_vs_cell_type.pdf", width = 4, height = 6)
#plot(table(as.character(igGDat$Specific), igGDat$fineGrainedCellType),
#     color = dColorVector(1:5, colorScale = "dark_rainbow"), main = "")
#dev.off()

spec_cell_df <- as.data.frame(table(as.character(igGDat$Specific), igGDat$fineGrainedCellType))
colnames(spec_cell_df) <- c("Specific", "B_subtype", "Freq")
spec_cell_df$B_subtype <- factor(spec_cell_df$B_subtype, c("ASC", "Mem", "Unswitched_mem", "Double_neg", "Naïve"))
ggplot(spec_cell_df, aes(x=Specific, y=Freq, fill=B_subtype)) +
    geom_bar(stat="identity") + theme_bw() + 
    scale_y_continuous(expand = expansion(mult=c(0,0.05))) +
    scale_fill_manual(values = c("#00204DFF","#838079FF", "#ABA074FF", "#D7C463FF", "#FFEA46FF"))
ggsave("Results/Figure_2_plots/Specificity_vs_fine-grained_cell_type.pdf", width = 5, height = 6)

#Now, the same without the finer details. 
spec_cell_df <- as.data.frame(table(as.character(igGDat$Specific), igGDat$cellType))
colnames(spec_cell_df) <- c("Specific", "Cell_type", "Freq")
spec_cell_df$B_subtype <- factor(spec_cell_df$Cell_type, c("ASC", "B"))
ggplot(spec_cell_df, aes(x=Specific, y=Freq, fill=Cell_type)) +
    geom_bar(stat="identity") + theme_bw() + 
    scale_y_continuous(limits = c(0,125),
                       expand = expansion(mult=c(0,0))) +
    scale_fill_manual(values = c("#00204DFF","red"))
ggsave("Results/Figure_2_plots/Specificity_vs_cell_type.pdf", width = 5, height = 6)



length(which(igGDat$Specific == "TRUE" & igGDat$fineGrainedCellType == "ASC"))
#116
length(which(igGDat$Specific == "TRUE" & igGDat$fineGrainedCellType != "ASC"))
#5
#So a tiny fraction in the specific compatment are non-ASC. 

length(which(igGDat$Specific == "FALSE" & igGDat$fineGrainedCellType == "ASC"))
#19
length(which(igGDat$Specific == "FALSE" & igGDat$fineGrainedCellType != "ASC"))
#12


#We will also make a fisher test here. 
fisherDf <- data.frame("Specific" = c(114,4),
                       "Not_specific" = c(17,10))
row.names(fisherDf) <- c("ASC", "B")

fisherDf
#    Specific Not_specific
#ASC      114           17
#B          4           10


fisher.test(fisherDf)
#p-value = 6.271e-06, i.e. the B-cell frequency is significantly lower in the specific group. 

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
ggsave("Results/Figure_2_plots/ASC_vs_B.pdf", width = 6, height = 4)
g + theme_bw() +theme(axis.text.x=element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.text.y=element_blank(),
                      axis.ticks.y = element_blank(),
                      legend.position = "none"
) +xlab("") + ylab("")
ggsave("Results/Figure_2_plots/ASC_vs_B_no_leg.pdf", width = 4, height = 4)

#And we save the fine grained cell types separately, as they are not used much downstream
#and it is unnecessary to save a completely new sce for this purpose. 
flowDatLoc$Cell <- colnames(csfSce)
saveRDS(flowDatLoc, "Data/Cytometry/Flow_data_and_cell_type_post_integration.rds")

