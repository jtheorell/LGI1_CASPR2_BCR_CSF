#Here, we are going to investigate those few but intriguing binder B-cells. 

fullBCRDb <- read.csv("Data/BCR_database_versions/4_IMGT_gapped_db_complete_post_clonality.csv")
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#Now, we import information also about the non-IgG, as we are now digging deeper
#into the B-cell compartment. 
fullBCRDbH <- fullBCRDb[which(fullBCRDb$LOCUS == "H" & fullBCRDb$CELL %in% colnames(aeSce)),
                        c("CELL", "V_CALL","JUNCTION",
                          "JUNCTION_LENGTH","ISOTYPE",
                          "HamClone","light_type",
                          "Clonal","All_mutations",
                          "Non_silent_mutations")]

identical(colnames(colData(aeSce))[which(colnames(colData(aeSce)) %in% c("V_CALL","JUNCTION",
                                                                         "JUNCTION_LENGTH",
                                                                         "ISOTYPE","HamClone",
                                                                         "light_type","Clonal",
                                                                         "All_mutations",
                                                                         "Non_silent_mutations"))], 
          colnames(fullBCRDbH)[-1])
#TRUE, so we will be able to transfer this in blocks
missingDf <- as.data.frame(matrix(NA, ncol(aeSce)-nrow(fullBCRDbH), 
                     ncol(fullBCRDbH)))
colnames(missingDf) <- colnames(fullBCRDbH)
missingDf$CELL <- colnames(aeSce)[-which(colnames(aeSce) %in% fullBCRDbH$CELL)]

BCRPlusEmpty <- rbind(fullBCRDbH, missingDf)
BCR_ordered <- BCRPlusEmpty[match(colnames(aeSce), BCRPlusEmpty$CELL),]

identical(BCR_ordered$CELL, colnames(aeSce))
#TRUE
BCR_ordered$CELL <- NULL

colData(aeSce)[,c("V_CALL","JUNCTION",
                  "JUNCTION_LENGTH",
                  "ISOTYPE","HamClone",
                  "light_type","Clonal",
                  "All_mutations",
                  "Non_silent_mutations")] <- BCR_ordered

#Now start the real analyses

aeSceB <- aeSce[,which(aeSce$cellType == "B")]


#Now, we are going to pick out the ones that 
aeSceBSpecKnown <- aeSceB[,which(aeSceB$Specific != "Not_tested")]

table(paste0(aeSceBSpecKnown$Specific, "_", aeSceBSpecKnown$donor), 
      aeSceBSpecKnown$Clonal, useNA = "ifany")
#                 FALSE TRUE
#FALSE_1166          3    1
#FALSE_1227          2    0
#FALSE_1284          6    0
#TRUE_1166           0    2
#TRUE_1284           1    2

#So much overrepresented is clonality and the cells here are luckily distributed among
#two donors and two disorders. 

#All in all, these cells show ASC traits, despite not having an ASC surface phenotype
library(ggplot2)

plotDat <-  data.frame("Specific" = aeSce$Specific, "cellType" = aeSce$cellType,
                       "CD27RNA" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "CD27"),],
                       t(normcounts(altExp(aeSce, "flowData"))),
                       "XBP1" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "XBP1"),],
                       "CD138RNA" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "SDC1"),],
                       "TNFRSF17RNA" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "TNFRSF17"),],
                       "CTCF" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "CTCF"),],
                       "IRF4" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "IRF4"),],
                       "BLIMP1" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "PRDM1"),],
                       "KI67" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "MKI67"),],
                       "HLADRA" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "HLA-DRA"),],
                       "CD19" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "CD19"),],
                       "CD79A" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "CD79A"),])

#Now, we remove the non-specific, as they only make noise at this stage. 
plotDat <- plotDat[-which(plotDat$Specific == "Not_tested"),]

#Now, we want the specifics to be on top
plotDat$Specific <- factor(plotDat$Specific, levels = c("TRUE", "FALSE"))

plotDatOrdered <- plotDat[order(plotDat$Specific),]

#We also introduce a size parameter, to make the specific B-cells come out
#clearer. 
plotDatOrdered$specB <- FALSE
plotDatOrdered$specB[which(plotDatOrdered$Specific == "TRUE" & plotDatOrdered$cellType == "B")] <-  TRUE

dir.create("Results/Figure_3_plots/Specific_B_cell_analysis")

ggplot(plotDatOrdered, aes(x = CD27RNA, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/CD27RNA_vs_CD20.pdf", width = 5, height = 6)

ggplot(plotDatOrdered, aes(x = CD79A, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/CD79A_vs_CD20.pdf", width = 5, height = 6)

ggplot(plotDatOrdered, aes(x = XBP1, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/XBP1_vs_CD20.pdf", width = 5, height = 6)

ggplot(plotDatOrdered, aes(x = IRF4, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/IRF4_vs_CD20.pdf", width = 5, height = 6)

ggplot(plotDatOrdered, aes(x = BLIMP1, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/BLIMP1_vs_CD20.pdf", width = 5, height = 6)

ggplot(plotDatOrdered, aes(x = HLADRA, y = CD20, 
                           color = Specific, 
                           shape = cellType,
                           size = specB)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("orange", "black")) +
    scale_shape_manual(values = c(17,15)) + 
    scale_size_manual(values = c(4, 7))
ggsave("Results/Figure_3_plots/Specific_B_cell_analysis/HLA-DRA_vs_CD20.pdf", width = 5, height = 6)

#ANd now some stats. We check if the expression of these markers
#is increased in specific B compared to non-specific (FIsher is 
#the other way around, hence "less")
bDatOrd <- plotDatOrdered[which(plotDatOrdered$cellType == "B"),]
fisherList <- lapply(c("XBP1", "BLIMP1", "CD27RNA", "IRF4"), function(x){
    locDat <- bDatOrd[,which(colnames(bDatOrd) %in% c(x, "Specific"))]
    locDat$expr <- "Neg"
    locDat$expr[which(locDat[,x] > 0)] <- "Pos"
    fisher.test(table(locDat$expr, locDat$Specific), alternative = "less")$p.value
})
names(fisherList) <- c("XBP1", "BLIMP1", "CD27RNA", "IRF4")

#Now, we are going to pick out the ones that are prePB: 
pcSce <- readRDS("Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")
prePBSce <- pcSce[,which(pcSce$llpcSingler == "prePB" | 
                             (pcSce$Specific == "TRUE" & 
                                  pcSce$cellType == "B"))]

xbp1Dat <- logcounts(prePBSce)[which(rowData(prePBSce)$hgnc_symbol == "XBP1"),]
cd20Dat <- normcounts(altExp(prePBSce, "flowData"))["CD20",]
xpb1Df <- data.frame(prePBSce$cellTypes, prePBSce$llpcSingler, 
                     prePBSce$Specific, xbp1Dat, cd20Dat)
xpb1Df <- xpb1Df[order(xpb1Df$xbp1Dat),]



