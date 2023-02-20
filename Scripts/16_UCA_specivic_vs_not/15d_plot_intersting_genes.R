#Here, we are going to produce individual gene plots for a few marker genes as 
#well as the three top genes in the differential expression analysis
library(SingleCellExperiment)
library(ggplot2)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

aeSceSpecKnown <- aeSce[,-which(aeSce$Specific == "Not_tested")]

#Now, we will remove all without a HGNC symbol and then change the names to them
dupHgnc <- rowData(aeSce)$hgnc_symbol[which(duplicated(rowData(aeSce)$hgnc_symbol))]

table(dupHgnc)
#      ABCF2    AHRR   ATXN7    GGT1   PINX1 POLR2J3    TBCE 
#606       1       1       1       1       1       1       1
aeSceHgnc <- aeSceSpecKnown[-which(rowData(aeSceSpecKnown)$hgnc_symbol %in% dupHgnc),]
row.names(aeSceHgnc) <- rowData(aeSceHgnc)$hgnc_symbol

topResults <- read.csv("Results/Specific_UCA_vs_not/Top_results.csv", row.names = 1)

plotExpression(aeSceHgnc, 
                    features = rownames(topResults), x = "Specific", ncol = 4,
                    colour_by = "Specific") + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggsave("Results/Specific_UCA_vs_not/Marker_expression_plots.pdf", width = 5, height = 10)


