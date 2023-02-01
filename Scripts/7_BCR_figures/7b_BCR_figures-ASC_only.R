library(ggnetwork)
library(viridis)
library(scatterpie)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

#Before doing anything else, we will exclude all cells with an unknown specificity, 
#as they are not useful in this analysis. 
BCR_spec <- BCR_all[-which(BCR_all$Specific == "Not_tested"),]

#We start by excluding all non-ASC, as we show in the previous subfigure that
#the specific cells are confined to this subset. 

csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

BCR_asc <- BCR_spec[which(BCR_spec$CELL %in% colnames(csfSce)),]


#For the visualisations, we do not need the light chains, as
#the heavy chains are already categorised according to light chain usage, and the
#light chain distances, that will be used for collapsing BCRs, are also present 
#for the heavy chains. 
#THis step is unnecessary here, as the light chains are not IgG. 
BCR_H <- BCR_asc[which(BCR_asc$LOCUS == "H"),]

#Now, we are going to start with anyway taking in information from the light chains: 
#we will add a complete mutational score for each cell
BCR_H$mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_asc$All_mutations[which(BCR_asc$CELL == x)])
})

#We also add a column with the CDR3 length
BCR_H$CDR3_length <- sapply(BCR_H$CDR3_IMGT, nchar)

#Here, we are going to add a dummy donor will all the info from all donors, 
#that will be displayed in the paper. 
BCR_H2 <- BCR_H
BCR_H2$Sample <- "All"
BCR_H_double <- rbind(BCR_H, BCR_H2)
#Now over to the plotting. 

#Now, we will create simple stacked bar graphs for the above parameters. 
barGraphDat <- BCR_H_double[,c("Sample","Clonal", "mutations_H_L", "All_mutations", "light_type", "ISOTYPE", "CDR3_length",
                        "Specific")]
barGraphDat$light_type <- factor(barGraphDat$light_type, levels = c("K", "L"))
barGraphDat$ISOTYPE <- factor(barGraphDat$ISOTYPE, levels = c("IGHG1", "IGHG2", "IGHG3", "IGHG4"))

barGraphDat$HLMutBinned <- "0"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 1:5)] <- "1-5"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 6:10)] <- "6-10"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 11:15)] <- "11-15"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 16:20)] <- "16-20"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 21:25)] <- "21-25"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 26:30)] <- "26-30"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L %in% 31:35)] <- "31-35"
barGraphDat$HLMutBinned[which(barGraphDat$mutations_H_L > 35)] <- ">35"
barGraphDat$HLMutBinned <- factor(barGraphDat$HLMutBinned, levels = c(">35", "31-35",
                                                                      "26-30", "21-25",
                                                                      "16-20", "11-15",
                                                                      "6-10", "1-5", "0"))

barGraphDat$MutBinned <- "0"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 1:5)] <- "1-5"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 6:10)] <- "6-10"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 11:15)] <- "11-15"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 16:20)] <- "16-20"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 21:25)] <- "21-25"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 26:30)] <- "26-30"
barGraphDat$MutBinned[which(barGraphDat$All_mutations %in% 31:35)] <- "31-35"
barGraphDat$MutBinned[which(barGraphDat$All_mutations > 35)] <- ">35"
barGraphDat$MutBinned <- factor(barGraphDat$MutBinned, levels = c(">35", "31-35",
                                                                  "26-30", "21-25",
                                                                  "16-20", "11-15",
                                                                  "6-10", "1-5", "0"))

barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length < 30)] <- "<30"
barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length %in% 30:39)] <- "30-39"
barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length %in% 40:49)] <- "40-49"
barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length %in% 50:59)] <- "50-59"
barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length %in% 60:69)] <- "60-69"
barGraphDat$CDR3Binned[which(barGraphDat$CDR3_length > 70)] <- ">70"
barGraphDat$CDR3Binned <- factor(barGraphDat$CDR3Binned, levels = c(">70", "60-69",
                                                                    "50-59",  "40-49",
                                                                    "30-39", "<30"))


colsAndNames <- list(list("Clonal", c("grey","black")),
                     list("HLMutBinned", rev(viridis(9))),
                     list("MutBinned", rev(viridis(9))),
                     list("light_type", c("#D40041", "#2E00B7")),
                     list("ISOTYPE", c("#B1CFFF", "#9E4532", "#541352", "#238A8D")),
                     list("CDR3Binned", rev(magma(6))))
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
    p + theme_void() + theme(legend.position="none") + scale_y_continuous(expand = c(0, 0))
    ggsave(paste0("Results/CloneBarGraphs/", i[[1]], "_no_legend.pdf"),
           width = 8, height = 6)
    
}


#FOR FIGURE: 
table(substr(BCR_H$CELL, 1,4))
#1166 1227 1284 
#89   13   35
