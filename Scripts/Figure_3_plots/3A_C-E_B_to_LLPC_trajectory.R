library(SingleCellExperiment)
library(SingleR)
library(scran)
library(scater)
library(ggplot2)
library(pheatmap)
library(ggforce)
library(epitools)
#Here, we are going to use LLPC differentiation data from Kassambara et al: 
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8102200/
#THis data has been downloaded in an xlsx format and then converted to .csv locally.

llpcDat <- read.csv("Data/Separation_of_PB_and_LLPC/41375_2021_1234_MOESM3_ESM.csv",
                    row.names = 1)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#First, we remove all genes that do not overlap between the two. 
commonEnsembl <- row.names(csfSce)[which(row.names(csfSce) %in% row.names(llpcDat))]
length(commonEnsembl)
#19234

#So now, we are going to create a singler model based on these four states, and
#then plot how our cells fare in that model, to see if there is a difference 
#between the specific and non-specific in this focus.
csfSceRed <- csfSce[which(row.names(csfSce) %in% commonEnsembl),]
llpcRed <- llpcDat[which(row.names(llpcDat) %in% commonEnsembl),]

llpcOrd <- llpcRed[match(row.names(csfSceRed), row.names(llpcRed)),]

identical(row.names(csfSceRed), row.names(llpcOrd)) #TRUE

#So now, the data is transformed
llpcSce <- SingleCellExperiment(assays = list(counts = llpcOrd[,9:20]),
                                      colData = data.frame("cellType" = gsub("|[1,2,3]", "", colnames(llpcOrd)[9:20])),
                                rowData = llpcOrd[,1:8])

#As this is not really a single cell dataset, we should probably use bulk normalisation
#methods. However, we are going to compare to single-cell data, and therefore we use the
#most similar methods. 

llpcSce <- computeSumFactors(llpcSce)
llpcSce <- logNormCounts(llpcSce)

assays(llpcSce) #counts logcounts

#Then, we cluster all the data using the singler method
singlerRes <- SingleR(logcounts(csfSce), logcounts(llpcSce), labels = llpcSce$cellType)

table(csfSce$Specific, singlerRes$labels, csfSce$cellType)
#ASC
#             MBC  PB  PC prePB
#  CASPR2       0  15  15     1
#  FALSE        0  16   1     0
#  LGI1         0  58  21     4
#  Not_tested   0 101  34     5

#B
#             MBC  PB  PC prePB
#  CASPR2       0   0   0     2
#  FALSE        9   0   0     1
#  LGI1         1   0   0     1
#  Not_tested 171   0   0    11

#So now, the data is transformed to fit into a heatmap
csfSce$llpcSingler <- singlerRes$labels
specSce <- csfSce[,-which(csfSce$Specific == "Not_tested")]
specSce$cellSpec <- paste0(specSce$cellType, "_", specSce$Specific)

specSce$llpcSingler <- factor(specSce$llpcSingler, levels = c("MBC", "prePB", "PB", "PC"))

dir.create("Results/Figure_3_plots/LLPC_analyses", recursive = TRUE)
#We also want to make a version of this that separates the donors. 
specSce$cellSpecDon <- factor(paste0(specSce$cellSpec, "_", specSce$donor),
                              levels = c("B_FALSE_1166", 
                                         "B_FALSE_1227", 
                                         "B_FALSE_1284",
                                         "B_LGI1_1166", 
                                         "B_LGI1_1227",
                                         "B_CASPR2_1284",
                                         "ASC_FALSE_1166", 
                                         "ASC_FALSE_1227", 
                                         "ASC_FALSE_1284",
                                         "ASC_LGI1_1166", 
                                         "ASC_LGI1_1227",
                                         "ASC_CASPR2_1284"))
specSinglerDonTab <- table(specSce$cellSpecDon, specSce$llpcSingler)

#And this is changed for a row-wise interpretation
specSinglerDonTabNorm <- apply(as.matrix(specSinglerDonTab), 1, function(x) x/max(x))
specSinglerDonTabLog <- t(log10(specSinglerDonTab+1))

pdf("Results/Figure_3_plots/LLPC_analyses/Heatmap_of_LLPC_dev_singler_per_donor.pdf", height = 2.5, width = 4)
pheatmap(specSinglerDonTabLog, cluster_rows = FALSE,
         cluster_cols = FALSE, display_numbers = t(as.matrix(specSinglerDonTab)),
         color = gray.colors(50, start =1, end = 0))
dev.off()

#Now, we do statistics. We will divide into LGI1 and CASPR2. 
bLGITab <- specSinglerDonTab[which(grepl("B", row.names(specSinglerDonTab)) &!
                                             grepl("1284", row.names(specSinglerDonTab))),1:2]
fisherDf <- data.frame("Negative" = c(bLGITab[1,1]+bLGITab[2,1], bLGITab[1,2]+bLGITab[2,2]),
                       "LGI1" = c(bLGITab[3,1]+bLGITab[4,1], bLGITab[3,2]+bLGITab[4,2]))
row.names(fisherDf) <- c("MBC", "prePB")
fisherDf
#      Negative LGI1
#MBC          5    1
#prePB        0    1

fisher.test(fisherDf, alternative = "greater") #p-value = 0.28

#And the CASPR2.
bCASPR2Tab <- specSinglerDonTab[which(grepl("B", row.names(specSinglerDonTab)) &
                                       grepl("1284", row.names(specSinglerDonTab))),1:2]

fisher.test(t(bCASPR2Tab), alternative = "greater") #p-value = 0.15

#And the combination: 
bTab <- specSinglerDonTab[grep("B", row.names(specSinglerDonTab)),1:2]

fisherDf <- data.frame("Negative" = c(sum(bTab[1:3,1]), sum(bTab[1:3,2])),
                       "Positive" = c(sum(bTab[4:6,1]), sum(bTab[4:6,2])))
row.names(fisherDf) <- c("MBC", "prePB")
fisherDf
#      Negative Positive
#MBC          9        1
#prePB        1        3
fisher.test(fisherDf, alternative = "greater") #p-value = 0.04

#####
#ASC
ascLGITab <- specSinglerDonTab[which(grepl("ASC", row.names(specSinglerDonTab)) &!
                                         grepl("1284", row.names(specSinglerDonTab))),2:4]
fisherDf <- data.frame("Negative" = c(sum(ascLGITab[1:2,1:2]), sum(ascLGITab[1:2,3])),
                       "LGI1" = c(sum(ascLGITab[3:4,1:2]), sum(ascLGITab[3:4,3])))
row.names(fisherDf) <- c("nonPC", "PC")
fisherDf
#      Negative LGI1
#nonPC       13   62
#PC           1   21

fisher.test(fisherDf, alternative = "greater") #p-value = 0.1196

#And the CASPR2.
ascCASPR2Tab <- specSinglerDonTab[which(grepl("ASC", row.names(specSinglerDonTab)) &
                                            grepl("1284", row.names(specSinglerDonTab))),2:4]

fishDf <- data.frame("Negative" = c(sum(ascCASPR2Tab[1,1:2]), sum(ascCASPR2Tab[1,3])),
                     "CASPR2" = c(sum(ascCASPR2Tab[2,1:2]), sum(ascCASPR2Tab[2,3])),
                     row.names = c("nonPC", "PC"))

#      Negative CASPR2
#nonPC        3     16
#PC           0     15

fisher.test(fishDf, alternative = "greater") #p-value = 0.1619

#And the combination: 
ascTab <- specSinglerDonTab[grep("ASC", row.names(specSinglerDonTab)),2:4]

fisherDf <- data.frame("Negative" = c(sum(ascTab[1:3,1:2]), sum(ascTab[1:3,3])),
                       "Positive" = c(sum(ascTab[4:6,1:2]), sum(ascTab[4:6,3])))
row.names(fisherDf) <- c("nonPC", "PC")
fisherDf
#      Negative Positive
#nonPC       16       78
#PC           1       36

fisher.test(fisherDf, alternative = "greater") #p-value = 0.02038

#Now, to correlate this to some internal characteristic of the data, we will use
#this data as input for the plotting of previously generated RNA velocity data
veloCytoRaw <- read.csv("Data/Velocity/velocyto_out/Raw_cells/obs.csv")
veloCytoRaw$Cell <- gsub("onefilepercell_JR...._._._..._and_others_.....:JR|.out.bam",
                         "", veloCytoRaw$CellID)
##Now, for plotting purposes, we will include the cell types here. 
veloCytoRaw$singler <- unlist(sapply(veloCytoRaw$Cell, function(x){
    if(x %in% colnames(csfSce)){
        as.character(csfSce$llpcSingler[which(colnames(csfSce) == x)])
    } else {
        "Undefined"
    }
}))

veloCytoRaw$excluded <- FALSE
veloCytoRaw$excluded[which(veloCytoRaw$singler == "Undefined")] <- TRUE
write.csv(veloCytoRaw, "Data/Velocity/velocyto_out/Raw_cells/obs_plus_cell_type.csv")

#Now, after considerable plotting in python, we import the cell cycle results, and
#overlay them here. 
pcPbDat <- specSce[,which(specSce$llpcSingler %in% c("PB", "PC"))]
veloCytoRes <- read.csv("Data/Velocity/velocyto_out/All_cells/obs.csv")

pcPbDat$veloCytoCycle <- sapply(colnames(pcPbDat), function(x){
    if(x %in% veloCytoRes$Cell){
        veloCytoRes$phase[which(veloCytoRes$Cell == x)]
    } else {
        NA
    }
})
table(pcPbDat$veloCytoCycle, pcPbDat$llpcSingler, useNA = "ifany")

#     MBC prePB PB PC
#G1     0     0 39 34
#G2M    0     0 23  0
#S      0     0 26  2
#<NA>   0     0  1  1
#So two cells are lost, as we lack information on their cell cycle stage. 

plotDatCycle <- as.data.frame(colData(pcPbDat))
plotDatCycle$llpcSingler <- factor(plotDatCycle$llpcSingler, levels = c("PB", "PC"))
plotDatCycle$MKI67 <- logcounts(pcPbDat)[which(rowData(pcPbDat)$hgnc_symbol == "MKI67"),]

#We sadly have to exclude two cells here, as they did not get a velocyto cycle stage 
#assigned to them. 
plotDatCycle <- plotDatCycle[-which(is.na(plotDatCycle$veloCytoCycle)),]

plotDatCycle$veloCytoCycle <- factor(plotDatCycle$veloCytoCycle, levels = c("G1", "S", "G2M"))
plotDatCycle$cycleStage <- "G1"
plotDatCycle$cycleStage[-which(plotDatCycle$veloCytoCycle == "G1")] <- "Other"
plotDatCycle$PB <- "nonPC"
plotDatCycle$PB[which(plotDatCycle$llpcSingler == "PC")] <- "PC"

fisherTab <- table(plotDatCycle$cycleStage, plotDatCycle$PB)
fisherTab
#        nonPC PC
#G1       39 34
#Other    49  2
fisher.test(fisherTab) #p-value = 5.793e-08

oddsratio(table(plotDatCycle$PB, plotDatCycle$cycleStage))$measure
#       odds ratio with 95% C.I.
#.       estimate       lower    upper
#nonPC 1.00000000          NA       NA
#PC    0.06317908 0.007297024 0.1831058

plotDatCycle$Ki67Pos <- "Neg"
plotDatCycle$Ki67Pos[which(plotDatCycle$MKI67> 0)] <- "Pos"

#Now, what about Ki67 in specifics and non?
pbDat <- plotDatCycle[which(plotDatCycle$llpcSingler == "PB"),]
table(pbDat$Specific, pbDat$Ki67Pos)
#        Neg Pos
#FALSE   7  9
#TRUE   31  41

fisher.test(table(pbDat$Specific, pbDat$Ki67Pos))
#p-value 1. 

#Now some plotting. 
locPlotDat <- as.data.frame(table(plotDatCycle$Ki67Pos, 
                                  plotDatCycle$llpcSingler, 
                                  plotDatCycle$Specific))
colnames(locPlotDat) <- c("Ki67Pos", "Cell_type", "Specific", "Freq")
locPlotDat$Specific <- factor(locPlotDat$Specific, levels = c("FALSE", "LGI1", "CASPR2"))
p <- ggplot(locPlotDat,                         
            aes(x = Cell_type,
                y = Freq,
                fill = Ki67Pos)) + 
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    facet_grid(~ Specific) +scale_fill_manual(values = c("grey", "#696969")) +
    scale_y_continuous(expand = c(0, 0, 0.05,0))
p
ggsave("Results/Figure_3_plots/LLPC_analyses/Ki67Pos_fractions.pdf", plot = p)

p + theme_void() + theme(legend.position="none",
                         panel.grid.major = element_line())
ggsave("Results/Figure_3_plots/LLPC_analyses/Ki67Pos_fractions_no_legend.pdf",
       width = 5, height = 5.5)

table(plotDatCycle$Ki67Pos, plotDatCycle$llpcSingler)
#      PB PC
#FALSE 38 33
#TRUE  50  3
fisher.test(table(plotDatCycle$Ki67Pos, plotDatCycle$llpcSingler))
#p-value 2.666e-07

#Strangely, the data should be transposed to work in the oddsratio function
oddsratio(table(plotDatCycle$llpcSingler, plotDatCycle$Ki67Pos))
#odds ratio with 95% C.I.
#.    estimate       lower     upper
#PB 1.00000000          NA        NA
#PC 0.0734367 0.01606675 0.2261905

#As these plots were very successful visually, we will make the same kind also for the cell cycle
locPlotDat <- as.data.frame(table(plotDatCycle$veloCytoCycle, as.character(plotDatCycle$llpcSingler), plotDatCycle$Specific))
colnames(locPlotDat) <- c("Cell_cycle", "Cell_type", "Specific", "Freq")
locPlotDat$Specific <- factor(locPlotDat$Specific, levels = c("FALSE", "LGI1", "CASPR2"))

p <- ggplot(locPlotDat,                         
            aes(x = Cell_type,
                y = Freq,
                fill = Cell_cycle)) + 
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    facet_grid(~ Specific) +scale_fill_manual(values = c("white", "#999999", "#303030")) +
    scale_y_continuous(expand = c(0, 0, 0.05,0))
p
ggsave(paste0("Results/Figure_3_plots/LLPC_analyses/Cell_cycle_distribution.pdf"), plot = p)

p + theme_void() + theme(legend.position="none",
                         panel.grid.major = element_line())
ggsave("Results/Figure_3_plots/LLPC_analyses/Cell_cycle_distribution_no_legend.pdf",
       width = 5, height = 5.5)

#Now, we will produce plots with all the donors separated. 
locPlotDat <- as.data.frame(table(plotDatCycle$Ki67Pos, 
                                  plotDatCycle$llpcSingler, 
                                  plotDatCycle$Specific, 
                                  plotDatCycle$donor))
colnames(locPlotDat) <- c("Ki67Pos", "Cell_type", "Specific", "Donor", "Freq")
locPlotDat$Specific <- factor(locPlotDat$Specific, levels = c("FALSE", "LGI1", "CASPR2"))
p <- ggplot(locPlotDat,                         
            aes(x = Cell_type,
                y = Freq,
                fill = Ki67Pos)) + 
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    facet_grid(~ Donor+Specific) +scale_fill_manual(values = c("grey", "#696969")) +
    scale_y_continuous(expand = c(0, 0, 0.05,0))
p
ggsave("Results/Figure_3_plots/LLPC_analyses/Ki67Pos_fractions_per_donor.pdf", plot = p)

p + theme_void() + theme(legend.position="none",
                         panel.grid.major = element_line())
ggsave("Results/Figure_3_plots/LLPC_analyses/Ki67Pos_fractions_per_donor_no_legend.pdf",
       width = 5, height = 5.5)

#Now per-donor plots for the cell cycle. 
locPlotDat <- as.data.frame(table(plotDatCycle$veloCytoCycle, 
                                  as.character(plotDatCycle$llpcSingler), 
                                  plotDatCycle$Specific, 
                                  plotDatCycle$donor))
colnames(locPlotDat) <- c("Cell_cycle", "Cell_type", "Specific", "Donor", "Freq")
locPlotDat$Specific <- factor(locPlotDat$Specific, levels = c("FALSE", "LGI1", "CASPR2"))

p <- ggplot(locPlotDat,                         
            aes(x = Cell_type,
                y = Freq,
                fill = Cell_cycle)) + 
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    facet_grid(~ Donor+Specific) +scale_fill_manual(values = c("white", "#999999", "#303030")) +
    scale_y_continuous(expand = c(0, 0, 0.05,0))
p
ggsave(paste0("Results/Figure_3_plots/LLPC_analyses/Cell_cycle_distribution_per_donor.pdf"), plot = p)

p + theme_void() + theme(legend.position="none",
                         panel.grid.major = element_line())
ggsave("Results/Figure_3_plots/LLPC_analyses/Cell_cycle_distribution_per_donor_no_legend.pdf",
       width = 5, height = 5.5)

#And we also save the specific SCE with the LLPC data included
saveRDS(specSce, "Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")
