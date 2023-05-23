library(SingleCellExperiment)
library(SingleR)
library(scran)
library(scater)
library(ggplot2)
library(pheatmap)
library(ggforce)
library(epitools)
#Here, we are going to use the same data as in 22, from 
# this study: https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8102200/
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

#So now, the data is transo
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
#FALSE        0  16   1     0
#Not_tested   0 101  34     5
#TRUE         0  73  36     5

#B
#             MBC  PB  PC prePB
#FALSE       9   0   0     1
#Not_tested 171   0   0    11
#TRUE         1   0   0     3

#So now, the data is transformed to fit into a heatmap
csfSce$llpcSingler <- singlerRes$labels
specSce <- csfSce[,-which(csfSce$Specific == "Not_tested")]
specSce$cellSpec <- paste0(specSce$cellType, "_", specSce$Specific)

specSce$llpcSingler <- factor(specSce$llpcSingler, levels = c("MBC", "prePB", "PB", "PC"))
specSce$cellSpec <- factor(specSce$cellSpec, levels = c("B_FALSE", "B_TRUE", "ASC_FALSE", "ASC_TRUE"))
specSinglerTab <- table(specSce$cellSpec, specSce$llpcSingler)

#And this is changed for a row-wise interpretation
specSinglerTabNorm <- apply(as.matrix(specSinglerTab), 1, function(x) x/max(x))
specSinglerTabLog <- t(log10(specSinglerTab+1))


dir.create("Results/Figure_3_plots/LLPC_analyses", recursive = TRUE)
pdf("Results/Figure_3_plots/LLPC_analyses/Heatmap_of_LLPC_dev_singler.pdf", height = 2, width = 2.5)
pheatmap(specSinglerTabLog, cluster_rows = FALSE,
         cluster_cols = FALSE, display_numbers = t(specSinglerTab),
         color = gray.colors(50, start =1, end = 0))
dev.off()

#Now, we do statistics. For the B first: 
BTab <- as.matrix(specSinglerTab[1:2,])
fisherDf <- data.frame("Non_spec" = c(BTab[1,1], sum(BTab[1,2:4])),
                       "Spec" = c(BTab[2,1], sum(BTab[2,2:4])))
row.names(fisherDf) <- c("MBC", "prePB")

fisherDf
#        Non_spec Spec
#MBC           9    1
#prePB        1    3
fisher.test(fisherDf, alternative = "greater") #p-value = 0.04096

#Now, we do statistics. For the B first: 
ASCTab <- as.matrix(specSinglerTab[3:4,])
fisherDf <- data.frame("Non_spec" = c(sum(ASCTab[1,1:3]), ASCTab[1,4]),
                       "Spec" = c(sum(ASCTab[2,1:3]), ASCTab[2,4]))
row.names(fisherDf) <- c("prePC", "PC")

fisherDf
#        Non_spec Spec
#prePC       16   78
#PC           1   36
fisher.test(fisherDf, alternative = "greater") #p-value = 0.02038


#We also want to make a version of this that separates the donors. 
specSinglerDonTab <- table(paste0(specSce$cellSpec, "_", specSce$donor), specSce$llpcSingler)

#And this is changed for a row-wise interpretation
specSinglerDonTabNorm <- apply(as.matrix(specSinglerDonTab), 1, function(x) x/max(x))
specSinglerDonTabLog <- t(log10(specSinglerDonTab+1))

pdf("Results/Figure_3_plots/LLPC_analyses/Heatmap_of_LLPC_dev_singler_per_donor.pdf", height = 2.5, width = 4)
pheatmap(specSinglerDonTabLog, cluster_rows = FALSE,
         cluster_cols = FALSE, display_numbers = t(as.matrix(specSinglerDonTab)),
         color = gray.colors(50, start =1, end = 0))
dev.off()

#Now, as this turned out to be interesting, we are going to plot and do statistics
#on the Ki67 expression and cell cycle position for the PB and PC populations. 
pcPbDat <- specSce[,which(specSce$llpcSingler %in% c("PB", "PC"))]

veloCytoRes <- read.csv("Data/Velocity/velocyto_out/All_cells/obs.csv")
veloCytoRes$Cell <- gsub("onefilepercell_JR...._._._..._and_others_.....:JR|.out.bam",
                         "", veloCytoRes$CellID)

pcPbDat$veloCytoCycle <- sapply(colnames(pcPbDat), function(x){
    if(x %in% veloCytoRes$Cell){
        veloCytoRes$phase[which(veloCytoRes$Cell == x)]
    } else {
        NA
    }
})


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
#G1       36 34
#Other    52  3
fisher.test(fisherTab) #p-value = 8.594e-08

oddsratio(table(plotDatCycle$PB, plotDatCycle$cycleStage))$measure
#       odds ratio with 95% C.I.
#.       estimate       lower    upper
#nonPC 1.00000000          NA       NA
#PC    0.06317908 0.01384882 0.1936517

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

locPlotDat <- as.data.frame(table(plotDatCycle$Ki67Pos, plotDatCycle$llpcSingler, plotDatCycle$Specific))
colnames(locPlotDat) <- c("Ki67Pos", "Cell_type", "Specific", "Freq")

p <- ggplot(locPlotDat,                         
            aes(x = Cell_type,
                y = Freq,
                fill = Ki67Pos)) + 
    geom_bar(stat = "identity", color = "black") + theme_bw() +
    facet_grid(~ Specific) +scale_fill_manual(values = c("grey", "#696969")) +
    scale_y_continuous(expand = c(0, 0, 0.05,0))
p
ggsave(paste0("Results/Figure_3_plots/LLPC_analyses/Ki67Pos_fractions.pdf"), plot = p)

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

#And we also save the specific SCE with the LLPC data included
saveRDS(specSce, "Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")
