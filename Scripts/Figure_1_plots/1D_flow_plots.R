library(ggplot2)
library(scater)
library(scran)

flowDat <- read.csv("Data/Cytometry/flowDataPlusIndex.csv")
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
TCRdb <- read.csv("Data/TCR_data/cell_data_all_together.csv")
TCRdb$Sample <- substr(TCRdb$cell_name, 3,8)
TCRdb$Cell <- substr(TCRdb$cell_name, 3,16)

#Here, we do not exclude possibly incomplete cells, as we are interested in getting an
#as rich count as possible. 
flowDataT_or_B <- flowDat[which(flowDat$Cell %in% TCRdb$Cell | 
                                 flowDat$Cell %in% colnames(aeSce)),]

dir.create("Results/Figure_1_plots/Flow_specific", recursive =TRUE)
#Now, we are going to follow up what we did in the beginning for the whole dataset
cellType <- rep("other", nrow(flowDataT_or_B))
cellType[which(flowDataT_or_B$CD3 > 1.25 & flowDataT_or_B$CD138 < 2.5)] <- "T.cell"
cellType[which(cellType == "T.cell" & flowDataT_or_B$CD19 > 2)] <- "other"

cellType[which(flowDataT_or_B$Cell %in% colnames(aeSce)[which(aeSce$cellType == "B")])] <- "B"
cellType[which(flowDataT_or_B$Cell %in% colnames(aeSce)[which(aeSce$cellType == "ASC")])] <- "ASC"

cellType[which(cellType == "T.cell")] <- "DNT"

cellType[which(cellType == "DNT" & flowDataT_or_B$CD4 > 2.2 & flowDataT_or_B$CD8 > 2.5)] <- "DPT"
cellType[which(cellType == "DNT" & flowDataT_or_B$CD4 > 2.2 & flowDataT_or_B$CD8 <= 2.5)] <- "CD4T"
cellType[which(cellType == "DNT" & flowDataT_or_B$CD4 <= 2.2 & flowDataT_or_B$CD8 > 2.5)] <- "CD8T"


flowDataT_or_B$cellType <- cellType

#For these plots, we will exclude the "others". 
flowDatRed <- flowDataT_or_B[-which(cellType %in% c("other", "CD8T", "DNT", "DPT")),]

#Now, there is one cell on the CD38 that is very, very much below zero, and it is
#therefore truncated
flowDatRed$CD38[which(flowDatRed$CD38 < -3)] <- -2.5

#We will also have to add CD27 RNA, which is not entirely simple, as the T-cell
#transcriptomes have not been analysed in any detail. Here, the same exclusions
#as used for the B cell data in the beginning is used. 
csfSce <- readRDS("../External/Data/csfSce_1b_preQC_1166_1227_1284_exp_1.rds")
#Identification of mitochondrial genes
mito <- which(rowData(csfSce)$chromosome_name=="MT")
#Now the preCellQCMetrcis are calculated. 
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))
colData(csfSce) <- cbind(colData(csfSce), stats)
lowCount <- fewFeatures <- highMito <- highERCC <- rep(FALSE, ncol(csfSce))
lowCount[which(csfSce$sum < 100000)] <- TRUE
fewFeatures[which(csfSce$detected < 1000)] <- TRUE
highMito[which(csfSce$subsets_Mt_percent > 10)] <- TRUE
#For the ERCC, there is a considerably larger population of cells with high percentages 
#in the JR1166 donor. 
highERCC[which(csfSce$altexps_ERCC_percent > 30)] <- TRUE
csfSce$lowQual <- FALSE
csfSce$lowQual[which(lowCount | fewFeatures | highMito | highERCC)] <- TRUE
csfSce <- csfSce[,!csfSce$lowQual]
csfSce <- computeSumFactors(csfSce) 
csfSce <- logNormCounts(csfSce)

CD27Dat <- logcounts(csfSce)[which(rowData(csfSce)$hgnc_symbol == "CD27"),]

flowDatRedRed <- flowDatRed[which(flowDatRed$Cell %in% colnames(csfSce)),]
#We also need to remove all A1. 
flowDatRedRed <- flowDatRedRed[-grep("A01", flowDatRedRed$Cell),]

names(CD27Dat) <- colnames(csfSce)
CD27DatRed <- CD27Dat[which(names(CD27Dat)  %in% flowDatRedRed$Cell)]
CD27DatRedOrd <- CD27DatRed[match(flowDatRedRed$Cell, names(CD27DatRed))]
identical(names(CD27DatRedOrd), flowDatRedRed$Cell) #TRUE
flowDatRedRed$CD27RNA <- CD27DatRedOrd

plotList <- list(c("CD3", "CD19"),
                 c("CD4", "CD8"),
                 c("CD20", "CD38"),
                 c("CD138", "CD27RNA"),
                 c("CD20", "CD138"),
                 c("FSC.A", "SSC.A"))

for(i in plotList){
    locDat <- flowDatRedRed[,which(colnames(flowDatRedRed) %in% c(i, "cellType"))]
    p <- ggplot(locDat, aes(x=locDat[,which(colnames(locDat) == i[1])], 
                            y=locDat[,which(colnames(locDat) == i[2])], 
                            shape = cellType, fill = cellType)) + 
        geom_point(size = 4) + theme_bw() + scale_fill_manual(values = c("black", "white", "grey")) +
        scale_shape_manual(values = c(24,22,21))
    p
    ggsave(paste0("Results/Figure_1_plots/Flow_specific/", i[1], "_vs_", i[2], ".pdf"))
    p + theme(legend.position = "none", 
              panel.grid.major = element_blank(), 
              panel.grid.minor = element_blank(),
              axis.text.x = element_blank(),
              axis.text.y = element_blank(),
              strip.text.x = element_blank(),
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
    ggsave(paste0("Results/Figure_1_plots/Flow_specific/", i[1], "_vs_", i[2], "_stripped.pdf"),
           width = 5, height = 5)
}

#For CD4/CD8, we also need to include a picture with the correct axes
i <- plotList[[2]]
locDat <- flowDatRedRed[,which(colnames(flowDatRedRed) %in% c(i, "cellType"))]
p <- ggplot(locDat, aes(x=locDat[,which(colnames(locDat) == i[1])], 
                        y=locDat[,which(colnames(locDat) == i[2])], 
                        shape = cellType, fill = cellType)) + 
    geom_point(size = 4) + theme_bw() + scale_fill_manual(values = c("black", "white", "grey")) +
    scale_shape_manual(values = c(24,22,21)) + scale_y_continuous(limits = c(-3, 9), expand = c(0,0))
p
ggsave(paste0("Results/Figure_1_plots/Flow_specific/", i[1], "_vs_", i[2], ".pdf"))
p + theme(legend.position = "none", 
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          axis.text.x = element_blank(),
          axis.text.y = element_blank(),
          strip.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
ggsave(paste0("Results/Figure_1_plots/Flow_specific/", i[1], "_vs_", i[2], "_stripped.pdf"),
       width = 5, height = 5)
