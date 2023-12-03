library(scran)
library(scater)
library(SingleCellExperiment)
library(ggplot2)
library(edgeR)
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
#We first need to go through some feature selection, to be able to run the below steps. 
sceVar <- modelGeneVar(aeSce, block = aeSce$cellType)
chosen <- getTopHVGs(sceVar, prop=0.2)
str(chosen) #[1:1450] ...

aeSceSpec <- aeSce[,-which(aeSce$Specific == "Not_tested")]

aeSceHgv <- aeSceSpec[chosen,]
row.names(aeSceHgv) <- rowData(aeSceHgv)$hgnc_symbol

aeSceHgv <- runUMAP(aeSceHgv )

plotReducedDim(aeSceHgv, dimred="UMAP", colour_by="Specific", 
               shape_by = "cellType") + scale_color_manual(values = c("black", "orange"))
library(ggplot2)
umapPlotDat <- data.frame(reducedDim(aeSceHgv, "UMAP"), "Specific" = aeSceHgv$Specific,
                          "Cell_type" = aeSceHgv$cellType)

ggplot(umapPlotDat, aes(x = X1, y = X2, 
                           color = Specific, 
                           shape = Cell_type)) + 
    geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("black", "orange")) +
    scale_shape_manual(values = c(17,15))
#ggsave("Results/UMAP_with_Spec_and_cellType.pdf", width = 5, height = 6)
#This umap is quite useless, so we will not show it. 

#Ok, so now the question is: can we find a direction through this data that is 
#by the specific B cells and the specific ASC? We will both look at the individual
#and the composite ASC vs B here. We will start with the sub-question - how
#strong is the correlation between the identified expression differences and 
#the mem-B vs pre-PB as well as PB vs PC expression differences, known from before?
llpcDat <- read.csv("Data/Separation_of_PB_and_LLPC/41375_2021_1234_MOESM3_ESM.csv",
                    row.names = 1)
llpcRed <- llpcDat[which(row.names(llpcDat) %in% rowData(aeSceHgv)$ensembl_gene_id),]
#So nine genes are lost here and we go from 1413 to 1404. This is accetable.

aeSceHgvRed <- aeSceHgv[which(rowData(aeSceHgv)$ensembl_gene_id %in% row.names(llpcDat)),]

llpcRedOrd <- llpcRed[match(rowData(aeSceHgvRed)$ensembl_gene_id, row.names(llpcRed)),]
identical(row.names(llpcRedOrd), rowData(aeSceHgvRed)$ensembl_gene_id) #TRUE
#Therefore, we can now change the row names.  
row.names(llpcRedOrd) <- row.names(aeSceHgvRed)

#And turn this into a singleCellExperiment
llpcRedOrdSce <- SingleCellExperiment(assays = list("counts" = llpcRedOrd[,grep("MBC|pre|PB|PC", colnames(llpcRedOrd))]), 
                                      rowData = llpcRedOrd[,-grep("MBC|pre|PB|PC", colnames(llpcRedOrd))],
                                      colData = data.frame("CellType" = c(rep("MBC", 3),
                                                                          rep("pre", 3),
                                                                          rep("PB", 3),
                                                                          rep("PC", 3))))

sceList <- list("bSce" = list("SingleCell" = aeSceHgvRed[,which(aeSceHgvRed$cellType == "B")],
                              "Bulk" = llpcRedOrdSce[,grep("MBC|pre", llpcRedOrdSce$CellType)],
                              "blocking" = NULL),
     "ascSce" = list("SingleCell" = aeSceHgvRed[,which(aeSceHgvRed$cellType == "ASC")],
                     "Bulk" = llpcRedOrdSce[,grep("PB|PC", llpcRedOrdSce$CellType)],
                     "blocking" = NULL),
     "allSce" = list("SingleCell" = aeSceHgvRed,
                     "Bulk" = llpcRedOrdSce[,grep("PB|PC", llpcRedOrdSce$CellType)],
                     "blocking" = aeSceHgvRed$cellType))

dir.create("Results/Figure_3_plots/3D-E_transcriptome")
hitMarkers <- lapply(names(sceList), function(i){
    x <- sceList[[which(names(sceList) == i)]]$SingleCell
    marker.info <- scoreMarkers(x, x$Specific, block = sceList[[which(names(sceList) == i)]]$blocking)
    chosen <- marker.info$"TRUE"
    #We also calculate Mann-Whitney p-values for all of the genes
    specDat <- logcounts(x[,which(x$Specific == "TRUE")])
    nonSpecDat <- logcounts(x[,which(x$Specific == "FALSE")])
    pVals <- sapply(seq_len(nrow(specDat)), function(y){
        wilcox.test(specDat[y,], nonSpecDat[y,], exact = FALSE, paired = FALSE)$p.value
    })
    
    #Now, we will also calculate the same values for the MBC vs prePB
    z <- sceList[[which(names(sceList) == i)]]$Bulk
    zDge <- DGEList(counts(z), samples=colData(z))

    zDge <- calcNormFactors(zDge)
    design <- model.matrix(~factor(CellType), zDge$samples)
    
    #design
    zDge <- estimateDisp(zDge, design)
    #summary(y$trended.dispersion)
    fit <- glmQLFit(zDge, design, robust=TRUE)
    #summary(fit$var.prior)
    res <- glmQLFTest(fit, coef=ncol(design))
    
    #And we retrieve the logFC values from this analysis for the plotting below
    
    #fdrCorrPVals <- p.adjust(pVals, method = "fdr")
    plotDat <- data.frame("LogFC" = chosen$median.logFC.detected, 
                          "negLog_non_adjusted_MW_p" = -log10(pVals),
                          "LogFC_in_bulk" = res$table$logFC)
    plotDat$Hit <- FALSE
    plotDat$Hit[which(abs(plotDat$LogFC) > 2)] <- TRUE
    p <- ggplot(plotDat, aes(x = LogFC, y = negLog_non_adjusted_MW_p, color = LogFC_in_bulk)) + 
        geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
        scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0,0)) +
        scale_y_continuous(limits = c(-0.3, 3.5), expand = c(0,0)) + 
        scale_color_viridis()
    ggsave(paste0("Results/Figure_3_plots/3D-E_transcriptome/", i, "_volcano.pdf"), plot = p,
           width = 5, height = 6)
    p <- ggplot(plotDat, aes(x = LogFC, y = negLog_non_adjusted_MW_p, color = Hit)) + 
        geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
        scale_color_manual(values = c("grey", "red")) + 
        scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0,0)) +
        scale_y_continuous(limits = c(-0.3, 3.5), expand = c(0,0))
    ggsave(paste0("Results/Figure_3_plots/3D-E_transcriptome/", i, "_volcano_hit_colors.pdf"), plot = p,
           width = 5, height = 6)
    p <- ggplot(plotDat, aes(x = LogFC, y = LogFC_in_bulk, color = LogFC_in_bulk)) + 
        geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
        scale_color_viridis()
    ggsave(paste0("Results/Figure_3_plots/3D-E_transcriptome/", i, "_FC_vs_FC.pdf"), plot = p,
           width = 5, height = 6)
    #Now, we are going to run a correlation analysis on the more significant genes. 
    quant95AbsFCBulk <- quantile(plotDat$LogFC_in_bulk, 0.95)
    quant95AbsFCSc <- quantile(plotDat$LogFC, 0.95)
    highDat <- plotDat[which(abs(plotDat$LogFC_in_bulk) > quant95AbsFCBulk &
                                 abs(plotDat$LogFC) > quant95AbsFCSc),]
    cor(highDat$LogFC, highDat$LogFC_in_bulk)
    
    list("Hit" = list("Up" = row.names(chosen)[which(chosen$median.logFC.detected > 2)],
                      "Down" = row.names(chosen)[which(chosen$median.logFC.detected < -2)],
                      "Up_FC1" = row.names(chosen)[which(chosen$median.logFC.detected > 1)],
                      "Down_FC1" = row.names(chosen)[which(chosen$median.logFC.detected < -1)]),
    "All" = data.frame("Gene_name" = row.names(chosen), 
                       "Sc_logFC" = chosen$median.logFC.detected,
                       "Bulk_logFC" = res$table$logFC),
    "Correlation_of_quantile_95_hits_to_bulk" = cor(highDat$LogFC, highDat$LogFC_in_bulk))
})

names(hitMarkers) <- names(sceList)

#Now, we first investigate if the identified hits at abs(logFC) > 2 are higher in the
#expected pattern in the bulk data.

hitMarkers$bSce$Hit$Up
#"FBH1"  "TPST2"
hitMarkers$bSce$Hit$Down
#"TSPAN3"
hitMarkers$ascSce$Hit$Up
#GRAP
hitMarkers$ascSce$Hit$Down
#None
hitMarkers$allSce$Hit$Up
#"GRAP"
hitMarkers$allSce$Hit$Down
#None

#So all in all, we are left with GRAP. 
#When looking at the bulk data: 
counts(llpcRedOrdSce)[which(row.names(llpcRedOrdSce) == "GRAP"),]
#     MBC1 MBC2 MBC3 prePB1 prePB2 prePB3 PB1 PB2 PB3 PC1 PC2 PC3
#GRAP  481  199  545    210    270    144 136 112  86  86  75  88
#It seems like there is a lowering of the protein along the differentiation
#axis. This expression therefore seems to go against the differentiation
#trend. 

#In the single cell data: 
aeSceLlpc <- readRDS("Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")
grapDat <- aeSceLlpc[which(rowData(aeSceLlpc)$hgnc_symbol == "GRAP"),]

singlerSpecGrap <- 
    data.frame("Cell_type" = c(rep("MBC", 2), rep("PB",2), rep("PC", 2), rep("prePB",2)),
               "Specific" = c(rep(c(FALSE, TRUE),4)),
        "value" = round(unlist(lapply(split(logcounts(grapDat), paste0(grapDat$llpcSingler, "_", grapDat$Specific)),mean)),2))

cellTypeSpecGrap <- 
    data.frame("Cell_type" = c(rep("B", 2), rep("ASC",2)),
               "Specific" = c(rep(c(FALSE, TRUE),2)),
               "value" = round(unlist(lapply(split(logcounts(grapDat), paste0(grapDat$cellType, "_", grapDat$Specific)),mean)),2))


#Ok, so now then: 
hitMarkers$bSce$All$Bulk_logFC[which(hitMarkers$bSce$All$Gene_name %in% hitMarkers$bSce$Hit$Up)]
#0.07511737 2.24619733. 
#When investigating FBH1 with the Monaco dataset, there is a 
#Shift from 28.7 to 162.8 in expression from switched memory B cells to plasmablasts,
#or a logFC of 0.76. 
#For TPST2, the same data shows 277.6 for the plasmablasts and 59.9 for switched mem, 
#in other words 0.66, so it is somewhat surprising that there is not a starker
#difference for FBH1. 
hitMarkers$bSce$All$Bulk_logFC[which(hitMarkers$bSce$All$Gene_name %in% hitMarkers$bSce$Hit$Down)]
#-2.576183, so that one is clearly down. The same trend is clear in the Monaco data, 
#where the logFC is  log10(153.2/504.3) or -0.52.

#Now, what about GRAP: 
hitMarkers$ascSce$All$Bulk_logFC[which(hitMarkers$ascSce$All$Gene_name %in% hitMarkers$ascSce$Hit$Up)]
#0.2281301. So this is not very convincing. Here there is no good alternative 
#golden standard. We can however of course look at the cells that have been classified
#as PB and PC to see if htere is an expression difference there. 

mean(logcounts(grapDat)[which(grapDat$llpcSingler == "PB")])
#0.9885177
mean(logcounts(grapDat)[which(grapDat$llpcSingler == "PC")])
#1.01413
#Median is null for both. So this gene does not seem to be well explained by 
#the PB-PC differentiation axis alone. 

#Now, what about the main comparison, i.e. B vs ASC? 
upHit <- hitMarkers[[2]]$Hit$Up_FC1[which(hitMarkers[[2]]$Hit$Up_FC1 %in% hitMarkers[[1]]$Hit$Up_FC1)]
#None

downHit <- hitMarkers[[2]]$Hit$Down_FC1[which(hitMarkers[[2]]$Hit$Down_FC1 %in% hitMarkers[[1]]$Hit$Down_FC1)]
#None

#And this is when the threshold is reduced to logFC 1. 

genePlotDf <- data.frame("Gene_names" = hitMarkers[[1]]$All$Gene_name, 
                         "B" = hitMarkers$bSce$All$Data, 
                         "ASC" = hitMarkers$ascSce$All$Data)
genePlotDf$Hit <- FALSE
genePlotDf$Hit[which(genePlotDf$Gene_names == upHit)] <- TRUE
genePlotDfOrdered <- genePlotDf[order(genePlotDf$Hit),]

ggplot(genePlotDfOrdered, aes(x = B, y = ASC, color = Hit)) + 
    geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("grey", "red")) +
    scale_x_continuous(limits = c(-2.5, 2.5), expand = c(0,0)) +
    scale_y_continuous(limits = c(-2.5, 2.5), expand = c(0,0))
ggsave("Results/Figure_3_plots/3D-E_transcriptome/Hit_in_B_vs_ASC.pdf", width = 5, height = 6)

