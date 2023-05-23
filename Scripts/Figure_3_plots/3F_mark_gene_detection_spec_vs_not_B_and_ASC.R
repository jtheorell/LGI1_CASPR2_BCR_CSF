library(scran)
library(scater)
library(SingleCellExperiment)
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


#CHANGE THIS TO PURE GGPLOT2
#Ok, so now the question is: can we find a direction through this data that is 
#by the specific B cells and the specific ASC? To do this, we need to tweak the
#algorithm slightly, as we are only interested in paired analyses. 
sceList <- list("bSce" = aeSceHgv[,which(aeSceHgv$cellType == "B")],
     "ascSce" = aeSceHgv[,which(aeSceHgv$cellType == "ASC")])

hitMarkers <- lapply(sceList, function(x){
    marker.info <- scoreMarkers(x, x$Specific)
    chosen <- marker.info$"TRUE"
    list("Hit" = list("Up" = row.names(chosen)[which(chosen$median.logFC.detected > 1)],
                      "Down" = row.names(chosen)[which(chosen$median.logFC.detected < -1)]),
    "All" = data.frame("Gene_name" = row.names(chosen), 
                       "Data" = chosen$median.logFC.detected))
})

upHit <- hitMarkers[[2]]$Hit$Up[which(hitMarkers[[2]]$Hit$Up %in% hitMarkers[[1]]$Hit$Up)]
#MFN2

downHit <- hitMarkers[[2]]$Hit$Down[which(hitMarkers[[2]]$Hit$Down %in% hitMarkers[[1]]$Hit$Down)]
#None

genePlotDf <- data.frame("Gene_names" = hitMarkers[[1]]$All$Gene_name, 
                         "B" = hitMarkers$bSce$All$Data, 
                         "ASC" = hitMarkers$ascSce$All$Data)
genePlotDf$Hit <- FALSE
genePlotDf$Hit[which(genePlotDf$Gene_names == upHit)] <- TRUE
genePlotDfOrdered <- genePlotDf[order(genePlotDf$Hit),]

ggplot(genePlotDfOrdered, aes(x = B, y = ASC, color = Hit)) + 
    geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("grey", "red"))
ggsave("Results/Figure_3_plots/Hit_in_B_vs_ASC.pdf", width = 5, height = 6)

#As the last hit is gone:
##And we calculate the fisher p-value for this individual one. 
#lapply(sceList, function(x){
#    nonSpecNonExpr <- length(which(counts(x)[which(rownames(x) == upHit),which(x$Specific == "FALSE")] == 0))
#    nonSpecExpr <- length(which(counts(x)[which(rownames(x) == upHit),which(x$Specific == "FALSE")] > 0))
#    specNonExpr <- length(which(counts(x)[which(rownames(x) == upHit),which(x$Specific == "TRUE")] == 0))
#    specExpr <- length(which(counts(x)[which(rownames(x) == upHit),which(x$Specific == "TRUE")] > 0))
#    fisherDf <-  data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
#                          "Spec_ASC" = c(specNonExpr, specExpr))
#    row.names(fisherDf) <- c("Not_expressed", "Expressed")
#    fisher.test(fisherDf)
#})
#
##bSce
##0.1912
##
##ascSce
##0.1149
#
##So, not significant, but apporaching. Now, we will plot it, to show this point loud and clear. 
##And now, the one gene that is identified as "specific"
#aeSceHgv$specificCellType <- paste0(aeSceHgv$cellType, "_", aeSceHgv$Specific)
#
#p <- plotExpression(aeSceHgv, features = "MFN2", colour_by = "Specific", x = "specificCellType",
#                    size_by = "lowQual") + 
#    scale_color_manual(values = c("black", "orange")) + scale_size_manual(values = 3) +
#    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
#ggsave(paste0("Results/MFN2_expression.pdf"), plot = p)
#
#p <- p + 
#    theme_void() + 
#    theme(axis.text.x = element_blank(), 
#          axis.ticks.x = element_blank(), 
#          axis.text.y = element_blank(),
#          axis.ticks.y = element_blank(), 
#          strip.text.x = element_blank(),
#          legend.position="none",
#          axis.title.x = element_blank(),
#          axis.title.y = element_blank())
#ggsave(paste0("Results/MFN2_expression_stripped.png"), plot = p, height = 5, width = 5)
#
#