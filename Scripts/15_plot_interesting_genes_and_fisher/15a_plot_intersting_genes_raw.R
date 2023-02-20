library(SingleCellExperiment)
library(scater)
#Here, we are going to produce individual gene plots for a few marker genes as 
#well as the three top genes in the differential expression analysis
#After much going back and forth, we will here show DPEP1 and TCN2
focusGenes <- c("CD19", "MS4A1", "XBP1", "SDC1", "IRF4", "PRDM1", "CD79A", "DPEP1", "TCN2",
                "HDAC6", "PLEKHO2", "PLCD1",  "H6PD",  "KLC4", "CD27")
# PRDM1 = BLIMP1
# MS4A1 = CD20
sceList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

#We need to import the data for the AE, as the data in the other dataset is incomplete
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
aeSceRed <- aeSce[which(rowData(aeSce)$hgnc_symbol %in% row.names(sceList$AE_CSF)),]
row.names(aeSceRed) <- rowData(aeSceRed)$hgnc_symbol
aeSceOrd <- aeSceRed[match(row.names(sceList$AE_CSF), row.names(aeSceRed)),]
sceList$AE_pos <- aeSceOrd[,which(aeSceOrd$Specific == "TRUE")]
sceList$AE_neg <- aeSceOrd[,which(aeSceOrd$Specific == "FALSE")]
sceList$AE_CSF <- NULL

sceListRed <- lapply(sceList, function(x) x[which(row.names(x) %in% focusGenes),])

plotDatList <- lapply(names(sceListRed), function(x){
    locDat <- sceListRed[[which(names(sceListRed) == x)]]
    list(logcounts(locDat), "colDat" = DataFrame("dataset" = rep(x, ncol(locDat)),
                                                 "cellType" = locDat$cellType))
})

plotSce <- SingleCellExperiment(assays = list("logcounts" = do.call("cbind", lapply(plotDatList, "[[", 1))),
                                              colData = do.call("rbind", lapply(plotDatList, "[[", 2)))

plotSce$dataCellType <- paste0(plotSce$dataset, ".", plotSce$cellType)

plotSce$dataCellType <- factor(plotSce$dataCellType, levels = c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                                                "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                                                "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                                                "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B",
                                                                "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                                                "Covid_PBMC.ASC", "Covid_PBMC.B",
                                                                "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                                                "Influensa_PBMC.ASC", "Influensa_PBMC.B"))

dir.create("Results/Gene_expression_plots/Raw")

for(i in focusGenes[-which(focusGenes %in% c("DPEP1", "TCN2", "HDAC6", "PLEKHO2", "PLCD1",  "H6PD",  "KLC4"))]){
    p <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + ylim(0, 15) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) 
    ggsave(paste0("Results/Gene_expression_plots/Raw/", i, ".pdf"), plot = p, height = 7, width = 8)
    
    q <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + 
        theme_void() + ylim(0, 17) + 
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              strip.text.x = element_blank(),
              legend.position="none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
    ggsave(paste0("Results/Gene_expression_plots/Raw/", i, ".png"), plot = q, height = 7, width = 8)
}

#We create a second one for the hit genes (deprecated as we do not see these any longer)
for(i in c("DPEP1", "TCN2",  "HDAC6", "PLEKHO2", "PLCD1",  "H6PD",  "KLC4")){
    p <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + ylim(0, 12) + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
    ggsave(paste0("Results/Gene_expression_plots/Raw/", i, ".pdf"), plot = p, height = 7, width = 4)
    
    q <- plotExpression(plotSce, features = i, colour_by = "cellType", x = "dataCellType") + 
        theme_void() + 
        theme(axis.text.x = element_blank(), 
              axis.ticks.x = element_blank(), 
              axis.text.y = element_blank(),
              axis.ticks.y = element_blank(), 
              strip.text.x = element_blank(),
              legend.position="none",
              axis.title.x = element_blank(),
              axis.title.y = element_blank())
    ggsave(paste0("Results/Gene_expression_plots/Raw/", i, ".png"), plot = q, height = 1, width = 7)
}


