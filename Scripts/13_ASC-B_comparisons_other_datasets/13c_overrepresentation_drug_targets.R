#See all explanations under 12. These are made for comparative reasons. 
#We will here focus on seeing if the five genes that were identified in the top,
#i.e. CD27 CD38 CSF2RB  IL6R SDC1, are enriched in any systematic way. 

fileList <- list.files("Data/ASC_vs_B", pattern = "EdgeR_results.rds", 
                       recursive = TRUE, full.names = TRUE)

resList <- lapply(fileList, readRDS)

names(resList) <- gsub("Data/ASC_vs_B/|/EdgeR_results.rds", "", fileList)

drugTargetFileNameList <- list.files("Data/STRING_data", full.names = TRUE, pattern = ".tsv")
directTargs <- c("CD27", "CD38", "CSF2RB", "IL6R", "SDC1")
#Now, we loop this
topAmongDirectTargets <- lapply(resList, function(x){
    resTab <- x$table
    topGenes <- row.names(resTab)[which(resTab$logFC > 2 & resTab$PValue < 0.05)]
    directTargs[which(directTargs %in% topGenes)]
})

#Now, this is made into a table
allTopTargets <- unique(unlist(topAmongDirectTargets))
topTargetTable <- do.call("cbind", lapply(allTopTargets, function(x){
    unlist(lapply(topAmongDirectTargets, function(y){
        if(x %in% y){
            TRUE
        } else {
            FALSE
        }
    }))
}))

colnames(topTargetTable) <- allTopTargets
topTargetTable
#                   CD27 CD38  IL6R SDC1 CSF2RB
#AE_CSF             TRUE TRUE  TRUE TRUE  FALSE
#Covid_PBMC         TRUE TRUE  TRUE TRUE   TRUE
#Healthy_PBMC       TRUE TRUE  TRUE TRUE   TRUE
#Influensa_PBMC     TRUE TRUE  TRUE TRUE   TRUE
#MS_Ramesh_CSF      TRUE TRUE  TRUE TRUE   TRUE
#MS_Ramesh_PBMC     TRUE TRUE  TRUE TRUE   TRUE
#MS_Shafflick_CSF  FALSE TRUE FALSE TRUE  FALSE
#MS_Shafflick_PBMC  TRUE TRUE  TRUE TRUE   TRUE

#So this shows that the results from the AE analysis is in analogy with the
#resto fo the datasets, with a few exceptions, all coming from the Shafflick CSF
#dataset. 

#Now we are going to try to plot this in a heatmap format. 
sceList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")
#We have to include also the few positive B cells and negative ASC
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
aeSceRed <- aeSce[which(rowData(aeSce)$hgnc_symbol %in% row.names(sceList$AE_CSF)),]
row.names(aeSceRed) <- rowData(aeSceRed)$hgnc_symbol
aeSceOrd <- aeSceRed[match(row.names(sceList$AE_CSF), row.names(aeSceRed)),]
sceList$AE_pos <- aeSceOrd[,which(aeSceOrd$Specific == "TRUE")]
sceList$AE_neg <- aeSceOrd[,which(aeSceOrd$Specific == "FALSE")]
sceList$AE_CSF <- NULL

heatDat <- do.call("cbind", lapply(sceList, function(x){
    locLogCounts <- logcounts(x)[which(rowData(x)$hgnc_symbol %in% directTargs),]
    locASCLogAv <- rowMeans(locLogCounts[,which(x$cellType == "ASC")])
    locBLogAv <- rowMeans(locLogCounts[,which(x$cellType == "B")])
    data.frame("ASC" = locASCLogAv/locASCLogAv, "B" = locBLogAv/locASCLogAv)
}))

heatDatRight <- heatDat[,match(c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                 "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                 "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                 "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B",
                                 "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                 "Covid_PBMC.ASC", "Covid_PBMC.B",
                                 "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                 "Influensa_PBMC.ASC", "Influensa_PBMC.B"),
                               colnames(heatDat))]

library(pheatmap)
pdf("Results/ASC_vs_B/Heatmap_of_drug_targets.pdf", width = 5, height = 6)
pheatmap(t(heatDatRight), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE)
dev.off()

#We also make an alternative version with the downsampled data, to avoid any other
#normalisations
sceList <- readRDS("Data/Comp_to_others/Downsampled_B_and_LGI1-C2_neg_included.rds")
heatDat <- do.call("cbind", lapply(sceList, function(x){
    locLogCounts <- logcounts(x)[which(rownames(x) %in% directTargs),]
    locASCLogAv <- rowMeans(locLogCounts[,which(x$cellType == "ASC")])
    locBLogAv <- rowMeans(locLogCounts[,which(x$cellType == "B")])
    data.frame("ASC" = locASCLogAv, "B" = locBLogAv)
}))

#Here, we divide each column by its max value, to make all variation more visible
heatDatDiv <- as.data.frame(t(apply(heatDat, 1, function(x) x/max(x))))

heatDatRight <- heatDatDiv[,match(c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                 "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                 "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                 "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B",
                                 "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                 "Covid_PBMC.ASC", "Covid_PBMC.B",
                                 "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                 "Influensa_PBMC.ASC", "Influensa_PBMC.B"),
                               colnames(heatDatDiv))]

library(pheatmap)
pdf("Results/ASC_vs_B/Heatmap_of_drug_targets_downsampled.pdf", width = 5, height = 6)
pheatmap(t(heatDatRight), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE)
dev.off()


#Now, we make a supplementary figure that we might not include, namely showing the overlap
#in expression pattern for the nine indirect drug targets that are expressed 
#among the significantly regulated genes, but that are not on a group-level 
#significant. 
indirectTargs <- c("BSG", "CCR5", "IFNAR2", "IL2RB", "MIF",  "PSMA5", "PSMB6", "PSMD8", "TNFRSF17")

heatDat <- do.call("cbind", lapply(sceList, function(x){
    locLogCounts <- logcounts(x)[which(rowData(x)$hgnc_symbol %in% indirectTargs),]
    locASCLogAv <- rowMeans(locLogCounts[,which(x$cellType == "ASC")])
    locBLogAv <- rowMeans(locLogCounts[,which(x$cellType == "B")])
    data.frame("ASC" = locASCLogAv/locASCLogAv, "B" = locBLogAv/locASCLogAv)
}))

heatDat[which(is.nan(heatDat) | is.infinite(heatDat))] <- 0

#A few are not expressed at all, so they will need to be zeroed here. 
heatDat2 <- do.call("cbind", lapply(seq_len(ncol(heatDat)), function(x){
    locCol <- heatDat[,x]
    if(length(which(is.nan(locCol) | is.infinite(locCol))) > 0){
        locCol[which(is.nan(locCol) | is.infinite(locCol))] <- 0
    } 
    locCol
    }))

dimnames(heatDat2) <- dimnames(heatDat )

library(pheatmap)
pdf("Results/ASC_vs_B/Heatmap_of_indirect_drug_targets.pdf", width = 5, height = 6)
pheatmap(t(heatDat2), 
         cluster_rows = FALSE, 
         cluster_cols = FALSE)
dev.off()

