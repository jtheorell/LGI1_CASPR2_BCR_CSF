library(SingleCellExperiment)
library(Biobase)
library(ggplot2)
library(scran)
library(scater)
library(pheatmap)
#Here, we are going to subtract the CD8 signals from the ASC signals from 
#all the external experiments, and then integrat with our own data. The hypothesis
#is that a substantial fraction of the imbalance noted between the upregulated and
#downregulated genes is due simply to our dataset being created with fresh cells
#and with deep coverage. Therefore, it is hypothetically so that the overrepresented genes
#to a high degree represent normal metabolism (it seemed so in the early analyses), and 
#therefore it is not unlikely that these "upregulated" genes are shared between unrelated cell
#types. This is the rationale for doing the following. 

#The first step is to import and integrate the 100 files. 
fileList <- list.files("Results/Specific_vs_others/Bootstrap_DE", full.names = TRUE)
ASCDatLong <- lapply(fileList, readRDS)

#I have previously checked that the row names are identical between all frames, so
#we can here run a simplified version. 
#We are also changing the signs here, so that a fold change increase means that our
#cells have an increase compared to the controls, and so on. 

integrDat <- lapply(names(ASCDatLong[[1]]), function(x){
    locDat <- lapply(ASCDatLong, "[[", which(names(ASCDatLong[[1]]) == x))
    locNames <- row.names(locDat[[1]])
    locResList <- lapply(locDat, function(y){
        locTab <- y$table
        list("logFC" = locTab$logFC, "p.val" = locTab$PValue,
             "F" = locTab$F)
    })
    medFC <- -rowMedians(do.call("cbind", lapply(locResList, "[[", 1)))
    medP <- rowMedians(do.call("cbind", lapply(locResList, "[[", 2)))
    medF <- rowMedians(do.call("cbind", lapply(locResList, "[[", 3)))
    data.frame("logFC" = medFC, 
               "p.value" = medP, 
               "F_stat" = medF,
               row.names = row.names(locDat[[1]]))
})
names(integrDat) <- names(ASCDatLong[[1]])

#Here, we make a median assay for the above. We do not include our data here, as
#it does not correlate with the rest. 
medFCDat <- rowMedians(do.call("cbind", lapply(integrDat, function(x) x$logFC)))
medPDat <- rowMedians(do.call("cbind", lapply(integrDat, function(x) x$p.value)))
medFDat <- rowMedians(do.call("cbind", lapply(integrDat, function(x) x$F_stat)))

integrDat$Median <- data.frame("logFC" = medFCDat, 
                               "p.value" = medPDat, 
                               "F_stat" = medFDat,
                               row.names = row.names(integrDat[[1]]))

saveRDS(integrDat, "Results/Specific_vs_others/All_DE_integrated_raw_downsamp.rds")
#integrDat <- readRDS("Results/Specific_vs_others/All_DE_integrated_raw_downsamp.rds")

############
#Investigation of highly expressed genes. 
#Previously, it has seemed like there was a strong correlation between the expression
#level and the likelihood of finding hits, so that almost all hits were among genes expressed
#by all cells. Therefore, we here investigate if that is still the case. 

allDatList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

#First, all genes not included are filtered out
allDatListRed <- lapply(allDatList, function(x){
    x[which(row.names(x) %in% row.names(integrDat$Covid_PBMC))]
})
#Are all identical? 
all(unlist(lapply(allDatListRed, function(x) identical(row.names(x), row.names(allDatListRed$AE_CSF)))))

numPosMat <- do.call("cbind", lapply(allDatListRed, function(x){
    locCounts <- counts(x)
    locNumPosVec <- apply(locCounts, 1, function(y) length(y[y>0]))
    locFracPosVec <- locNumPosVec/ncol(x)
}))
minValVec <- rowMin(numPosMat)

#Now, we plot a very interesting metric, namely the percentage of expression as 
#a function of the Cohen value. We will also add the ribosomal genes as a color

#As an add-on to this, we clear out the ribosomal genes
library(biomaRt)
ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl") #uses human ensembl annotations
#GO:0005840 = ribosome
#GO:0003735 = ribosomal protein
#GO:0006412 = translation
ribosomal<- getBM(attributes=c('hgnc_symbol', 'go_id'),
                  filters = 'go', values = c('GO:0005840', "GO:0003735", "GO:0006412"), mart = ensembl)
#This catches 6748 genes. 

colVec <- rep("black", length(minValVec))
colVec[which(row.names(integrDat$Median) %in% ribosomal$hgnc_symbol)] <- "blue"
pdf("Results/Specific_vs_others/DE_logFC_vs_minimal_percent_expression_blu_is_ribosomal.pdf")
plot(minValVec, integrDat$Median$logFC, xlab = "Lowest fraction of expressing cells, %",
     ylab = "Median Cohen-corrected logFC", pch = 16, col = colVec)
dev.off()
#And now, we remove the ribosomally-associated genes here. 
integrDatNoRib <- lapply(integrDat, function(x){
    x[-which(row.names(x) %in% ribosomal$hgnc_symbol),]
})

#Still 8932 genes. 

saveRDS(integrDatNoRib, "Results/Specific_vs_others/All_DE_integrated_downsamp_reduced.rds")
#integrDatNoRib <- readRDS("Results/Specific_vs_others/All_DE_integrated_downsamp_reduced.rds")

#############
#COMPARISONS
#############
#We start by identifying the most up-and downregulated genes in our spec-vs-not and 
#these spec-vs-other datasets
specVsNot <- read.csv("Results/Specific_vs_not/Top_results.csv", row.names = 1)
specVsNot$p.value <- specVsNot$PValue
integrDatNoRib$Spec_vs_not <- specVsNot

upGeneList <- lapply(integrDatNoRib[-which(names(integrDatNoRib) == "Median")], function(x){
    row.names(x)[which(x$p.value < 0.05 & x$logFC > 2)]
})

#Are there any shared by all?

upGeneCommon <- Reduce(intersect, upGeneList)
upGeneCommon
#Nope. 

#Spec vs PBMC?
upGenePbmcAll <- Reduce(intersect, upGeneList[-grep("CSF", names(upGeneList))])
upGenePbmc <- upGenePbmcAll[-which(upGenePbmcAll %in% upGeneCommon)]
upGenePbmc
#Nope

#If the MS cells are removed? 
upGeneNonMsAll <- Reduce(intersect, upGeneList[-grep("MS", names(upGeneList))])
upGeneNonMs <- upGeneNonMsAll[-which(upGeneNonMsAll %in% c(upGeneCommon, upGenePbmcAll))]
upGeneNonMs
#This adds:
#None

#And what if we instead remove the non-MS?
upGeneMsAll <- Reduce(intersect, upGeneList[grep("MS|Spec", names(upGeneList))])
upGeneMs <- upGeneMsAll[-which(upGeneMsAll %in% upGeneCommon)]
upGeneMs
#None

#We also make the further separation into the MS PBMC and MS CSF compartments
upGeneMsCsfAll <- Reduce(intersect, upGeneList[grep("CSF|Spec", names(upGeneList))])
upGeneMsCsf <- upGeneMsCsfAll[-which(upGeneMsCsfAll %in% upGeneMsAll)]
upGeneMsCsf
#None

#And the PBMC?
upGeneMsPbmcAll <- Reduce(intersect, upGeneList[grep("MS_.+_PBMC|Spec", names(upGeneList))])
upGeneMsPbmc <- upGeneMsPbmcAll[-which(upGeneMsPbmcAll %in% c(upGeneMsAll, upGenePbmcAll))]
#None

#Now: downregulated genes?
downGeneList <-  lapply(integrDatNoRib[-which(names(integrDatNoRib) == "Median")], function(x){
    row.names(x)[which(x$FDR.p < 0.05 & x$logFC < -2)]
})

downGeneCommon <- Reduce(intersect, downGeneList)
downGeneCommon
#None

#The PBMC?
downGenePbmcAll <- Reduce(intersect, downGeneList[-grep("CSF", names(downGeneList))])
downGenePbmc <- downGenePbmcAll[-which(downGenePbmcAll %in% downGeneCommon)]
downGenePbmc
#None

#Now, if MS is excluded?
downGeneNonMsAll <- Reduce(intersect, downGeneList[-grep("MS", names(upGeneList))])
downGeneNonMs <- downGeneNonMsAll[-which(downGeneNonMsAll %in% downGeneCommon)]
downGeneNonMs
#None

#Now, what about the MS?
downGeneMsAll <- Reduce(intersect, downGeneList[grep("MS|Spec", names(upGeneList))])
downGeneMs <- downGeneMsAll[-which(downGeneMsAll %in% downGeneCommon)]
downGeneMs
#None here either. 

#And when the CSF and PBMC compartments among the CSF are separated? 
downGeneMsCsfAll <- Reduce(intersect, downGeneList[grep("CSF|Spec", names(downGeneList))])
downGeneMsCsf <- downGeneMsCsfAll[-which(downGeneMsCsfAll %in% downGeneMsAll)]
downGeneMsCsf
#None

#And the PBMC?
downGeneMsPbmcAll <- Reduce(intersect, downGeneList[grep("MS_.+_PBMC|Spec", names(downGeneList))])
downGeneMsPbmc <- downGeneMsPbmcAll[-which(downGeneMsPbmcAll %in% c(downGeneMsAll, downGenePbmcAll))]
downGeneMsPbmc
#None

#So, summing this up: 
#All: "DPEP1"
#PBMC: "TCN2"
#Non-MS: "HDAC6" "PLEKHO2"
#MS-CSF: PLCD1"
#MS-PBMC: "H6PD" "KLC4"

#Old, when we had hits. 
#And we create a table of this. 
#hitList <- c("DPEP1", "TCN2", "HDAC6", "PLEKHO2", "PLCD1", "H6PD", "KLC4")
#
#hitPattern <- do.call("rbind", lapply(upGeneList, function(x){
#    hitList %in% x
#}))
#colnames(hitPattern) <- hitList
#row.names(hitPattern) <- paste0(row.names(hitPattern), ".ASC")

#############
#VOLCANO
#############
#THis is now plotted as a volcano plot
volcDf <- readRDS("Results/Specific_vs_not/EdgeR_result.rds")$table
volcDf$p.value <- volcDf$PValue
volcDf$colors <- "None"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneCommon, downGeneCommon))] <- "Common"
volcDf$colors[which(row.names(volcDf) %in% c(upGenePbmc, downGenePbmc))] <- "PBMC"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneNonMs, downGeneNonMs))] <- "Non-MS"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneMs, downGeneMs))] <- "MS"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneMsPbmc, downGeneMsPbmc))] <- "MS-PBMC" 
volcDf$colors[which(row.names(volcDf) %in% c(upGeneMsCsf, downGeneMsCsf))] <- "MS-CSF"
volcDf$hit <- FALSE
volcDf$hit[-which(volcDf$colors == "None")] <- TRUE
volcDf$colors <- factor(volcDf$colors, levels = c("None", "Common", "PBMC", "Non-MS", "MS", "MS-PBMC", "MS-CSF"))
ggplot(data=volcDf, aes(x=logFC, y=-log10(p.value), 
                        col=colors, size = hit)) + 
    geom_point() + 
    theme_bw() + scale_size_manual(values =c(1,3)) +
    scale_color_manual(values = c("grey", "#8B0AA5FF", "#ED7953FF", "#FDB32FFF",  "#DB5C68FF", "black")) +
    xlim(-15, 15) + ylim(0, 3)
ggsave("Results/Specific_vs_others/Volcano_top_hits_downsamp.pdf", width = 6, height = 5)

###############
#PLOTTING OF HITS
###############
#Here we will first get the data for the downsampled datasets, and then also integrate
#the B-cells for all the datasets
###############
heatmapSceList <- readRDS("Data/Comp_to_others/Downsampled_B_and_LGI1-C2_neg_included.rds")

#Next, we will restrict the rows to the hits
heatmapSceListRed <- lapply(heatmapSceList, function(x){
    x[which(row.names(x) %in% hitList),]
    x[match(hitList, row.names(x)),]
})

heatDat <- do.call("cbind", lapply(heatmapSceListRed, function(x){
    locBLogAv <- rowMeans(logcounts(x[,which(x$cellType == "B")]))
    locAscLogAv <- rowMeans(logcounts(x[,which(x$cellType == "ASC")]))
    data.frame("B" = locBLogAv, 
               "ASC" = locAscLogAv)
}))


#Now, we normalise this to the highest value in each row. 
heatDatFrac <- t(apply(heatDat, 1, function(x) x/max(x)))
#Now, we reorder the columns as we would like them
heatDatFracOrd <- heatDatFrac[,match(c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                       "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                       "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                       "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B", 
                                       "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                       "Covid_PBMC.ASC", "Covid_PBMC.B",
                                       "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                       "Influensa_PBMC.ASC", "Influensa_PBMC.B"), 
                                     colnames(heatDatFrac))]

#In addition tot he actual heatmap, we will also include data on if the value
#was identified as significant in the comparison at hand. Red border color will
#indicate this, whereas grey will not. 
edgeResFileList <- list.files("Data/ASC_vs_B", recursive = TRUE, full.names = TRUE,
                              pattern = "EdgeR_results.rds")

ascVsBResList <- lapply(edgeResFileList, readRDS)

names(ascVsBResList) <- gsub("Data/ASC_vs_B/|/EdgeR_results.rds", 
                             "", edgeResFileList)

ascVsBTabList <- lapply(ascVsBResList, function(x) x$table)

#Now, we identify the hits in each of these, with the same criterion as above. 
upGeneListAscVsB <- lapply(ascVsBTabList, function(x){
    row.names(x)[which(x$PValue < 0.05 & x$logFC > 2)]
})

#Now, are any of these hits present among these, and if so, how does it look?

#Which have which? 
hitPatternB <- do.call("rbind", lapply(upGeneListAscVsB, function(x){
    hitList %in% x
}))
colnames(hitPatternB) <- hitList
row.names(hitPatternB) <- paste0(row.names(hitPatternB), ".B")

fullHitPattern <- rbind(hitPattern, hitPatternB)
row.names(fullHitPattern)[which(row.names(fullHitPattern) == "AE_CSF.B")] <- "AE_neg.B"
fullHitPatternOrd <- fullHitPattern[match(c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC", "AE_neg.B", 
                                            "MS_Ramesh_CSF.ASC", "MS_Ramesh_CSF.B",
                                            "MS_Ramesh_PBMC.ASC", "MS_Ramesh_PBMC.B",
                                            "MS_Shafflick_CSF.ASC", "MS_Shafflick_CSF.B", 
                                            "MS_Shafflick_PBMC.ASC", "MS_Shafflick_PBMC.B", 
                                            "Covid_PBMC.ASC", "Covid_PBMC.B",
                                            "Healthy_PBMC.ASC", "Healthy_PBMC.B", 
                                            "Influensa_PBMC.ASC", "Influensa_PBMC.B"), row.names(fullHitPattern)),]
row.names(fullHitPatternOrd)[1:3] <- c("AE_pos.ASC", "AE_pos.B", "AE_neg.ASC")
fullHitPatternOrd

#                      DPEP1  TCN2 HDAC6 PLEKHO2 PLCD1  H6PD  KLC4
#AE_pos.ASC               NA    NA    NA      NA    NA    NA    NA
#AE_pos.B                 NA    NA    NA      NA    NA    NA    NA
#AE_neg.ASC               NA    NA    NA      NA    NA    NA    NA
#AE_neg.B               TRUE  TRUE FALSE   FALSE FALSE FALSE FALSE
#MS_Ramesh_CSF.ASC      TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE
#MS_Ramesh_CSF.B        TRUE FALSE FALSE   FALSE FALSE FALSE FALSE
#MS_Ramesh_PBMC.ASC     TRUE  TRUE FALSE   FALSE FALSE  TRUE  TRUE
#MS_Ramesh_PBMC.B       TRUE FALSE FALSE   FALSE FALSE FALSE FALSE
#MS_Shafflick_CSF.ASC   TRUE FALSE FALSE   FALSE  TRUE FALSE  TRUE
#MS_Shafflick_CSF.B    FALSE FALSE FALSE   FALSE FALSE FALSE FALSE
#MS_Shafflick_PBMC.ASC  TRUE  TRUE  TRUE    TRUE  TRUE  TRUE  TRUE
#MS_Shafflick_PBMC.B    TRUE FALSE FALSE   FALSE FALSE FALSE FALSE
#Covid_PBMC.ASC         TRUE  TRUE  TRUE    TRUE FALSE  TRUE FALSE
#Covid_PBMC.B           TRUE  TRUE FALSE   FALSE FALSE FALSE FALSE
#Healthy_PBMC.ASC       TRUE  TRUE  TRUE    TRUE FALSE FALSE FALSE
#Healthy_PBMC.B         TRUE FALSE FALSE   FALSE  TRUE FALSE FALSE
#Influensa_PBMC.ASC     TRUE  TRUE  TRUE    TRUE  TRUE  TRUE FALSE
#Influensa_PBMC.B       TRUE FALSE FALSE   FALSE FALSE FALSE FALSE


#This information is now used for the heatmap
fullHitPatternVec <- as.vector(fullHitPatternOrd)
fullHitPatternCols <- sapply(as.character(fullHitPatternVec), switch, "TRUE" = "red", 
                             "FALSE" = "blue", "NA" = "grey60")

pdf("Results/Specific_vs_others/Heatmap_of_hits.pdf", width = 5, height = 6)
pheatmap(t(heatDatFracOrd), 
         border_color = fullHitPatternCols,
         cluster_rows = FALSE, 
         cluster_cols = FALSE)
dev.off()



