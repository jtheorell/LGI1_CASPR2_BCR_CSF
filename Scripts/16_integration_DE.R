library(SingleCellExperiment)
library(Biobase)
library(ggplot2)
library(scran)
library(scater)
#Here, we are going to subtract the CD8 signals from the ASC signals from 
#all the external experiments, and then integrat with our own data. The hypothesis
#is that a substantial fraction of the imbalance noted between the upregulated and
#downregulated genes is due simply to our dataset being created with fresh cells
#and with deep coverage. Therefore, it is hypothetically so that the overrepresented genes
#to a high degree represent normal metabolism (it seemed so in the early analyses), and 
#therefore it is not unlikely that these "upregulated" genes are shared between unrelated cell
#types. This is the rationale for doing the following. 

#The first step is to import and integrate the 100 files. 
fileList <- list.files("Results/Comp_to_others/Bootstrap_DE", full.names = TRUE)
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

saveRDS(integrDat, "Results/Comp_to_others/All_DE_integrated_raw_downsamp.rds")
#integrDat <- readRDS("Results/Comp_to_others/All_DE_integrated_raw_downsamp.rds")

############
#Investigation of highly expressed genes. 
#Previously, it has seemed like there was a strong correlation between the expression
#level and the likelihood of finding hits, so that almost all hits were among genes expressed
#by all cells. Therefore, we here investigate if that is stil the case. 

allDatList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

    
#Here, we remove the non-specific cells, as that dataset is the smallest and of lowest
#importance, as we find no differences with it. 
numPosMat <- do.call("cbind", lapply(allDatList, function(x){
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
pdf("Results/Comp_to_others/DE_logFC_vs_minimal_percent_expression_blu_is_ribosomal.pdf")
plot(minValVec, integrDat$Median$logFC, xlab = "Lowest fraction of expressing cells, %",
     ylab = "Median Cohen-corrected logFC", pch = 16, col = colVec)
dev.off()
#And now, we remove the ribosomally-associated genes here. 
integrDatNoRib <- lapply(integrDat, function(x){
    x[-which(row.names(x) %in% ribosomal$hgnc_symbol),]
})

#Still 7016 genes. 

#Now, we calculate the FDR-corrected p-values
integrDatNoRib <- lapply(integrDatNoRib, function(x){
    x$FDR.p <- p.adjust(x$p.value, method = "fdr")
    x
})
#############
#COMPARISONS
#############
#We start by identifying the most up-and downregulated genes. 

upGeneList <- lapply(integrDatNoRib[-which(names(integrDatNoRib) == "Median")], function(x){
    row.names(x)[which(x$FDR.p < 0.01 & x$logFC > 3)]
})

#Are there any shared by all?

upGeneCommon <- Reduce(intersect, upGeneList)

#Yes: 
#"DDIT4"
#That are shared.

#If the MS-CSF cells are removed? 

upGeneNonMsAll <- Reduce(intersect, upGeneList[-which(names(upGeneList) %in% 
                                                        c("MS_Ramesh", "MS_Schafflick"))])
upGeneNonMs <- upGeneNonMsAll[-which(upGeneNonMsAll %in% upGeneCommon)]
#This adds:
#"DHCR24" "FADS1"  "IL16"   "IL21R"  "INSIG1" "PHGDH"  "SCD"    "SPAG4" 

#And what if we instead remove the non-MS?
upGeneMsAll <- Reduce(intersect, upGeneList[which(names(upGeneList) %in% 
                                                        c("MS_Ramesh", "MS_Schafflick"))])
upGeneMs <- upGeneMsAll[-which(upGeneMsAll %in% upGeneCommon)]

#"ADM"      "ADORA2A"  "BCKDHA"   "BOLA2B"   "BSCL2"    "C19orf48" "COPG1"    "DDX39A"   "EWSR1"   
#"IFI30"    "IGLL5"    "LIME1"    "MAN2B1"   "MCM2"     "MCM3"     "MCM4"     "MCM5"     "MCM7"    
#"NOMO2"    "PNP"      "RRM2"     "SEC13"    "SELPLG"   "SLX1A"    "STT3A"    "TBL3"     "TK1"     
#"TM9SF1"   "TPI1"     "TXNDC5"   "UQCRC1"   "USP11"   

#Now: downregulated genes?
downGeneList <-  lapply(integrDatNoRib[-which(names(integrDatNoRib) == "Median")], function(x){
    row.names(x)[which(x$FDR.p < 0.01 & x$logFC < -3)]
})

downGeneCommon <- Reduce(intersect, downGeneList)

# "HLA-B" "HLA-C" "PSMB8"

#Now, if MS is excluded?
downGeneNonMsAll <- Reduce(intersect, downGeneList[-which(names(downGeneList) %in% 
                                                        c("MS_Ramesh", "MS_Schafflick", "Specific"))])
downGeneNonMs <- downGeneNonMsAll[-which(downGeneNonMsAll %in% downGeneCommon)]
#"NTAN1"  "PLAC8"  "SMCHD1" "VPS28" 

#Now, what about the MS?
downGeneMsAll <- Reduce(intersect, downGeneList[which(names(downGeneList) %in% 
                                                            c("MS_Ramesh", "MS_Schafflick"))])
downGeneMs <- downGeneMsAll[-which(downGeneMsAll %in% downGeneCommon)]

#"C12orf75" "CHMP4B"   "HLA-DRA"  "HLA-DRB1" "PPDPF"    "TGFB1"    "TPT1"     
#Are unique

#############
#VOLCANO
#############
#THis is now plotted as a volcano plot
volcDf <- integrDatNoRib$Median
volcDf$colors <- "None"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneCommon, downGeneCommon))] <- "Common"
volcDf$hit <- FALSE
volcDf$hit[-which(volcDf$colors == "None")] <- TRUE
ggplot(data=volcDf, aes(x=logFC, y=-log10(FDR.p), 
                        col=colors, size = hit)) + 
    geom_point() + 
    theme_bw() + scale_size_manual(values =c(1,3)) +
    scale_color_manual(values = c("darkviolet", "grey"))
ggsave("Results/Comp_to_others/Volcano_top_common_hits_downsamp.pdf", width = 6, height = 5)

#And now we make an extended version with the subsets
volcDf$colors[which(row.names(volcDf) %in% c(upGeneNonMs, downGeneNonMs))] <- "Non-MS"
volcDf$colors[which(row.names(volcDf) %in% c(upGeneMs, downGeneMs))] <- "MS"
volcDf$hit <- FALSE
volcDf$hit[-which(volcDf$colors == "None")] <- TRUE
ggplot(data=volcDf, aes(x=logFC, y=-log10(FDR.p), 
                               col=colors, size = hit)) + 
    geom_point() + 
    theme_bw() + scale_size_manual(values =c(1,3)) +
    scale_color_manual(values = c("darkviolet", "red", "blue", "grey"))
ggsave("Results/Comp_to_others/Volcano_top_hits_all_downsamp.pdf", width = 6, height = 5)

###############
#PLOTTING OF HITS
###############
rawSceList <- readRDS("Data/Comp_to_others/Downsampled_specificity_included_boot_100.rds")

#First, all the common
allUp <- c(upGeneCommon, upGeneNonMs, upGeneMs)
allDown <- c(downGeneCommon, downGeneNonMs, downGeneMs)

upAndDown <- c(allUp, allDown)

allLogCounts <- do.call("cbind", lapply(rawSceList, logcounts))
allSce <- SingleCellExperiment(assays = list("logcounts" = allLogCounts),
                                     rowData = rowData(rawSceList$Influensa))
allSce <- allSce[which(row.names(allSce) %in% unique(c(upAndDown))),]
allSce$group <- unlist(lapply(rawSceList, function(x) x$group))

#Now visualisation. 

plotExpression(allSce, features=c(upGeneCommon, downGeneCommon), 
               x="group", colour_by="group", ncol = 5)
dev.copy(pdf,'Results/Comp_to_others/Common_hit_genes_downsamp.pdf', height = 4, width = 6)
dev.off()

plotExpression(allSce, features=c(upGeneNonMs, downGeneNonMs), 
               x="group", colour_by="group", ncol = 5)
dev.copy(pdf,'Results/Comp_to_others/Non-Ms_intersecting_hit_genes_downsamp.pdf', height = 12, width = 12)
dev.off()

plotExpression(allSce, features=c(upGeneMs, downGeneMs), 
               x="group", colour_by="group", ncol = 5)
dev.copy(pdf,'Results/Comp_to_others/Ms_intersecting_hit_genes_downsamp.pdf', height = 27, width = 12)
dev.off()

#So, there are essentially two really clear hits: HLA-B and HLA-C. Are these present
#in the CD8s? Yes. 

#And here these genes are saved
allUp <- c(upGeneCommon, upGeneNonMs, upGeneMs)
allDown <- c(downGeneCommon, downGeneNonMs, downGeneMs)

upAndDownCommon <- data.frame("Gene" = c(upGeneCommon, downGeneCommon),
                              "Direction" = c(rep("Up", length(upGeneCommon)), 
                                              rep("Down", length(downGeneCommon))))
upAndDownCommon$Data <- "Common"

upAndDownNonMs <- data.frame("Gene" = c( upGeneNonMs, downGeneNonMs),
                              "Direction" = c(rep("Up", length( upGeneNonMs)), 
                                              rep("Down", length(downGeneNonMs))))
upAndDownNonMs$Data <- "NonMS"

upAndDownMs <- data.frame("Gene" = c(upGeneMs, downGeneMs),
                              "Direction" = c(rep("Up", length(upGeneMs)), 
                                              rep("Down", length(downGeneMs))))
upAndDownMs$Data <- "MS"

upAndDownDf <- rbind(upAndDownCommon, upAndDownNonMs, upAndDownMs)

write.csv(upAndDownDf, "Results/Comp_to_others/Top_DE_hits_downsamp.csv", row.names = FALSE)

#And here, we also export a list containing both the data and the names for all
#datasets.
saveRDS(integrDatNoRib, "Results/Comp_to_others/All_DE_integrated_downsamp.rds")

