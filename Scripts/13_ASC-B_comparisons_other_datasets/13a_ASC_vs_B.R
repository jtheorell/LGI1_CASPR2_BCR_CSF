#We have already done this now for our cells, but the question is: is there anything
#particular that differs between our ASC and others? We compare therefore to a considerable
#number of examples from the literature
library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(edgeR)

sceList <- readRDS("Data/Comp_to_others/All_sce_common_genes_pre_all.rds")

do.call("rbind", lapply(sceList, function(x) table(x$cellType)))

#.                  ASC     B
#AE_CSF             117   223
#Covid_PBMC        1048 19504
#Healthy_PBMC        35  3754
#Influensa_PBMC      47  6062
#MS_Ramesh_CSF       38   746
#MS_Ramesh_PBMC     100  7184
#MS_Shafflick_CSF   186   425
#MS_Shafflick_PBMC   36  2339

#We also check what the sizes are of the ASC compartment when divided per donor
numDonWithASCn <- do.call("rbind", lapply(c(3:8), function(i){
    unlist(lapply(sceList, function(x){
        length(which(table(x$donor[which(x$cellType == "ASC")]) > i))
    }))
}))
row.names(numDonWithASCn) <- paste0("ASC_n>", c(3:8))
numDonWithASCn
#        AE_CSF Covid_PBMC Healthy_PBMC Influensa_PBMC MS_Ramesh_CSF MS_Ramesh_PBMC MS_Shafflick_CSF MS_Shafflick_PBMC
#ASC_n>3      3         22            4              6             4              8                5                 4
#ASC_n>4      3         19            4              5             2              6                3                 4
#ASC_n>5      3         18            4              3             2              6                3                 4
#ASC_n>6      3         17            3              3             2              6                3                 3
#ASC_n>7      3         16            2              3             2              6                2                 2
#ASC_n>8      3         14            2              3             2              6                2                 2

#This excercise sadly shows that we need to get down to four cells to get 
#sufficient numbers of donors from all datasets. As can also be seen, the other
#datasets also start dropping of soon thereafter. 

#We will now the genes that were picked up as expressed in our spec vs non-spec comparison. 
expressedHGNC <- row.names(readRDS("Results/Specific_vs_not/EdgeR_result.rds"))
sceListRed <- lapply(sceList, function(x) x[which(row.names(x) %in% expressedHGNC),])

#ARe the row names now identical between all these experiments? 
all(unlist(lapply(sceListRed, function(x) identical(row.names(x), row.names(sceListRed[[1]])))))
#TRUE

#Now, we loop this. All documentation is identical to the 8_comparisons_AE/8b_ASC_vs_B_comparisons_AE/1_ASC_vs_B.R script
topResList <- lapply(names(sceListRed), function(x){
    locDat <- sceListRed[[x]]
    current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = counts(locDat))), 
                                    ids=colData(locDat)[,c("cellType", "donor")])
    y <- DGEList(counts(current), samples=colData(current))
    discarded <- current$ncells < 4
    y <- y[,!discarded]
    y <- calcNormFactors(y)
    design <- model.matrix(~factor(cellType), y$samples)
    y <- estimateDisp(y, design)
    fit <- glmQLFit(y, design, robust=TRUE)
    res <- glmQLFTest(fit, coef=ncol(design))
    resMinus <- res
    #We change sign here, to get the ASC and not the B-cell perspective
    resMinus$table$logFC <- -resMinus$table$logFC
    dir.create(paste0("Data/ASC_vs_B/", x))
    saveRDS(resMinus, paste0("Data/ASC_vs_B/", x, "/EdgeR_results.rds"))
    topRows <- which(row.names(locDat) %in% row.names(topTags(resMinus)))
    hgncSymbs <- rowData(locDat)$hgnc_symbol[topRows]
    locRowNames <- row.names(locDat)[topRows]
    hgncSymbsOrd <- hgncSymbs[match(row.names(topTags(resMinus)), locRowNames)]
    topRes <- topTags(resMinus)[[1]]
    row.names(topRes) <- hgncSymbsOrd
    topRes
})

##Any overlaps in the very top between all, or the subgroups?
#topResNameList <- lapply(topResList, row.names)
#names(topResNameList) <- names(sceList)
#
#topGeneCommon <- Reduce(intersect, topResNameList)
##None, but 
## - MZB1, HSP90B1 and FKBP11 is present in all but Healthy, 
## - XBP1 is present in all but Covid
#
#nonMsTopGenesCommon <- Reduce(intersect, topResNameList[-grep("MS", names(topResNameList))])
##None
#
#msTopGenesCommon <- Reduce(intersect, topResNameList[grep("MS", names(topResNameList))])
##"XBP1"    "MZB1"    "FKBP11"  "HSP90B1"
#
##XBP1 is a plasma-cell transcription factor. 
##FKBP11 is a plasma-cell specific antibody-folding catalyst
##MZB1 is called Marginal Zone B And B1 Cell Specific Protein. According to 
##the primary cell atlas at biogps: http://biogps.org/#goto=genereport&id=51237, 
##it is highly ASC specific
#
#csfTopGenesCommon <- Reduce(intersect, topResNameList[grep("CSF", names(topResNameList))])
##"XBP1"    "MZB1"    "FKBP11"  "HSP90B1" "CD38" 
#
##Now, we import all these files, and see what hits that we get that are specific to
##our data
#
#edgeResFileList <- list.files("Data/ASC_vs_B", recursive = TRUE, full.names = TRUE,
#                              pattern = "EdgeR_results.rds")
#
#edgeResList <- lapply(edgeResFileList, readRDS)
#
#names(edgeResList) <- gsub("Data/ASC_vs_B/|/EdgeR_results.rds", 
#                         "", edgeResFileList)
#
##We remove our data here
#aeEdge <- edgeResList[[grep("AE", names(edgeResList))]]
#edgeResList <- edgeResList[-grep("AE", names(edgeResList))]
#
##Now, we identify the hits for all of these. We include 100 on each side, to
##make sure that we do not lose any hits on the thresholds. 
#
#hitGenesUp <- lapply(edgeResList, function(x){
#    locTab <- x$table
#    locSig <- locTab[which(locTab$PValue < 0.01),]
#    row.names(locSig)[which(locSig$logFC > 2)]
#})
#
#hitGenesDown <- lapply(edgeResList, function(x){
#    locTab <- x$table
#    locSig <- locTab[which(locTab$PValue < 0.01),]
#    row.names(locSig)[which(locSig$logFC < -2)]
#})
#
##Now, which of our genes are unique then? 
#rawHits <- read.csv( "Results/ASC_vs_B/AE/Raw_top_results.csv", row.names = 1)
#
##Here, we divide this further into an up and a down compartment
#rawHitsUp <- rawHits[which(rawHits$logFC > 0),]
#rawHitsDown <- rawHits[which(rawHits$logFC < 0),]
#
##Now to the comparisons
#aeSharedUp <- unique(unlist(lapply(hitGenesUp, function(x){
#    x[which(x %in% row.names(rawHitsUp))]
#})))
#length(aeSharedUp)
##75
#
#aeSharedDown <- unique(unlist(lapply(hitGenesDown, function(x){
#    x[which(x %in% row.names(rawHitsDown))]
#})))
#length(aeSharedDown)
##23
#
##And the rest? 
#uniqueAETab <- rawHits[-which(row.names(rawHits) %in% c(aeSharedUp, aeSharedDown)),]
#
#uniqueAETabOrd <- uniqueAETab[order(uniqueAETab$logFC),]
#
#volcDf <- aeEdge$table
#volcDf$Direction <- "None"
#volcDf$Direction[which(row.names(volcDf) %in% row.names(uniqueAETab) & volcDf$logFC > 0)] <- "Up"
#volcDf$Direction[which(row.names(volcDf) %in% row.names(uniqueAETab) & volcDf$logFC < 0)] <- "Down"
#
#volcDf$Direction <- factor(volcDf$Direction, levels = c("Up", "None", "Down"))
#ggplot(data=volcDf, aes(x=logFC, y=-log10(PValue), 
#                        col=Direction)) + 
#    geom_point() + 
#    theme_bw() +
#    scale_color_manual(values = c("red", "grey", "blue")) + xlim(-20, 20) + ylim(0, 8)
#ggsave("Results/ASC_vs_B/AE/Volcano_top_hits_post_filtering.pdf", width = 6, height = 5)
#
#write.csv(uniqueAETabOrd, "Results/ASC_vs_B/AE/Filtered_top_results.csv")
#
##We also want to investigate what percentage of these top genes that were 
##represented in each dataset. This is done here. 
#aeSharedUpFraction <- unlist(lapply(hitGenesUp, function(x){
#    length(x[which(x %in% row.names(rawHitsUp))])/nrow(rawHitsUp)
#}))
#
#aeSharedDownFraction <- unlist(lapply(hitGenesDown, function(x){
#    length(x[which(x %in% row.names(rawHitsDown))])/nrow(rawHitsDown)
#}))
#
#ggDat <- data.frame("Dataset" = c(names(aeSharedUpFraction), names(aeSharedDownFraction)),
#                    "Direction" = c(rep("Up", length(aeSharedUpFraction)),
#                                    rep("Down", length(aeSharedDownFraction))),
#                    "Fraction" = c(unname(aeSharedUpFraction), unname(aeSharedDownFraction)))
#
#ggplot(data=ggDat, aes(x=Dataset, y=Fraction, fill=Direction)) +
#    geom_bar(stat="identity", position=position_dodge()) + scale_fill_manual(values = c("blue", "red")) +
#    theme_bw() + scale_y_continuous(expand = c(0, 0)) + ylab("Fraction of hits among AE top hits")
#ggsave("Results/ASC_vs_B/Shared_top_hit_fraction.pdf")
