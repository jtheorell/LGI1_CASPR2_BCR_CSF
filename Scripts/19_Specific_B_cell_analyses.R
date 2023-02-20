#Here, we are going to investigate those few but intriguing binder B-cells. 

fullBCRDb <- read.csv("Data/BCR_database_versions/4_IMGT_gapped_db_complete_post_clonality.csv")
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#Now, we import information also about the non-IgG, as we are now digging deeper
#into the B-cell compartment. 
fullBCRDbH <- fullBCRDb[which(fullBCRDb$LOCUS == "H" & fullBCRDb$CELL %in% colnames(aeSce)),
                        c("CELL", "V_CALL","JUNCTION",
                          "JUNCTION_LENGTH","ISOTYPE",
                          "HamClone","light_type",
                          "Clonal","All_mutations",
                          "Non_silent_mutations")]

identical(colnames(colData(aeSce))[which(colnames(colData(aeSce)) %in% c("V_CALL","JUNCTION",
                                                                         "JUNCTION_LENGTH",
                                                                         "ISOTYPE","HamClone",
                                                                         "light_type","Clonal",
                                                                         "All_mutations",
                                                                         "Non_silent_mutations"))], 
          colnames(fullBCRDbH)[-1])
#TRUE, so we will be able to transfer this in blocks
missingDf <- as.data.frame(matrix(NA, ncol(aeSce)-nrow(fullBCRDbH), 
                     ncol(fullBCRDbH)))
colnames(missingDf) <- colnames(fullBCRDbH)
missingDf$CELL <- colnames(aeSce)[-which(colnames(aeSce) %in% fullBCRDbH$CELL)]

BCRPlusEmpty <- rbind(fullBCRDbH, missingDf)
BCR_ordered <- BCRPlusEmpty[match(colnames(aeSce), BCRPlusEmpty$CELL),]

identical(BCR_ordered$CELL, colnames(aeSce))
#TRUE
BCR_ordered$CELL <- NULL

colData(aeSce)[,c("V_CALL","JUNCTION",
                  "JUNCTION_LENGTH",
                  "ISOTYPE","HamClone",
                  "light_type","Clonal",
                  "All_mutations",
                  "Non_silent_mutations")] <- BCR_ordered

#Now start the real analyses

aeSceB <- aeSce[,which(aeSce$cellType == "B")]


#Now, we are going to pick out the ones that 
aeSceBSpecKnown <- aeSceB[,which(aeSceB$Specific != "Not_tested")]

table(paste0(aeSceBSpecKnown$Specific, "_", aeSceBSpecKnown$donor), 
      aeSceBSpecKnown$Clonal, useNA = "ifany")
#                 FALSE TRUE
#FALSE_1166          3    1
#FALSE_1227          2    0
#FALSE_1284          6    0
#TRUE_1166           0    2
#TRUE_1284           1    2

#So much overrepresented is clonality and the cells here are luckily distributed among
#two donors and two disorders. 

############
#sPLS-DA
###########
#Now, if this does in fact mean something: which genes define these cells?
#We of course only look at the ones with a true, known specificity. 
library(mixOmics)

#First, we select genes with meaningful variance
aeSceBSpecKnownRes <- modelGeneVarWithSpikes(aeSceBSpecKnown, "ERCC")
hvg.csf.var <- getTopHVGs(aeSceBSpecKnownRes, n=1000)
str(hvg.csf.var)

redAESceB <- aeSceBSpecKnown[which(row.names(aeSceBSpecKnown) %in% hvg.csf.var),]
row.names(redAESceB) <- rowData(redAESceB)$hgnc_symbol

plsRes <- splsda(t(logcounts(redAESceB)), 
                 factor(redAESceB$Specific), ncomp = 1, keepX = 10)

signLoads <- plsRes$loadings$X[which(abs(plsRes$loadings$X[,1]) > 0),]
signLoadsOrd <- data.frame(signLoads[order(abs(signLoads), decreasing = TRUE)])
colnames(signLoadsOrd) <- "loading"
round(signLoadsOrd, 3)
#        loading
#TPST2     0.604
#DOK3      0.579
#JCHAIN    0.352
#XBP1      0.295
#ZFP36L1  -0.189
#STIP1     0.169
#MRPL55    0.133
#TSPO      0.084
#SNX17     0.011
#ANXA6     0.001

signLoadsOrd$protName <- c("Tyrosylprotein Sulfotransferase 2",
                           "Docking protein 3",
                           "Joining chain of multimeric IgA and IgM",
                           "X-box binding protein 1",
                           "ZFP36 ring finger protein like 1",
                           "Stress induced phosphoprotein 1",
                           "Mitochondrial ribosomal protein L55",
                           "Translocator protein",
                           "Sorting nexin 17",
                           "Annexin A6")

signLoadsOrd$minMonacoASCBRatio <- c(277.6/70.3, 
                                     94.6/125.8, 
                                     130626.8/6942.2, 
                                     2875/131.9, 
                                     32.5/250.5,
                                     182.6/154.2,
                                     78.9/50.7,
                                     70.7/71.3,
                                     73.1/80.4,
                                     780.3/665.9)

#All in all, these cells show ASC traits, despite not having an ASC surface phenotype
library(ggplot2)

plotDat <-  data.frame("Specific" = aeSce$Specific, "cellType" = aeSce$cellType,
                       "CD27RNA" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "CD27"),],
                       t(normcounts(altExp(aeSce, "flowData"))),
                       "XBP1" = logcounts(aeSce)[which(rowData(aeSce)$hgnc_symbol == "XBP1"),])
#Now, we want the specifics to be on top
plotDat$Specific <- factor(plotDat$Specific, levels = c("Not_tested", "FALSE", "TRUE"))

plotDatOrdered <- plotDat[order(plotDat$Specific),]

dir.create("Results/Specific_B_cell_analysis")
ggplot(plotDatOrdered, aes(x = CD27RNA, y = CD38, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/CD27RNA_vs_CD38.pdf")

ggplot(plotDatOrdered, aes(x = CD27, y = CD38, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/CD27_vs_CD38.pdf")
#So this shows that the cells are essentially ASC. We now have to check what happens
#if the protein and mRNA are combined, as the CD27 protein clearly is very poor.
CD27combo <- sapply(seq_along(plotDatOrdered$CD27RNA), function(x){
    max(c(plotDatOrdered$CD27RNA[x], plotDatOrdered$CD27[x]))
})
plotDatOrdered$CD27Combo <- CD27combo
ggplot(plotDatOrdered, aes(x = CD27Combo, y = CD38, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/CD27RNAplusProt_vs_CD38.pdf")
#This does not look any different. 
ggplot(plotDatOrdered, aes(x = CD27RNA, y = CD20, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/CD27RNA_vs_CD20.pdf")

ggplot(plotDatOrdered, aes(x = CD138, y = CD20, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/CD138_vs_CD20.pdf")

ggplot(plotDatOrdered, aes(x = XBP1, y = CD20, color = Specific, size = cellType)) + 
    geom_point() + theme_bw() + theme(aspect.ratio=1) + scale_color_manual(values = c("grey","black", "orange"))
ggsave("Results/Specific_B_cell_analysis/XBP1_vs_CD20.pdf")

#Now, what it turns out, all-in-all, is that the cells are ASC with CD20 expression
#And lacking CD138. 















#Now, how good is this model for prediction purposes? 

aeSceBRed <- aeSceB[which(row.names(aeSceB) %in% hvg.csf.var),]
row.names(aeSceBRed) <- rowData(aeSceBRed)$hgnc_symbol
plsPred <- predict(plsRes, t(logcounts(aeSceBRed)))
table(plsPred$MajorityVote$centroids.dist, aeSceBRed$Specific)
#      FALSE Not_tested TRUE
#FALSE    12        199    0
#TRUE      0         12    5

#So the model seems to give meaningful results. However, the cells picked up
#in the not_tested group are comparatively few. But if these are investigated, 
#we can get an educated guess on if they seem to make sense
aeSceBRed$specPred <- plsPred$MajorityVote$centroids.dist
notTestSce <- aeSceBRed[,which(aeSceBRed$Specific == "Not_tested")]

table(notTestSce$specPred, notTestSce$ISOTYPE, useNA = "ifany")

#      IGHA1 IGHA2 IGHD IGHG1 IGHG2 IGHG3 IGHG4 IGHM None Unknown <NA>
#FALSE    33     4    3    39     5     7     7   53   28       1   19
#TRUE      0     1    0     4     0     0     0    4    0       0    3

#This does not necessarily strengthen the argument here. And what about clonality?
table(paste0("SpecPred_", notTestSce$specPred), notTestSce$Clonal, useNA = "ifany")

#               FALSE TRUE <NA>
#SpecPred_FALSE   172    8   19
#SpecPred_TRUE      9    0    3

#Not this either - cells predicted to bind are not IgG4 and not clonal. 




















#Now, if we were to run a k-means thing on all the b-cells, using the true and false
#as reference, what would we pick up, i.e. how uncommon is the phenotype?

glmCentSpec <- colMeans(reducedDims(aeSceBSpecKnown[,which(aeSceBSpecKnown$Specific == "TRUE")])[[1]])

glmCentNonSpec <- colMeans(reducedDims(aeSceBSpecKnown[,which(aeSceBSpecKnown$Specific == "FALSE")])[[1]])

#Now, we cluster all, to see how well these known ones separate and also how the others divide
library(FNN)
phenoClust <- knnx.index(t(data.frame("Spec" = glmCentSpec,
                      "Unspec" = glmCentNonSpec)), 
           reducedDims(aeSceB)[[1]], k = 1)
phenoClustNamed <- sapply(phenoClust, switch, "1" = "Binder", "2" = "Non-binder")
table(phenoClustNamed, aeSceB$Specific)

#phenoClust  FALSE Not_tested TRUE
# Binder         1         35    4
# Non-binder    11        176    1

#So this shows, almost as one would have hoped, that the specific phenotype is highly
#enriched among the specific cells, but not exclusive to it, and conversely, one 
#ouf of five specific cells do not classify as specific in this way, so an 80% sensitivity
#and an 80% specificity as well. Given this, we can assume that at least 28 of the cells
#in the not tested group that classify as  are specific. We would estimate that 60 cells 
#should be specific, but these 35 should under all circumstances be highly enriched for them
#How do they distribute donor-wise?
table(phenoClustNamed, aeSceB$donor)

#phenoClustNamed 1166 1227 1284
#     Binder        3    4   33
#     Non-binder   30   43  115

#So intriguingly, desppite the truly binding cells are not more common in 1284 than in 1166, 
#there is a vast difference in the likelihood of being classified here if being a 1284 sample. 

