#Here comes now a combined analysis. Previous analyses have shown that, at least
#when using Batchelor, the batch correction completely removes all differences 
#between the ASC in the CSF and in the periphery, possibly as they are so comparatively
#homogenous that any nearest neighbor analysis will anchor the whole dataset to itself.
#Now the strategy here is instead to use two separate peripheral dataset: one from
#healthy individuals and one from Dengue infected ditto, and not do any batch
#correction but only consider the hits that are consistent between the two 
#datasets. This of course builds on the assumption that there are as large batch
#differences between the dengue and healthy datasets as it is between our dataset
#and the two others, but that will have to do.
library(SingleCellExperiment)
library(scran)
library(scater)
library(edgeR)
#In these analyses, we will stop calling our data "csfSce" as there are two such
#datasets now. Instead aeSce for "Autoimmune encephalitis"
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")


aeSceSpecKnown <- aeSce[,-which(aeSce$Specific %in% c("Not_tested", "FALSE"))]

table(aeSceSpecKnown$Specific_UCA, aeSceSpecKnown$donor)
#      1166 1227 1284
#FALSE   56    5   27
#TRUE    19    8    4
#So this analysis will work technically. 

#Now, we will remove all without a HGNC symbol and then change the names to them
dupHgnc <- rowData(aeSceSpecKnown)$hgnc_symbol[which(duplicated(rowData(aeSceSpecKnown)$hgnc_symbol))]

table(dupHgnc)
#      ABCF2    AHRR   ATXN7    GGT1   PINX1 POLR2J3    TBCE 
#606       1       1       1       1       1       1       1
aeSceHgnc <- aeSceSpecKnown[-which(rowData(aeSceSpecKnown)$hgnc_symbol %in% dupHgnc),]
row.names(aeSceHgnc) <- rowData(aeSceHgnc)$hgnc_symbol

#Here, we exclude all cells 
##################
#EdgeR analysis
##################

current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = counts(aeSceHgnc))), 
                                ids=colData(aeSceHgnc)[,c("Specific_UCA", "donor")])
current
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))
y

#Now, we are bringing this in for official reasons, to make the script comparible
#also for the non-AE data, but we will not need to throw out any sample in our case
discarded <- current$ncells < 4
y <- y[,!discarded]
summary(discarded)

keep <- filterByExpr(y, group=current$Specific_UCA)
y <- y[keep,]
summary(keep)

#   Mode   FALSE    TRUE 
#logical.  10165    9055

y <- calcNormFactors(y)
design <- model.matrix(~factor(Specific_UCA), y$samples)

#design
y <- estimateDisp(y, design)
#summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))

dir.create("Results/Specific_UCA_vs_not")

saveRDS(res, "Results/Specific_UCA_vs_not/EdgeR_result.rds")

topTags(res)

resSig <- res$table[which(abs(res$table$logFC) > 2 & res$table$PValue < 0.05),]

#Now, as we are dealing with very, very low cell numbers, we are also going to make
#a further filtering step, namely exclusion of all genes that are not significantly
#regulated in a simple Fisher test. 
aeUcaSig <- aeSceHgnc[which(row.names(aeSceHgnc) %in% row.names( resSig)),]

specAll <-length(which(aeUcaSig$Specific_UCA == "TRUE"))
nonSpecAll <-length(which(aeUcaSig$Specific_UCA == "FALSE"))
fishSig <- sapply(row.names(aeUcaSig), function(x){
    specExpr <- length(which(aeUcaSig$Specific_UCA == "TRUE" & counts(aeUcaSig)[which(row.names(aeUcaSig) == x),] > 0))
    specNonExpr <- specAll-specExpr
    nonSpecExpr <- length(which(aeUcaSig$Specific_UCA == "FALSE" & counts(aeUcaSig)[which(row.names(aeUcaSig) == x),] > 0))
    nonSpecNonExpr <- nonSpecAll-nonSpecExpr
    fisherDf <- data.frame("Non_spec" = c(nonSpecNonExpr, nonSpecExpr),
                           "Spec_ASC" = c(specNonExpr, specExpr))
    row.names(fisherDf) <- c("Not_expressed", "Expressed")
    if( resSig$logFC[which(row.names(resSig) == x)] > 0){
        fisher.test(fisherDf, alternative = "greater")$p.value
    } else {
        fisher.test(fisherDf, alternative = "less")$p.value
    }
})

resSigReal <- resSig[which(fishSig < 0.05),]

resOrd <- resSig[order(abs(resSigReal$logFC), decreasing = TRUE),]
volcDf <- res$table
volcDf$Direction <- "None"
volcDf$Direction[which(row.names(volcDf) %in% row.names(resOrd) & volcDf$logFC > 0)] <- "Up"
volcDf$Direction[which(row.names(volcDf) %in% row.names(resOrd) & volcDf$logFC < 0)] <- "Down"

volcDf$Direction <- factor(volcDf$Direction, levels = c("Up", "None", "Down"))
ggplot(data=volcDf, aes(x=logFC, y=-log10(PValue), 
                        col=Direction)) + 
    geom_point() + 
    theme_bw() +
    scale_color_manual(values = c("red", "grey", "blue")) + xlim(-20, 20) + ylim(0, 8)
ggsave("Results/Specific_UCA_vs_not/Volcano_top_hits.pdf", width = 6, height = 5)

write.csv(resOrd, "Results/Specific_UCA_vs_not/Top_results.csv")

