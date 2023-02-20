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
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_3_ASC_only.rds")

aeSceSpecKnown <- aeSce[,-which(aeSce$Specific == "Not_tested")]

#Now, we will remove all without a HGNC symbol and then change the names to them
dupHgnc <- rowData(aeSce)$hgnc_symbol[which(duplicated(rowData(aeSce)$hgnc_symbol))]

table(dupHgnc)
#      ABCF2    AHRR   ATXN7    GGT1   PINX1 POLR2J3    TBCE 
#606       1       1       1       1       1       1       1
aeSceHgnc <- aeSceSpecKnown[-which(rowData(aeSceSpecKnown)$hgnc_symbol %in% dupHgnc),]
row.names(aeSceHgnc) <- rowData(aeSceHgnc)$hgnc_symbol

##################
#EdgeR analysis
##################

current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = counts(aeSceHgnc))), 
                                ids=colData(aeSceHgnc)[,c("Specific", "donor")])
current
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))
y

#Now, we are bringing this in for official reasons, to make the script comparible
#also for the non-AE data, but we will not need to throw out any sample in our case
discarded <- current$ncells < 4
y <- y[,!discarded]
summary(discarded)
#   Mode   FALSE 
#logical       5

keep <- filterByExpr(y, group=current$Specific)
y <- y[keep,]
summary(keep)

#   Mode   FALSE    TRUE 
#logical    9246    9974

y <- calcNormFactors(y)
design <- model.matrix(~factor(Specific), y$samples)

#design
y <- estimateDisp(y, design)
#summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))

dir.create("Results/Specific_vs_not")

saveRDS(res, "Results/Specific_vs_not/EdgeR_result.rds")

resSig <- res$table[which(abs(res$table$logFC) > 2 & res$table$PValue < 0.05),]

#After investigating this a bit further, we do on top of the above criteria
#also include an additional criterion: as we have so few cells, we require the hits
#to be significantly different also from a fisher count persective. 
aeSig <- aeSceHgnc[which(row.names(aeSceHgnc) %in% row.names( resSig)),]

specAll <-length(which(aeSig$Specific == "TRUE"))
nonSpecAll <-length(which(aeSig$Specific == "FALSE"))

fishSig <- sapply(row.names(aeSig), function(x){
    specExpr <- length(which(aeSig$Specific == "TRUE" & counts(aeSig)[which(row.names(aeSig) == x),] > 0))
    specNonExpr <- specAll-specExpr
    nonSpecExpr <- length(which(aeSig$Specific == "FALSE" & counts(aeSig)[which(row.names(aeSig) == x),] > 0))
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

resOrd <- resSigReal[order(abs(resSigReal$logFC), decreasing = TRUE),]

write.csv(resOrd, "Results/Specific_vs_not/Top_results.csv")

