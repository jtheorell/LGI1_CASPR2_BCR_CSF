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
library(BiocParallel)
#In these analyses, we will stop calling our data "csfSce" as there are two such
#datasets now. Instead aeSce for "Autoimmune encephalitis"
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

aeSceSpecKnown <- aeSce[,-which(aeSce$Specific == "Not_tested")]

marker.info <- scoreMarkers(aeSceSpecKnown, aeSceSpecKnown$Specific)
marker.info
chosen <- marker.info[["TRUE"]]

#And now, we slightly controversially calculate the Mann-Whitney p-values. 
#The reason this is controversial is that by doing this, we assume that each 
#observation, i.e. cell, is idepentent, which is not true, as all cells that
#come from one individual are by definition dependent. But we will do so anyway
#as this is the norm in the field. 
#We start by removing all that have fewer than four cells with expression. 
#This is necessary to remove all NaNs from the data. 
numPosCells <- apply(counts(aeSceSpecKnown), 1, function(x) sum(x != 0))

aeSceNonZero <- aeSceSpecKnown[which(numPosCells > 3),]
spec <- logcounts(aeSceNonZero[,which(aeSceNonZero$Specific == "TRUE")])
nonSpec <- logcounts(aeSceNonZero[,which(aeSceNonZero$Specific == "FALSE")])

locP <- lapply(seq_along(row.names(aeSceNonZero)), function(y) {
    wilcox.test(spec[y,], nonSpec[y,],
                exact = FALSE)$p.value
})

#Here, we check that we can integrate the data
identical(row.names(chosen), row.names(aeSceSpecKnown))
#TRUE

pDf <- data.frame("HGNC" = rowData(aeSceNonZero)$hgnc_symbol, "P_val" = unlist(locP),
                  "logFC" = chosen$median.logFC.detected[which(row.names(chosen) %in% row.names(aeSceNonZero))],
                  row.names = row.names(aeSceNonZero))

pDf$FDR <- p.adjust(pDf$P_val, method = "fdr")

#THis is now plotted as a volcano plot
ggplot(data=pDf, aes(x=logFC, y=-log10(FDR))) + 
    geom_point() + 
    theme_bw()
ggsave("Results/Volcano_spec_vs_non-spec_no_findings.pdf", width = 6, height = 5)

#And we also plot the top hits here, just for the sake of it, really showing the non-differences
topGenes <- which(-log10(pDf$FDR) > 0.6)
#topGenes <- which(pDf$logFC < -1.6)
smallDataSet <- aeSceNonZero[topGenes,]
row.names(smallDataSet) <- rowData(smallDataSet)$hgnc_symbol
plotExpression(smallDataSet, features=row.names(smallDataSet), 
               x="Specific", colour_by="Specific", ncol = 3)
dev.copy(pdf,'Results/Hit_genes_that_are_not_hit_genes.pdf', height = 4, width = 6)
dev.off()





