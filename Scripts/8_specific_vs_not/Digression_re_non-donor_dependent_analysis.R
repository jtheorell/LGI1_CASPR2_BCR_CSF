


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
ggsave("Results/Specific_vs_not/Volcano_spec_vs_non-spec_no_findings.pdf", width = 6, height = 5)

#And we also plot the top hits here, just for the sake of it, really showing the non-differences
topGenes <- which(-log10(pDf$FDR) > 0.6)
#topGenes <- which(pDf$logFC < -1.6)
smallDataSet <- aeSceNonZero[topGenes,]
row.names(smallDataSet) <- rowData(smallDataSet)$hgnc_symbol
plotExpression(smallDataSet, features=row.names(smallDataSet), 
               x="Specific", colour_by="Specific", ncol = 3)
dev.copy(pdf,'Results/Specific_vs_not/Hit_genes_that_are_not_hit_genes.pdf', height = 4, width = 6)
dev.off()

