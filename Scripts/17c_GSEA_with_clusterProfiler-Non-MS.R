resList <- readRDS("Results/Comp_to_others/All_DE_integrated_downsamp.rds")

#Now we will run a clusterProfiler gene set enrichment analysis
#See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#We start by getting the average Cohen-corrected fold change values for all markers. 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(DOSE)
library(Biobase)

resList <- resList[-which(names(resList) %in% c("Median",
                                                "MS_Ramesh","MS_Schafflick"))]
    
targetList <- lapply(1:length(resList), function(x){
    locDf <- data.frame(row.names(resList[[x]]),
                        resList[[x]]$logFC)
    colnames(locDf) <- c("Name", paste0(names(resList)[[x]], "_logFC"))
    locDf <- locDf[order(locDf$Name),]
})

all(sapply(targetList, function(x){
    identical(targetList[[1]]$Name, x$Name)
}))
#TRUE

targetDf <- data.frame(do.call("cbind", lapply(targetList, "[[", 2)),
                       row.names = targetList[[1]]$Name)

colnames(targetDf) <- names(resList)

targetDf$medianFC <- rowMedians(as.matrix(targetDf))

gene_list <- targetDf$medianFC
names(gene_list) <- row.names(targetDf)

gene_list <- sort(gene_list, decreasing = TRUE)
#We also remove all the ones with a value of 0
gene_list <- gene_list[-which(gene_list == 0)]

library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
set.seed(101)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             eps = 0)

#Warning messages:
#    1: In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                      There are ties in the preranked stats (0.01% of the list).
#                                  The order of those tied genes will be arbitrary, which may produce unexpected results.
#                                  2: In fgseaMultilevel(...) :
#                                      There were 45 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
#                                  3: In fgseaMultilevel(...) :
#                                      For some of the pathways the P-values were likely overestimated. For such pathways log2err is set to NA.
                                  
saveRDS(gse, "Results/Gene_set_enrichment_analysis/GS_non-MS.rds")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/Gene_set_enrichment_analysis/Gene_set_enrichment_dotplot_downsamp_non-MS.pdf", width = 8, height = 12)

gseTermSim <- pairwise_termsim(gse)
emapplot(gseTermSim, showCategory = 10)
ggsave("Results/Gene_set_enrichment_analysis/Emap_downsamp_non-MS.pdf", width = 8, height = 7)

