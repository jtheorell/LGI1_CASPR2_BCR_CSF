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

gene_list <- resList$Median$logFC
names(gene_list) <- row.names(resList$Median)

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
#                                      There were 46 pathways for which P-values were not calculated properly due to unbalanced (positive and negative) gene-level statistic values. For such pathways pval, padj, NES, log2err are set to NA. You can try to increase the value of the argument nPermSimple (for example set it nPermSimple = 10000)
#                                  3: In fgseaMultilevel(...) :
#                                      For some of the pathways the P-values were likely overestimated. For such pathways log2err is set to NA.                                  
saveRDS(gse, "Results/Gene_set_enrichment_analysis/GS_common.rds")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/Gene_set_enrichment_analysis/Gene_set_enrichment_dotplot_downsamp.pdf", width = 8, height = 12)

gseTermSim <- pairwise_termsim(gse)
emapplot(gseTermSim, showCategory = 10)
ggsave("Results/Gene_set_enrichment_analysis/Emap_downsamp.pdf", width = 8, height = 7)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list)

#ridgeplot(gse) + labs(x = "enrichment distribution")


#Now, we are investigating whether genes with an appealing and known drugable profile
#are overrepresented among these terms. 
targetsAE <- c("CD27", "CD38", "CXCR4", "CXCR5", "MS4A1", "CD19", "SDC1", "CD22", "TNFSF13B",
               "TNFRSF13C", "ANP32B", "TNFRSF13B", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
               "IFNGR1","IFNGR2", "TNF", "TNFRSF1A", "TNFRSF1B", "IL10", "IL10RA","IL10RB", 
               "CSF2", "CSF2RA", "CSF2RB", "IL6", "IL6R", "BTK", "ITGA4", "ITGB1", "CD40")

#We start by excluding all genes that are downregulated in our data compared to the
#others, as such genes are by definition less druggable
posResult <- gse@result[which(gse@result$enrichmentScore >0),]

#This removes 87 downregulated GO terms

#Let us do this again for the full set. 
allTopGenes <- unlist(unname(sapply(posResult$core_enrichment, function(x) strsplit(x, "\\/"))))

targHits <- length(which(unique(allTopGenes) %in% targetsAE))
#Here we get 2, but we have 4 occurrences of them. They dfistribute like this: 
table(allTopGenes[which(allTopGenes %in% targetsAE)])
#CD19  SDC1 
#1     3

nonTarghits <- length(unique(allTopGenes))-targHits
#1063
targGenes <- length(which(names(gene_list) %in% targetsAE))
#22
nonTargGenes <- length(gene_list)-targGenes
#7029

fisherDf <- data.frame("Non_target" = c(nonTargGenes-nonTarghits, nonTarghits),
                       "Target" = c(targGenes-targHits, targHits))
row.names(fisherDf) <- c("Non_hit", "Hit")
#So a fisher on this
fisher.test(fisherDf, alternative = "greater")
#p-value = 0.7858

#What about the top 20?

topResults <- posResult[1:20,]

topGenes <- unlist(unname(sapply(topResults$core_enrichment, function(x) strsplit(x, "\\/"))))

targHits <- length(which(topGenes %in% targetsAE))
#1
#CD19

nonTarghits <- length(unique(topGenes))-targHits
#582

fisherDf <- data.frame("Non_target" = c(nonTargGenes-nonTarghits, nonTarghits),
                       "Target" = c(targGenes-targHits, targHits))
row.names(fisherDf) <- c("Non_hit", "Hit")
#So a fisher on this
fisher.test(fisherDf, alternative = "greater")
#p-value = 0.9143


