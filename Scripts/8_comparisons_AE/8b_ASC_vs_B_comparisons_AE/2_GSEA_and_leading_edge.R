res <- readRDS("Data/ASC_vs_B/AE/EdgeR_results.rds")

#Now we will run a clusterProfiler gene set enrichment analysis
#See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/

#We start by getting the average Cohen-corrected fold change values for all markers. 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(DOSE)
library(Biobase)
library(GSEAmining)

gene_list <- res$table$logFC
names(gene_list) <- row.names(res)

gene_list <- sort(gene_list, decreasing = TRUE)
#We also remove all the ones with a value of <5
gene_list <- gene_list[-which(gene_list == 0)]

library(org.Hs.eg.db)
#keytypes(org.Hs.eg.db)
set.seed(101)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "ENSEMBL", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             eps = 0)

#Warning messages:
#Warning message:
#    In preparePathwaysAndStats(pathways, stats, minSize, maxSize, gseaParam,  :
#                                   There are ties in the preranked stats (5.14% of the list).
#                               The order of those tied genes will be arbitrary, which may produce unexpected results.                                  

dir.create("Results/ASC_vs_B/AE", recursive = TRUE)
saveRDS(gse, "Data/ASC_vs_B/AE/Gene_set_enrichment_GO_terms.rds")
#gse <- readRDS("Data/ASC_vs_B/Gene_set_enrichment_GO_terms.rds")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/ASC_vs_B/AE/Gene_set_enrichment_dotplot_ASC_vs_B.pdf", width = 8, height = 12)

gseTermSim <- pairwise_termsim(gse)
emapplot(gseTermSim, showCategory = 10)
ggsave("Results/ASC_vs_B/AE/Emap_downsamp.pdf", width = 8, height = 7)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list)

#ridgeplot(gse) + labs(x = "enrichment distribution")
