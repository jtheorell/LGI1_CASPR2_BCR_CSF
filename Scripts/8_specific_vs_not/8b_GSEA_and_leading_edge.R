res <- readRDS("Results/Specific_vs_not/EdgeR_result.rds")

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

saveRDS(gse, "Results/Specific_vs_not/Gene_set_enrichment_GO_terms.rds")
#gse <- readRDS("Data/ASC_vs_B/Gene_set_enrichment_GO_terms.rds")


dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/Specific_vs_not/Gene_set_enrichment_dotplot_ASC_vs_B.pdf", width = 8, height = 12)

gseTermSim <- pairwise_termsim(gse)
emapplot(gseTermSim, showCategory = 10)
ggsave("Results/Specific_vs_not/Emap_downsamp.pdf", width = 8, height = 7)
