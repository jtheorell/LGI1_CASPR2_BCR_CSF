res <- readRDS("Results/Specific_UCA_vs_not/EdgeR_result.rds")

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

#Here, 0 enriched terms were found. Clean. 
saveRDS(gse, "Results/Specific_UCA_vs_not/Gene_set_enrichment_GO_terms.rds")
