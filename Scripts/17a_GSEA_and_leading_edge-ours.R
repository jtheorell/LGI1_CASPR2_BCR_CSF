res <- readRDS("Data/ASC_vs_B/EdgeR_results.rds")

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

dir.create("Results/ASC_vs_B")
saveRDS(gse, "Data/ASC_vs_B/Gene_set_enrichment_GO_terms.rds")
#gse <- readRDS("Data/ASC_vs_B/Gene_set_enrichment_GO_terms.rds")

dotplot(gse, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/ASC_vs_B/Gene_set_enrichment_dotplot_ASC_vs_B.pdf", width = 8, height = 12)

gseTermSim <- pairwise_termsim(gse)
emapplot(gseTermSim, showCategory = 10)
ggsave("Results/ASC_vs_B/Emap_downsamp.pdf", width = 8, height = 7)

cnetplot(gse, categorySize="pvalue", foldChange=gene_list)

#ridgeplot(gse) + labs(x = "enrichment distribution")


#Now, we are investigating whether genes with an appealing and known druggable profile
#are overrepresented among these terms. 
targetsAE <- c("CD27", "CD38", "CXCR4", "CXCR5", "MS4A1", "CD19", "SDC1", "CD22", "TNFSF13B",
               "TNFRSF13C", "ANP32B", "TNFRSF13B", "IL12A", "IL12B", "IL12RB1", "IL12RB2", "IFNG",
               "IFNGR1","IFNGR2", "TNF", "TNFRSF1A", "TNFRSF1B", "IL10", "IL10RA","IL10RB", 
               "CSF2", "CSF2RA", "CSF2RB", "IL6", "IL6R", "BTK", "ITGA4", "ITGB1", "CD40")

#We start by excluding all genes that are downregulated in our data compared to the
#others, as such genes are by definition less druggable

posResult <- gse@result[which(gse@result$enrichmentScore >0),]

#This removes 151 downregulated GO terms, leaving us with 109. 

#Let us do this again for the full set. 
allTopGenesEnsembl <- unlist(unname(sapply(posResult$core_enrichment, function(x) strsplit(x, "\\/"))))

#These now need to be exchanged for their HGNC symbols
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")
allTopGenes <- rowData(aeSce)$hgnc_symbol[which(row.names(aeSce) %in% allTopGenesEnsembl)]

targHits <- length(which(unique(allTopGenes) %in% targetsAE))
#Here we get 14, distributing like this: 
table(allTopGenes[which(allTopGenes %in% targetsAE)])
#CD19      CD22     CXCR4     CXCR5    IFNGR1     IL12B       IL6     ITGB1     MS4A1       TNF TNFRSF13C 
#1         1         1         1         1         1         1         1         1         1         1 
#TNFRSF1A  TNFRSF1B  TNFSF13B 
#1         1         1 

nonTarghits <- length(unique(allTopGenes))-targHits
#1362

targGenes <- length(which(rowData(aeSce)$hgnc_symbol[which(row.names(aeSce) %in% names(gene_list))] %in% targetsAE))
#33
nonTargGenes <- length(gene_list)-targGenes
#15463

fisherDf <- data.frame("Non_target" = c(nonTargGenes-nonTarghits, nonTarghits),
                       "Target" = c(targGenes-targHits, targHits))
row.names(fisherDf) <- c("Non_hit", "Hit")
#So a fisher on this
fisher.test(fisherDf, alternative = "greater")
#p-value = 2.884e-07, so extremely significant, of course. 

#What about the top 20?

topResults <- posResult[1:20,]

topGenesEnsembl <- unlist(unname(sapply(topResults$core_enrichment, function(x) strsplit(x, "\\/"))))


topGenes <- rowData(aeSce)$hgnc_symbol[which(row.names(aeSce) %in% topGenesEnsembl)]

targHits <- length(which(topGenes %in% targetsAE))
#12

nonTarghits <- length(unique(topGenes))-targHits
#582

fisherDf <- data.frame("Non_target" = c(nonTargGenes-nonTarghits, nonTarghits),
                       "Target" = c(targGenes-targHits, targHits))
row.names(fisherDf) <- c("Non_hit", "Hit")
#So a fisher on this
fisher.test(fisherDf, alternative = "greater")
#p-value = 2.654e-08, so even significanter. 


###############
#GSEA mining
#Here, we are dependent on the tutorial: 
#http://www.bioconductor.org/packages/devel/bioc/vignettes/GSEAmining/inst/doc/GSEAmining.html#installation
#We start by reducing the number of gene sets, based on the normalised enrichment score values. 


gs.filt <- gm_filter(gse@result, 
                     p.adj = 0.05, 
                     neg_NES = 2, 
                     pos_NES = 2)

gs.cl <- gm_clust(gs.filt)

gm_dendplot(gs.filt, 
            gs.cl)

gs_enrich <- gm_enrichcores(gs.filt, gs.cl)

#Now, are our genes enriched among these clustered genes?
gs_enrich_res <- gs_enrich$data

targetsAEEnsembl <- row.names(aeSce)[which(rowData(aeSce)$hgnc_symbol %in% targetsAE)]

length(which(gs_enrich_res$lead_token %in% targetsAEEnsembl))
#1, namely TNF, so this will never be significant.

fisherDf <- data.frame("Non_target" = c(nonTargGenes-161, 161),
                       "Target" = c(targGenes-1, 1))
row.names(fisherDf) <- c("Non_hit", "Hit")
#So a fisher on this
fisher.test(fisherDf, alternative = "greater")
#p-value = 0.2933, not very surprisingly. 




