#Now we will run a clusterProfiler gene set enrichment analysis
#See https://learn.gencore.bio.nyu.edu/rna-seq-analysis/gene-set-enrichment-analysis/
#This script is also identical to the one for the AE dataset, but just looped. 

#We start by getting the average Cohen-corrected fold change values for all markers. 
library(clusterProfiler)
library(enrichplot)
library(ggplot2)
library(ggnewscale)
library(DOSE)
library(Biobase)
library(GSEAmining)
library(org.Hs.eg.db)

fileList <- list.files("Data/ASC_vs_B", pattern = "EdgeR_results.rds", 
                       recursive = TRUE, full.names = TRUE)

#Here, we remove the AE file. 
fileList <- fileList[-grep("/AE/", fileList)]

resList <- lapply(fileList, readRDS)

names(resList) <- gsub("Data/ASC_vs_B/|/EdgeR_results.rds", "", fileList)

lapply(names(resList), function(x){
    locDat <- resList[[which(names(resList) == x)]]
    gene_list <- locDat$table$logFC
    names(gene_list) <- row.names(locDat)
    gene_list <- sort(gene_list, decreasing = TRUE)
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
    saveRDS(gse, paste0("Data/ASC_vs_B/", x, "/Gene_set_enrichment_GO_terms.rds"))
})

#Now, these are re-imported, together with ours, and we identify which of our
#significant GO terms that are unique to our cells. 

goResFileList <- list.files("Data/ASC_vs_B", recursive = TRUE, full.names = TRUE,
                            pattern = "Gene_set_enrichment_GO_terms.rds")

goResList <- lapply(goResFileList, readRDS)

names(goResList) <- gsub("Data/ASC_vs_B/|/Gene_set_enrichment_GO_terms.rds", 
                         "", goResFileList)

#goResSig <- lapply(goResList, function(x){
#    locRes <- x@result
#    #Here, we are only including the terms that are on the top, i.e.
#    #more than an absolute enrichment score of 2
#    locRes[which(abs(locRes$NES) > 2),]
#})

#What are the overlaps here?
goResTerms <- lapply(goResList, function(x) x$ID)
#Now many overlap in all? 
allGO <- Reduce(intersect, goResTerms)
#"GO:0008137" "GO:0034663" "GO:0072599" "GO:0045047" "GO:0003954" "GO:0050136" "GO:0018196"

aeMsGO <- Reduce(intersect, goResTerms[grep("MS|AE", names(goResTerms))])

aeCsfGO <- Reduce(intersect, goResTerms[grep("CSF|AE", names(goResTerms))])

aeNonMsGO <- Reduce(intersect, goResTerms[-grep("MS", names(goResTerms))])

aeShared <- unique(unlist(lapply(goResTerms[-1], function(x){
    x[which(x %in% goResTerms$AE)]
})))

aeUnique <- goResTerms$AE[-which(goResTerms$AE %in% aeShared)]

#Now, we try to plot the top ones in this subgroup
aeRes <- goResList$AE

aeUniqueRes <- aeRes

aeUniqueRes@result <- aeRes@result[which(goResList$AE@result$ID %in% aeUnique),]

dotplot(aeUniqueRes, showCategory=10, split=".sign") + facet_grid(.~.sign)
ggsave("Results/ASC_vs_B/AE/Gene_set_enrichment_dotplot_unique_GOs.pdf", width = 8, height = 12)


#What are the contributions/overlaps to the different groups? We of course need to
#take the different asay qualities into account here, so we will normalise to the
#total number of significant GO terms. 

aeShared <- unlist(lapply(goResTerms[-1], function(x){
    length(x[which(x %in% goResTerms$AE)])/length(x)
}))

aeShared
#Covid_PBMC      Healthy_PBMC    Influensa_PBMC     MS_Ramesh_CSF    MS_Ramesh_PBMC  MS_Shafflick_CSF 
#0.5917266         0.3378378         0.6105263         0.5977011         0.6917293         0.4909366 
#MS_Shafflick_PBMC 
#0.6947368


#So, as expected after the gene investigateions, we see the lowest similarity to the healthy dataset. 
#interestingly, there is not celar trent towards a clearer overlap to the MS CSF cells. Infact, 
#the Shafflick CSF dataset has the lowest overlap percentage-wise. 

ggDat <- data.frame("Dataset" = names(aeShared),
                    "Fraction" = unname(aeShared))

ggplot(data=ggDat, aes(x=Dataset, y=Fraction)) +
    geom_bar(stat="identity", position=position_dodge()) +
    theme_bw() + scale_y_continuous(expand = c(0, 0)) + ylab("Fraction of hits among AE top GO")
ggsave("Results/ASC_vs_B/Shared_top_GO_fraction.pdf")
