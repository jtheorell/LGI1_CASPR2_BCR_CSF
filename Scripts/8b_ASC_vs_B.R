#As the ASC seem to be more or less identical, and as there is such a clear divide
#of the cell types, where the non-specific ASC make up only 13 percent of the total
#pool, and the specific B-cells, even in the IgG compartment, make up only 22% of
#the cells, we will now make the slightly dum comparison of ASC vs B, to get to 
#questions about druggability, etc. 

library(SingleCellExperiment)
library(scran)
library(scuttle)
library(scater)
library(BiocParallel)
library(edgeR)

aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")
#Now, the exclusions: 
aeSce <- aeSce[,which(aeSce$clustersLouvain != "CD8T")]
aeSce <- aeSce[,-which(aeSce$clustersLouvain == "ASC" & aeSce$Specific != "TRUE")]

table(aeSce$donor, aeSce$clustersLouvain)
#     ASC   B
#1166  84  36
#1227  19  48
#1284  40 149

#This is the most sensible comparison yet, from a size perspective. 
current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = counts(aeSce))), 
                                                     ids=colData(aeSce)[,c("clustersLouvain", "donor")])
current
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))
y
#summary(discarded)
#THis is for iteration 1, dataset 1
#   Mode   FALSE    TRUE
# logical      16      16
y <- calcNormFactors(y)
design <- model.matrix(~factor(current$clustersLouvain), y$donor)

#design
y <- estimateDisp(y, design)
#summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))


#Which are the most different genes? 
hgncSymbs <- rowData(aeSce)$hgnc_symbol[which(row.names(aeSce) %in% row.names(topTags(res)))]
locRowNames <- row.names(aeSce)[which(row.names(aeSce) %in% row.names(topTags(res)))]
hgncSymbsOrd <- hgncSymbs[match(row.names(topTags(res)), locRowNames)]
topRes <- topTags(res)[[1]]
row.names(topRes) <- hgncSymbsOrd

topRes
#This is seen from the B-cell side. 
#logFC    logCPM        F       PValue          FDR
#PYCR1   -15.859232  7.482947 313.6332 2.638084e-08 0.0005235541
#SDC1    -15.520284  8.271722 228.3726 1.212212e-07 0.0012028782
#LY86     13.988939  6.424052 198.1698 1.903606e-07 0.0012592989
#KCNN3    -9.714018  7.090833 166.1168 4.469034e-07 0.0019563803
#CNR2      9.052579  6.964555 156.7432 5.348300e-07 0.0019563803
#PPIB     -4.840455 11.879392 151.9892 5.914684e-07 0.0019563803
#BCL11A    8.577639  8.106796 145.8075 8.178125e-07 0.0023186154
#HSP90B1  -4.539244 10.716997 128.0773 1.223022e-06 0.0030340110
#CAV1    -15.054732  6.668163 121.4085 1.761624e-06 0.0038076347
#NIBAN3    7.874243  6.334814 113.8277 2.014696e-06 0.0038076347

#So this is the result we are after. 
dir.create("Data/ASC_vs_B")
saveRDS(res, "Data/ASC_vs_B/EdgeR_results.rds")

