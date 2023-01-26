#The method to use here is slinshot. See: 
#https://blog.bioturing.com/2022/06/13/single-cell-rna-seq-trajectory-analysis-review/

#See also: 
#https://www.bioconductor.org/packages/release/bioc/vignettes/slingshot/inst/doc/vignette.html

library(slingshot)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_4_BCR_plus_all_others.rds")

csfSceB <- csfSce[,-which(csfSce$clustersLouvain == "CD8T")]
bFlow <- normcounts(altExp(csfSceB, "flowData"))
plot(bFlow["IgD",], bFlow["CD27",])

# filter genes down to potential cell-type markers
# at least M (3) reads in at least N (10) cells
geneFilter <- apply(assays(csfSceB)$counts,1,function(x){
    sum(x >= 3) >= 10
})
csfSceB <- csfSceB[geneFilter, ]

csfSceB <- slingshot(csfSceB, clusterLabels = 'clustersLouvain', reducedDim = 'flowAndGLMPCA')

dir.create("Results/Slingshot")
plotReducedDim(csfSceB, "combinedUmap", colour_by = "slingPseudotime_1", shape_by = "clustersLouvain")
ggsave("Results/Slingshot/Slingshot_pseudotime_on_combined_UMAP.pdf")
DepecheR::dColorPlot(bFlow["IgD",], xYData = reducedDim(csfSceB, "combinedUmap"),
                     plotDir = "Results/Slingshot/", plotName = "IgD_on_combined_Umap")



#So this now makes perfect sense with the ASC being older than the B. 

#Now, what about specificity vs slingshot?
nonNASce <- csfSceB[,-which(csfSceB$Specific == "Not_tested" | is.na(csfSceB$Specific))]
densDat <- data.frame("Specificity" = nonNASce$Specific[which(nonNASce$clustersLouvain == "ASC")],
                      "Pseudotime" = nonNASce$slingPseudotime_1[which(nonNASce$clustersLouvain == "ASC")])
ggplot(densDat, aes(x=Pseudotime, fill=Specificity)) +
    geom_density(alpha=0.4) + theme_bw()
ggsave("Results/Slingshot/Pseudotime_vs_specificity_ASC.pdf")






#Over to the question of plotting the trajectories on the trees. 

#We of course also need to overlay this on the hierarchical trees.
trees <- readRDS("Data/BCR_auxiliaries/IgphyML_trees.rds")

#Now, we bin the pseudotime. 
floatPseudo <- csfSceB$slingPseudotime_1
definedPseudo <- .bincode(floatPseudo, seq(0,max(floatPseudo), by = max(floatPseudo)/10), include.lowest = T)
csfSceB$binnedPseudo <- as.factor(definedPseudo)
library(viridis)

trees$data <- lapply(trees$data, function(x){
    x@data$pseudotime <- unlist(sapply(x@data$sequence_id, function(y){
        locName <- gsub("|H", "", y)
        if(locName %in% colnames(csfSceB)){
            csfSceB$binnedPseudo[which(colnames(csfSceB) == locName)]
        } else {
            "Not_present"
        }
    }))
    x
})

plots = plotTrees(trees, tips="pseudotime", tipsize = 3, common_scale = TRUE,
                  labelsize = 8, tip_palette = c("10" ="#FCFFA4FF",
                                                 "9" = "#F6D645FF",
                                                 "8" = "#FCA50AFF",
                                                 "7" = "#F3771AFF",
                                                 "6" = "#DD513AFF",
                                                 "5" ="#BB3754FF",
                                                 "4" = "#932667FF",
                                                 "3" = "#6B186EFF",
                                                 "2" = "#420A68FF",
                                                 "1" = "#170C3AFF",
                                                 "Germline" = "#000004FF",
                                                 "Not_present" = "grey"))

treesToPDF(plots, file="Results/Trees/BCR_trees_with_slingshot.pdf",
           nrow=2, ncol=2)

#Now, do trees with high intraclonal celocuty diversity differ in specificity from the others?
cloneSce <- csfSceB[,which(csfSceB$Clonal)]

perCloneSlingDivergence <- unlist(lapply(split(cloneSce$slingPseudotime_1, cloneSce$HamClone), function(x){
    if(any(is.na(x))){
        locRange <- range(x[-which(is.na(x))])
    } else {
        locRange <- range(x)
    }
    locRange[2]-locRange[1]
}))

densDat <- data.frame("Specificity" = unlist(lapply(split(cloneSce$Specific, cloneSce$HamClone), "[[",1)),
                      "Pseudotime_divergence" = perCloneSlingDivergence)
ggplot(densDat, aes(y=Pseudotime_divergence, x = Specificity, fill=Specificity)) +
    geom_dotplot(binaxis='y', stackdir='center') + theme_bw()
ggsave("Results/Slingshot/Intraclonal_pseudotime_divergence.pdf", width = 8, height = 5)                            


