library(scRNAseq)
library(SingleCellExperiment)
library(scran)
library(scater)
#This dataset is from the study https://www.nature.com/articles/s41591-020-0769-8
#Broad immune activation underlies shared set point signatures for vaccine 
#responsiveness in healthy individuals and disease activity in patients with lupus
#By cotliarov et al. 
influensaSce <- KotliarovPBMCData(mode = c("rna", "adt"), ensembl = TRUE, location = TRUE)

#Here, we are going to use the protein data to define the cells according to our 
#definitions. But first the normal stuff. 
influensaSce  <- getBMFeatureAnnos(influensaSce, filters = "ensembl_gene_id", 
                               attributes = c("ensembl_gene_id", "hgnc_symbol",
                                              "chromosome_name", "gene_biotype", 
                                              "start_position", "end_position",
                                              "transcript_length"), 
                               dataset = "hsapiens_gene_ensembl")

saveRDS(influensaSce, "../External/Data/Influensa/Influensa/0_raw.rds")

mito <- which(rowData(influensaSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(influensaSce, subsets=list(Mt=mito))

colData(influensaSce) <- cbind(colData(influensaSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

colSums(as.matrix(reasons))
#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#4392                    4448                    4896                    5891

#So 10% are lost. Reasonable. 
influensaSce$discard <- reasons@listData$discard

dir.create("Diagnostics/Comp_to_others/Influensa", recursive = TRUE)

influensaSce$donor <- gsub(".+_.....|", "", colnames(influensaSce))

gridExtra::grid.arrange(
    plotColData(influensaSce, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(influensaSce, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(influensaSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/Influensa/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

influensaSceRed <- influensaSce[,-which(influensaSce$discard)]

#NORMALIZATION
influensaSce <- computeSumFactors(influensaSce) 
influensaSce <- logNormCounts(influensaSce)

#Now, we exclude the rest of the cells
saveRDS(influensaSce, "../External/Data/Influensa/1_logNorm.rds")

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- influensaSce[which(rowData(influensaSce)$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
dir.create("Data/Comp_to_others/Influensa", recursive = TRUE)

saveRDS(pred, "Data/Comp_to_others/Influensa/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/Influensa/MonacoSingler.rds")

unique(pred$pruned.labels)
#[1] "Natural killer cells"          "Naive CD8 T cells"             "Naive CD4 T cells"            
#[4] "Classical monocytes"           "Terminal effector CD8 T cells" "Th17 cells"                   
#[7] "Naive B cells"                 "Th1/Th17 cells"                "Th1 cells"                    
#[10] "Vd2 gd T cells"                "MAIT cells"                    "Intermediate monocytes"       
#[13] "Effector memory CD8 T cells"   "Central memory CD8 T cells"    "Terminal effector CD4 T cells"
#[16] NA                              "T regulatory cells"            "Th2 cells"                    
#[19] "Exhausted B cells"             "Non-switched memory B cells"   "Non classical monocytes"      
#[22] "Myeloid dendritic cells"       "Follicular helper T cells"     "Non-Vd2 gd T cells"           
#[25] "Plasmacytoid dendritic cells"  "Low-density basophils"         "Switched memory B cells"      
#[28] "Progenitor cells"              "Plasmablasts"                  "Low-density neutrophils"   

influensaSce$cellType <- pred$pruned.labels
influensaSceB <- influensaSce[,which(pred$pruned.labels %in% c("Exhausted B cells",
                                                                "Naive B cells",
                                                                "Non-switched memory B cells",
                                                                "Switched memory B cells",
                                                                "Plasmablasts"))]
influensaSceB$cellType[grep("B", influensaSceB$cellType)] <- "B"
influensaSceB$cellType[grep("Plasmablasts", influensaSceB$cellType)] <- "ASC"
table(influensaSceB$cellType)
#ASC    B 
#66 6132 

saveRDS(influensaSceB, "../External/Data/Influensa/2_normalised_post_monaco.rds")



