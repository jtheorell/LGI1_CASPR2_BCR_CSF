library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(data.table)

#This data comes from the study https://www.nature.com/articles/s41467-019-14118-w#Sec10
#Integrated single cell analysis of blood and cerebrospinal fluid leukocytes in multiple sclerosis

#Data from https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138266

matrixNameList <- list.files("../External/Data/MS_Schafflick_PBMC/GSE138266_PBMC", pattern = "matrix", 
                         full.names = TRUE)

matrixList <- lapply(matrixNameList, readMM)

#For the gene names, they should all be identical, so we will check if that is the case for two

identical(fread("../External/Data/MS_Schafflick_PBMC/GSE138266_PBMC/GSM4104134_MS19270_PBMCs_GRCh38_genes.tsv.gz"),
          fread("../External/Data/MS_Schafflick_PBMC/GSE138266_PBMC/GSM4104135_MS71658_PBMCs_GRCh38_genes.tsv.gz"))
#TRUE

#So then we combine it all into one big dataframe after giving the cells some 
#colnames and then all get these rownames. 
fileNameList <- list.files("../External/Data/MS_Schafflick_PBMC/GSE138266_PBMC", pattern = "matrix")

matrixListNamed <- lapply(seq_along(matrixList), function(x){
    locMat <- matrixList[[x]]
    colnames(locMat) <- paste0(substr(fileNameList, 1, 10)[x],
                               "_cell_", seq(1,ncol(locMat)))
    locMat
})

matrixComined <- do.call("cbind", matrixListNamed)

rowDataMS <-  fread("../External/Data/MS_Schafflick_PBMC/GSE138266_PBMC/GSM4104134_MS19270_PBMCs_GRCh38_genes.tsv.gz",
                    header = FALSE)

row.names(matrixComined) <- rowDataMS$V1

msSchafSce<- SingleCellExperiment(assays = list(counts = matrixComined))

rowData(msSchafSce)$original_hgnc <- rowDataMS$V2

msSchafSce  <- getBMFeatureAnnos(msSchafSce, filters = "ensembl_gene_id", 
                             attributes = c("ensembl_gene_id", "hgnc_symbol",
                                            "chromosome_name", "gene_biotype", 
                                            "start_position", "end_position"), 
                             dataset = "hsapiens_gene_ensembl")

saveRDS(msSchafSce, "../External/Data/MS_Schafflick_PBMC/0_raw.rds")

############
#EXCLUSIONS
###########
mito <- which(rowData(msSchafSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(msSchafSce, subsets=list(Mt=mito))

#At this stage, we need to remove all completely empty droplets (sum 0), that make up a very considerable
#portion of this dataset (about 82%)

colData(msSchafSce) <- stats

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

msSchafSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#.        62                     600                     471                    1058 

#Here, there are none to discard. 
msSchafSce$donor <- substr(colnames(msSchafSce), 1, 10)

#Here we make a subset to plot
msSchafSceSub <- msSchafSce[,sample(1:ncol(msSchafSce), 10000)]

dir.create("Diagnostics/Comp_to_others/MS_Schafflick_PBMC", recursive = TRUE)

gridExtra::grid.arrange(
    plotColData(msSchafSceSub, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(msSchafSceSub, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(msSchafSceSub, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/MS_Schafflick_PBMC/Exclusion_metrics_per_donor_automatic.pdf', height = 5, width = 15)
dev.off()

msSchafSce<- msSchafSce[,-which(msSchafSce$discard)]
saveRDS(msSchafSce, "../External/Data/MS_Schafflick_PBMC/1_lowQual_excluded.rds")

msSchafSce<- computeSumFactors(msSchafSce)
msSchafSce<- logNormCounts(msSchafSce)

assayNames(msSchafSce)

#Here, a group column is introduced
groupCode <- read.csv("../External/Data/MS_Schafflick_PBMC/Group_code.csv")

msSchafSce$Group <- FALSE
for(i in groupCode$Four_last){
    msSchafSce$Group[which(substr(colnames(msSchafSce), 7,10) == i)] <- 
        groupCode$Group[which(groupCode$Four_last == i)]
}

table(msSchafSce$Group)
#  MS 
#24773 
#Good, no controls still in the system. 

#Now, we will head straight for the singler analysis. 
library(SingleR)
library(celldex)
#This ref is taken from: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE107011
ref <- MonacoImmuneData()
singlerSce <- msSchafSce [which(rowData(msSchafSce )$hgnc_symbol %in% row.names(ref)),]

rownames(singlerSce) <- rowData(singlerSce)$hgnc_symbol
pred <- SingleR(test=singlerSce, ref=ref, labels=ref$label.fine)
table(pred$pruned.labels)
dir.create("Data/Comp_to_others/MS_Schafflick_PBMC/", recursive = TRUE)

saveRDS(pred, "Data/Comp_to_others/MS_Schafflick_PBMC/MonacoSingler.rds")
#pred <- readRDS("Data/Comp_to_others/MS_Schafflick/MonacoSingler.rds")
unique(pred$pruned.labels)
#[1] "Naive CD4 T cells"             "Classical monocytes"           "T regulatory cells"           
#[4] "Th1/Th17 cells"                "Th2 cells"                     "Non-switched memory B cells"  
#[7] NA                              "Intermediate monocytes"        "Th1 cells"                    
#[10] "Natural killer cells"          "Plasmacytoid dendritic cells"  "Naive B cells"                
#[13] "Switched memory B cells"       "Vd2 gd T cells"                "Naive CD8 T cells"            
#[16] "Central memory CD8 T cells"    "MAIT cells"                    "Terminal effector CD8 T cells"
#[19] "Th17 cells"                    "Non-Vd2 gd T cells"            "Effector memory CD8 T cells"  
#[22] "Non classical monocytes"       "Myeloid dendritic cells"       "Terminal effector CD4 T cells"
#[25] "Follicular helper T cells"     "Plasmablasts"                  "Progenitor cells"             
#[28] "Exhausted B cells"             "Low-density basophils"   

msSchafSce$cellType <- pred$pruned.labels
msSchafSceB <- msSchafSce[,which(msSchafSce$cellType %in% c("Exhausted B cells",
                                                               "Naive B cells",
                                                               "Non-switched memory B cells",
                                                               "Switched memory B cells",
                                                               "Plasmablasts"))]
msSchafSceB$cellType[grep("B", msSchafSceB$cellType)] <- "B"
msSchafSceB$cellType[grep("Plasmablasts", msSchafSceB$cellType)] <- "ASC"
table(msSchafSceB$cellType)
#ASC    B 
#36 2385 

saveRDS(msSchafSceB, "../External/Data/MS_Schafflick_PBMC/2_normalised_post_monaco.rds")




