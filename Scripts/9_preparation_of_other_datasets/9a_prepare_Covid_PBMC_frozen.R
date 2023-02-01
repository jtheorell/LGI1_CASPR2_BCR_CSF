library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(data.table)
library(biomaRt)
#We start by running a short python script, to convert the anndata object to 
#csv files
#In terminal: 
#python3
#import scanpy
#import os
#os.chdir("/Users/jakob.theorell/Labbet/2022/220818_full_LGI1_B-cell_analysis/External/Data/Covid")
#adata = scanpy.read_h5ad("meyer_nikolic_covid_pbmc_raw.h5ad")
#adata.write_csvs('metaData')

library(zellkonverter)

#This time, the data comes from the study: 
#https://www.medrxiv.org/content/10.1101/2020.11.20.20227355v1
#Single cell profiling of COVID-19 patients: an international data resource from multiple tissues
#Chan Zuckerberg Initiative Single-Cell COVID-19 Consortia, Esteban Ballestar, 
#Donna L. Farber, Sarah Glover, Bruce Horwitz, Kerstin Meyer, Marko Nikolić, 
#Jose Ordovas-Montanes, Peter Sicovid, Alex Shalek, Niels Vandamme, 
#Linos Vandekerckhove, Roser Vento-Tormo, Alexandra Chloe Villani
#The data is from: 
#https://covid19.cog.sanger.ac.uk/submissions/release2/meyer_nikolic_covid_airway_raw.h5ad

covidSce <- readH5AD("../External/Data/Covid/meyer_nikolic_covid_pbmc_raw.h5ad",
                     reader = "R")
#Warning messages:
#    1: In value[[3L]](cond) :
#    setting 'colData' failed for '../External/Data/Covid/meyer_nikolic_covid_pbmc_raw.h5ad': HDF5.
#Links. Can't get value.
#2: In value[[3L]](cond) :
#  setting 'rowData' failed for '../External/Data/Covid/meyer_nikolic_covid_pbmc_raw.h5ad': HDF5.
#  Links. Can't get value.

#This is a craaaaaaaazy dataset with more than 400000 cells. 
#As can be seen above, we have lost the col- and rowData in the import. They are
#included here. 
colDataCovid <- read.csv("../External/Data/Covid/metaData/obs.csv")
colData(covidSce) <- DataFrame(colDataCovid)
rowDataCovid <- read.csv("../External/Data/Covid/metaData/var.csv")
rowData(covidSce) <- DataFrame(rowDataCovid)
rownames(covidSce) <- rowData(covidSce)$name

#is an internal annotation, so we can extract the plasmablasts and B-cells immediately. 
covidSce <- covidSce[,which(covidSce$annotation_broad %in% c("B", "Plasma"))]

#We also remove all pediatric patients
covidSce <- covidSce[,-which(covidSce$Group == "Paediatric")]

#This dataset does not use ensembl ids as row names (probably due to some silly
#thing in the cellranger version used. They have annotated to GRCh37, so the annotation
#is older than 10 years). But this is true also for another dataset of the
#comparison ones, so we will stick with the set of common hgnc symbols. 19029
#are overlapping between this dataset and the csfSce.
#We can however still use the information from the rowData to retrieve the 
#mitochondrial genes, etc. 

covidSce  <- getBMFeatureAnnos(covidSce, filters = "hgnc_symbol", 
                              attributes = c("ensembl_gene_id", "hgnc_symbol",
                                             "chromosome_name", "gene_biotype", 
                                             "start_position", "end_position"), 
                              dataset = "hsapiens_gene_ensembl")

#We also rename the assay from "X" to "counts". 
counts(covidSce) <- assay(covidSce, "X")
assay(covidSce, "X") <- NULL

dir.create("../External/Data/Covid", recursive = TRUE)
saveRDS(covidSce, "../External/Data/Covid/1_raw.rds")

mito <- which(rowData(covidSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(covidSce, subsets=list(Mt=mito))

colData(covidSce) <- cbind(colData(covidSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

covidSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#        665                     718                    1284                    2000 

covidSce$donor <- gsub("CV001_|-.+", "", covidSce$patient_id)

dir.create("Diagnostics/Comp_to_others/Covid", recursive = TRUE)

gridExtra::grid.arrange(
    plotColData(covidSce, x="donor", 
                y="sum", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Total count"),
    plotColData(covidSce, x="donor", 
                y="detected", colour_by = "discard") + scale_y_log10() + 
        ggtitle("Detected features"),
    plotColData(covidSce, x="donor", 
                y="subsets_Mt_percent", colour_by = "discard") + 
        ggtitle("Mito percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Comp_to_others/Covid/Exclusion_metrics_per_donor.pdf', height = 5, width = 20)
dev.off()

#NORMALIZATION
covidSce <- computeSumFactors(covidSce) 
covidSce <- logNormCounts(covidSce)

assayNames(covidSce)
#Now, for comparative reasons, we will change the label here to ASC from plasma
covidSce$cellType <-covidSce$annotation_broad
covidSce$cellType[which(covidSce$cellType == "Plasma")] <- "ASC"

table(covidSce$cellType)
#ASC     B 
#1229 19708 

#Now, we exclude the rest of the cells
saveRDS(covidSce, "../External/Data/Covid/2_normalised.rds")


