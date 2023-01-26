library(SingleCellExperiment)
library(scater)
library(Matrix)
library(scran)
library(data.table)

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

#This is a craaaaaaaazy dataset with more than 400000 cells. But there
#is an internal annotation, so we can extract the plasmablasts immediately. 
covidSce <- covidSce[,which(covidSce$annotation_detailed == "Plasmablasts")]

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

dir.create("Data/Comp_to_others/Covid/SingleCellExpFiles", recursive = TRUE)
saveRDS(covidSce, "Data/Comp_to_others/Covid/SingleCellExpFiles/0_raw.rds")

mito <- which(rowData(covidSce)$chromosome_name=="MT")

stats <- perCellQCMetrics(covidSce, subsets=list(Mt=mito))

colData(covidSce) <- cbind(colData(covidSce), stats)

reasons <- perCellQCFilters(stats, 
                            sub.fields=c("subsets_Mt_percent"))

covidSce$discard <- reasons@listData$discard
colSums(as.matrix(reasons))

#low_lib_size          low_n_features high_subsets_Mt_percent                 discard 
#4                       5                       21                       24

covidSce$donor <- gsub("CV001_|-.+", "", colnames(covidSce))

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

#Actually, this quite expressedly shows that this data already has gone through an exlusion
#event, so no cells will be excluded. 

#NORMALIZATION
covidSce <- computeSumFactors(covidSce) 
covidSce <- logNormCounts(covidSce)

assayNames(covidSce)

#Now, we exclude the rest of the cells
saveRDS(covidSce, "Data/Comp_to_others/Covid/SingleCellExpFiles/2_plasmablasts.rds")


