library(SingleCellExperiment)
library(scran)
library(scater)

#Here the data from step 1 is added.
csfSce <- readRDS("SingleCellExpFiles/csfSce_1_preQC.rds")

#Identification of mitochondrial genes
mito <- which(rowData(csfSce)$chromosome_name=="MT")

#At this stage, any cell with less than 100000 reads is excluded.
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))

#Now, we start with excluding cells that would otherwise screw up the statistics
#In this case, it is a matter of excluding cells that have failed to be sorted
tooFewReads <- which(stats$sum < 250000)

#And the number of lost cells is: 
length(tooFewReads)
#571. High, but acceptable.

#Which is split up as:
tooFew <- paste0(csfSce$donor, "_", csfSce$exp_number)[tooFewReads]
table(tooFew)
#0022_1 0051_2 0147_1 0253_1 1124_1 1166_1 1227_1 1227_2 1284_1 
#92    126     41     18     60     61     83     15     75 

#Which is a clear improvement; we have gained 211 cells since the last attempt,
#and then we have included more cells than then.

csfSce <- csfSce[,-tooFewReads]
unfilt2 <- csfSce
#Now the preCellQCMetrcis are calculated. 
stats <- perCellQCMetrics(csfSce, subsets=list(Mt=mito))

qc <- quickPerCellQC(stats, percent_subsets=c("subsets_Mt_percent",
                                              "altexps_ERCC_percent"), batch=csfSce$donor)

colData(unfilt2) <- cbind(colData(unfilt2), stats)
unfilt2$plate <- factor(unfilt2$plate)
unfilt2$discard <- qc$discard

unfilt2$donorNumber <- paste0(unfilt2$donor, "_", unfilt2$exp_number)

dir.create("Diagnostics")
gridExtra::grid.arrange(
    plotColData(unfilt2, x="donorNumber", y="sum", 
                colour_by="discard") + scale_y_log10() + ggtitle("Total count"),
    plotColData(unfilt2, x="donorNumber", y="detected", 
                colour_by="discard") + scale_y_log10() + ggtitle("Detected features"),
    plotColData(unfilt2, x="donorNumber", y="subsets_Mt_percent", 
                colour_by="discard") + ggtitle("Mito percent"),
    plotColData(unfilt2, x="donorNumber", y="altexps_ERCC_percent", 
                colour_by="discard") + ggtitle("ERCC percent"),
    nrow=2,
    ncol=2
)
dev.copy(pdf,'Diagnostics/Metrics_per_plate.pdf', height = 5, width = 10)
dev.off()

gridExtra::grid.arrange(
    plotColData(unfilt2, x="sum", y="subsets_Mt_percent", 
                colour_by="discard") + scale_x_log10(),
    plotColData(unfilt2, x="altexps_ERCC_percent", y="subsets_Mt_percent",
                colour_by="discard"),
    ncol=2
)
dev.copy(pdf,'Diagnostics/Mito_per_reads_and_ERCC.pdf')
dev.off()

csfSce <- csfSce[,!qc$discard]

#Here, we check how many cells that have been discarded and why
colSums(as.matrix(qc))
#low_lib_size:67  
#low_n_features: 16
#high_subsets_Mt_percent: 50 
#high_altexps_ERCC_percent: 39             
#discard: 152

#Here, we make the extended check, including the donorNumbers
qc$donorNumber <- paste0(unfilt2$donor, "_", unfilt2$exp_number)

qcTab <- do.call("rbind", lapply(split(as.data.frame(qc)[,1:5], qc$donorNumber), function(x){
    colSums(x)
}))

#In addition to this, we need to drop all A1, as it is a positive control containing
#ten cells. 
csfSce <- csfSce[,-grep("A01", colnames(csfSce))]

dim(csfSce)
#Which still leaves us with the quite impressive 3110 cells!

#Now, for integration with the flow data, we export the experiments, 
#plates and plate positions to a file. 
write.csv(colData(csfSce), "Successful_sequencing.csv")

#Now, we save the object again 
saveRDS(csfSce, file="SingleCellExpFiles/csfSce_2_postCellQC.rds")

