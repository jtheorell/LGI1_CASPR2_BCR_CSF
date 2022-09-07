library(ensembldb)
library(tximport)
library(SingleCellExperiment)
library(scran)
library(scater)

#I use the file Homo_sapiens.GRCh38.98.gtf.gz and combine it with the ERCC information
#In the following way.
#If not run before: 
#DB <- ensDbFromGtf(gtf="~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/Homo_sapiens.GRCh38.98.gtf.gz")

EDB <- EnsDb("~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/Homo_sapiens.GRCh38.98.sqlite")
tx2geneId <- transcripts(EDB, columns=c("tx_id", "gene_id"), return.type="DataFrame")
#Now add the ERCC spike-ins
ERCCgtf <- read.table('~/Labbet/Saved_transcriptome_files/HomoS_GRCh38.98_cdna_ncrna_ERCC/ERCC92.gtf.gz', 
                      header = FALSE, sep = '\t')
ERCCNamed <- data.frame("tx_id"=as.character(ERCCgtf[,1]), 
                        "gene_id"=as.character(ERCCgtf[,1]), 
                        stringsAsFactors = FALSE)
tx2geneComplete <- rbind(tx2geneId, ERCCNamed)

#And now we import this.
dirs <- list.files("Salmon_results")
fullDirs <- file.path("Salmon_results", dirs, "quant.sf") 

txiFileGenes <- tximport(files=fullDirs, type = "salmon",  
                         tx2gene = tx2geneComplete,
                         ignoreTxVersion=TRUE)

dirs <- list.files("Salmon_results")
fullDirs <- file.path("Salmon_results", dirs, "quant.sf") 

txiFileGenes <- tximport(files=fullDirs, type = "salmon",  
                         tx2gene = tx2geneComplete,
                         ignoreTxVersion=TRUE)

csfSce <- SingleCellExperiment(assays = 
                                   list(counts = txiFileGenes$counts, 
                                        cpm = txiFileGenes$abundance))

colnames(csfSce) <- dirs

#Add some meta information embedded in the cell names
csfSce$donor <- gsub("|_._._...", "\\1", colnames(csfSce))
csfSce$exp_number <- gsub("...._|_._...", "\\1", colnames(csfSce))
csfSce$plate <- gsub("...._._|_...", "\\1", colnames(csfSce))
csfSce$plate_position <- gsub("...._._._|", "\\1", colnames(csfSce))

#Here, rowData is added, to get information not only about the gene symbol, but
#also the chromosome position, gene biotype, etc. 
library(scater)
csfSce <- getBMFeatureAnnos(csfSce, filters = "ensembl_gene_id", 
                            attributes = c("ensembl_gene_id", "hgnc_symbol",
                                           "chromosome_name", "gene_biotype", 
                                           "start_position", "end_position"), 
                            dataset = "hsapiens_gene_ensembl")

#Now, the ERCC are defined as such
is.spike <- grepl("ERCC", rownames(csfSce))
csfSce <- splitAltExps(csfSce, ifelse(is.spike,"ERCC", "gene"), ref = "gene")

#Now, we exclude all cells that are not included among the BCRs included in
#this analysis. 
BCR_all <- read.csv("Data/BCR_database_versions/6_Specificity_included.csv")
csfSceB <- csfSce[,which(colnames(csfSce) %in% BCR_all$CELL)]

#Here, the file is saved
dir.create("SingleCellExpFiles")

#Here, we exclude all cells that do not have a complete BCR
#Three cells are not represented among the transcriptome-containing cells, 
#but these are likely to disappear downstream in the filtering steps to come anyway.
saveRDS(csfSceB, file="Data/SingleCellExpFiles/csfSce_1_preQC.rds")
