#Here, we include a new set of mutations analyses, as it has turned out
#that the standard setting does not show the CDR3 mutations. 
library(shazam)
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")

#With this analysis, we do not only include the FWR4 and the CDR3, but we also
#split it all up into the individual regions. 
BCR_full_mutations <- observedMutations(BCR_all, sequenceColumn = "SEQUENCE_IMGT",
                                   germlineColumn = "GERMLINE_IMGT",
                                   regionDefinition = IMGT_VDJ_BY_REGIONS,
                                   juncLengthColumn = "JUNCTION_LENGTH",
                                   nproc = 7)
BCR_mutations_red <- BCR_full_mutations[,-grep("_cdr_|_fwr_", colnames(BCR_full_mutations))]

BCR_mutations_red$All_mutations <- apply(BCR_mutations_red[,grep("cdr|fwr", colnames(BCR_mutations_red))], 1, sum)
BCR_mutations_red$Non_silent_mutations <- apply(BCR_mutations_red[,grep("_r", colnames(BCR_mutations_red))], 1, sum)

BCR_mutations_red$CDR_mutations <- apply(BCR_mutations_red[,grep("cdr", colnames(BCR_mutations_red))], 1, sum)
BCR_mutations_red$Non_silent_CDR_mutations <- apply(BCR_mutations_red[,grep("cdr._r", colnames(BCR_mutations_red))], 1, sum)


#And this is saved
write.csv(BCR_mutations_red, 
          "Data/BCR_database_versions/8_new_mutations_added.csv")