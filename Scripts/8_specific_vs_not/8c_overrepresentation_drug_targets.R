#We are interested in clarifying whether any of the drug targets of interest are enriched
#when comparing our specific cells, i.e. the ASC, to the non-specifics, i.e. the B-cells. 
#To do this, we are going to use a widely establshed method, namely over-representation
#analysis. It is implemented in the clusterProfiler package, in a function there called 
#enricher. However, there are formatting problems with this function. The essential
#function here is however investigations of the hypergeometric distribution, which 
#can be done either using the function phyper with the lower tail = FALSE or a one-tailed
#Fisher exact test, focusing on the upper tail. Introductions to gene over-representation
#analysis can be found here: 
#https://alexslemonade.github.io/refinebio-examples/03-rnaseq/pathway-analysis_rnaseq_01_ora.html#4_Over-Representation_Analysis_with_clusterProfiler_-_RNA-seq
#And here, section 5.1:
#https://yulab-smu.top/biomedical-knowledge-mining-book/enrichment-overview.html
#And here: 
#https://yulab-smu.top/biomedical-knowledge-mining-book/universal-api.html

#To do this, we need to define a subset of genes that are the top upregulated ones
#Then, we need to define a gene universe, i.e. the genetic background to compare
#our hit genes to. The third important component is of course the gene sets that
#we are interested in. As most of the drugs we are after are biologicals, we
#have had to resort to not checking for known pharmoacological targets, but rather
#the interactome for each individual drug target. These interactomes have been 
#identified in the STRING database: 
#https://string-db.org/cgi/input?sessionId=bSQFytfwqJ0e&input_page_show_search=on. 
#Here, we have searched for each individual protein listed below, changed two settings: 
#the minimum required interaction score has been increased to a confidence of >0.7, 
#and the number of interactors have been increased to a max of 20. 

#The list of targets is: 
#Now, we are going to make a list of all drugs that we think are of interest
#in this context. 
#General: prednisolone, dexamethasone
#CD27: varlilumab (agonist)
#CD38: Daratumumab
#CXCR4: Plerixafor, BL-8040, LY2510924, ulocuplumab, CX-01 (dociparstat sodium)
#CXCR5: can only find CAR-T cells
#MS4A1: rituximab, ofatumumab, ocrelizumab, Ublituximab, Obinutuzumab, Ocaratuzumab, Veltuzumab, Tositumomab
#CD19: tafasitamab, loncastuximab tesirine, Inebiluzimab, Obexelimab
#SDC1: Indatuximab ravtansine
#CD22: epratuzumab, inotuzumab ozogamicin
#BAFF/BAFFR: Belimumab, Ianalumab, Telitacicept, Atacicept (not efficient in MS) (2 last are TACI-Ig-fusion proteins),
#APRIL/TACI: See BAFF above
#IL-21/IL-21R: NN 9828, NNC0114-0005, One called Ab-01, but seems far from clinical
#IL-6R: tocilizumab, satriluzimab
#IL6: siltuximab
#IL-10 unclear
#IL-10R unclear
#IFN-g HuZAF/unclear
#IGNGR unclear
#GMCSF Gimsilumab, lenzilumab, namilumab, and otilimab
#GMCSFR Mavrilimumab
#a4B1 integrin (ITGA4 is the target): natalizumab
#CD40 (agonists): dacetuzumab, CDX-1140, APX005M, CP-870,893, ADC-1013

#These are also included in small molecular screening efforts. 
#S1PR (all five): fingolimod, siponimod
#proteosome pathways (chosen target PSMC4): bortezomib, carfilzomib
#BTK: ibrutinib, Tolebrutinib, Evobrutinib
#JAK-STAT (JAK1/2): baricitinib , tofacitinib, upadacitinib 

#Taken out, due to lack of clinically used drugs for these disorders: "IL10", "IL10RA","IL10RB", "IFNG",
#"IFNGR1","IFNGR2", "TNF"  (infliximab, adalimumab, etanercept, golimumab, and certolizumab),
#"TNFAR" infliximab, adalimumab, etanercept, golimumab, and certolizumab


#The names that might be hard to get from the list are: 
#MS4A1 = CD20
#SDC1 = CD138
#TNFSF13B = BAFF
#TNFRSF13C = BAFFR
#ANP32B = APRIL
#TNFRSF13B = TACI
#TNFRSF1A = TNF receptor 1a
#TNFRSF1B = TNF receptor 1b
#CSF2 = GM-CSF
#CSF2RA = GM-CSF-RA
#CSF2RB = GM-CSF-RB

res <- readRDS("Results/Specific_vs_not/EdgeR_result.rds")
#Now, we identify the gene names of the genes that are most significantly upregulated
resTab <- res$table
topGenes <- row.names(resTab)[which(resTab$logFC > 2 & resTab$PValue < 0.01)]
length(topGenes)
#10

bCellUniverse <- row.names(resTab)
length(bCellUniverse)
#9976

#The drug tagets are now imported
drugTargetFileNameList <- list.files("Data/STRING_data", full.names = TRUE, pattern = ".tsv")
drugTargetFileList <- lapply(drugTargetFileNameList, read.table)

#Now, as column 1 and 2 list interactions, we will identify all the unique ones
#in both, and then serve these as the interesting ones. 
drugTargetList <- lapply(drugTargetFileList, function(x){
    unique(unlist(x[,1:2]))
})

names(drugTargetList) <- gsub("Data/STRING_data/|_string_interactions_short.tsv", 
                              "", drugTargetFileNameList)
length(unique(unlist(drugTargetList)))
#327

##########################
#Fisher-based over-representation analysis
##########################
#These will be needed for all below

nonTargUni <- bCellUniverse[-which(bCellUniverse %in% unlist(drugTargetList))]
length(nonTargUni)
#9780

nonTargetHit <- length(which(topGenes %in% nonTargUni))
nonTargetHit
#10

#################
#Direct targets
#################

#To check this, we need to start by changing CD138 and CD20 to their gene names
targHit <- length(which(topGenes %in% names(drugTargetList)))

fisherDf <- data.frame("Non_target" = c(length(nonTargUni)-nonTargetHit, nonTargetHit),
                       "Target" = c(length(drugTargetList)-targHit, targHit))
row.names(fisherDf) <- c("Non_hit", "Hit")

fisherDf
#        Non_target Target
#Non_hit       9770     27
#Hit             10      0

fisher.test(fisherDf, alternative = "greater")
#p-value = 1, so no hits there. 

#Now, what about all genes associated to all targets?

targHit <- length(which(topGenes %in% unlist(drugTargetList)))
targHit
#0

#So that is that. 
