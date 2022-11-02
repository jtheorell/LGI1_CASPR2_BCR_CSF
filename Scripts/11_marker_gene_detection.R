#Here, we are doing a basic marker gene detection analysis of the cells with
#known specificities. We will not take the donor information into consideration
#as the donors and the cells are so few. Instead, we will do the downsteam
#analyses based on the idea that the specific and non-specific cells represent
#discernably different "cell types"

#We are, as always, following the OSCA project. This time: 
#http://bioconductor.org/books/3.15/OSCA.basic/marker-detection.html
library(scran)
library(scater)
library(GOfuncR)
library(Homo)
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_6_ASC_only.rds")

#We now remove all cells without known specificity.
specSce <- csfSce[,-which(csfSce$Specific == "Not_tested")]

table(specSce$Specific)
#FALSE  TRUE 
#31   101

marker.info <- scoreMarkers(specSce, specSce$Specific)
marker.info
#This can be done with any number of cell types, and as we in this case
#have only two, we can go for the specific one, and take the analysis from there
chosen <- marker.info[["TRUE"]]

##Here, we change the axis to more informative names. 
identical(row.names(chosen), row.names(specSce))
#TRUE

row.names(chosen) <- paste0(seq(1,nrow(specSce)),"_",rowData(specSce)$hgnc_symbol)
specSceNewNames <- specSce
row.names(specSceNewNames) <- row.names(chosen)
#We also need to make a dataset with the same row names


#Now, as a value in 0.5 corresponds to no signal, we make it possible to order 
#according to the largest divergence from this in both directions. As the
#values here generated are made considering multiple clusters at the same time, 
#and we only have two, median, min, max, mean etc are all giving the same value.
normalizedAUC <- chosen$median.AUC-0.5
    
summary(normalizedAUC)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.201693 -0.006547  0.000000  0.004721  0.017886  0.195465


ordered <- chosen[order(abs(normalizedAUC), decreasing=TRUE),]
dir.create("Results/Marker_gene_detection")
plotExpression(specSceNewNames, features=rownames(ordered)[1:20], 
               x="Specific", colour_by="Specific")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_features_by_AUC.pdf', height = 20, width = 6)
dev.off()
dir.create("Data/Marker_gene_detection")
saveRDS(ordered, "Data/Marker_gene_detection/Genes_in_AUC_order.rds")


gene_ids <- gsub(".+_|", "", row.names(ordered)[1:20])
#Here we identify all expressed genes
expressedGenes <- rowData(csfSce)$hgnc_symbol[-which(rowSums(counts(csfSce)) == 0)]
all_genes <- unique(expressedGenes)
all_non_selected <- all_genes[-which(all_genes %in% unique(gsub(".+_|", "", row.names(ordered)[1:100])))]
input_hyper = data.frame(gene_ids, is_candidate=1)
non_input <- data.frame("gene_ids" = all_non_selected, is_candidate=0)
full_input <- rbind(input_hyper, non_input)
full_input <- full_input[-which(is.na(full_input$gene_ids)),]
res_hyper = go_enrich(full_input, n_randset=100)

stats = res_hyper[[1]]
head(stats)
by(stats, stats$ontology, head, n=3)
#This shows that there are only a very small number that actually seem to differ at all. 

top_gos_hyper = res_hyper[[1]][1:3, 'node_id']
top_gos_hyper

#"GO:0031967" "GO:0031975" "GO:0005737" "GO:1990332" "GO:1990597"
pdf("Results/Marker_gene_detection/GO_annotation_scores_top_genes.pdf")
plot_anno_scores(res_hyper, top_gos_hyper)
dev.off()

#Now, there are three potential ways of analysing this: with AUC, Cohen adjusted
#log fold changes or log fold changes directly. 
#We now go for the other main option, the Cohen values
ordered <- chosen[order(chosen$median.logFC.cohen, decreasing=TRUE),]
plotExpression(specSceNewNames, features=rownames(ordered)[1:20], 
               x="Specific", colour_by="Specific")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_features_by_Cohen.pdf', height = 20, width = 6)
dev.off()



#In all cases of a standard analyis, it is not allowed
#to directly use p-values, for the simple reason that traditionally, the cell types
#are defined by the data. In our case, however, we have two "cell types" that are
#in fact a result of a completely independent experiment, namely the specificity 
#testing. So we can calculate Wilcoxon values instead. 

wilcoxRes <- apply(logcounts(specSceNewNames), 1, function(x) 
    wilcox.test(x[which(specSceNewNames$Specific == "TRUE")],
                x[which(specSceNewNames$Specific == "FALSE")])$p.value)
write.csv(wilcoxRes, "Results/Marker_gene_detection/Wilcoxon_specific_vs_non-specific.csv")

wilcoxOrdered <- wilcoxRes[order(wilcoxRes, decreasing=FALSE)]

plotExpression(specSceNewNames, features=names(wilcoxOrdered)[1:20], 
               x="Specific", colour_by="Specific")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_features_by_Wilcox_p.pdf', height = 20, width = 6)
dev.off()
#This was a very bad strategy. 

#Now, we are going to do an even narrower analysis along the same lines, just focusing
#on the cells with known specificity and comparing the ones with a specific
#UCA and the ones without it. 
UCADat <- csfSce[,which(csfSce$Specific == "TRUE" & csfSce$Specific_UCA != "Not_tested")]

#This gives us 100 cells with a 7/3 spec/non-spec cell ratio. 
marker.info.UCA <- scoreMarkers(UCADat, UCADat$Specific_UCA)
chosen.UCA <- marker.info.UCA[["TRUE"]]
row.names(chosen.UCA) <- row.names(chosen)
UCADatNewNames <- UCADat
row.names(UCADatNewNames) <- row.names(chosen)
normalizedAUC <- chosen.UCA$median.AUC-0.5

summary(normalizedAUC)
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#-0.173234 -0.006793  0.000000  0.012721  0.026268  0.252264 

ordered <- chosen.UCA[order(abs(normalizedAUC), decreasing=TRUE),]
dir.create("Results/Marker_gene_detection")
plotExpression(UCADatNewNames, features=rownames(ordered)[1:20], 
               x="Specific_UCA", colour_by="Specific_UCA")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_UCA_features_by_AUC.pdf', height = 20, width = 6)
dev.off()

ordered <- chosen.UCA[order(chosen.UCA$median.logFC.cohen, decreasing=TRUE),]
plotExpression(UCADatNewNames, features=rownames(ordered)[1:20], 
               x="Specific_UCA", colour_by="Specific_UCA")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_UCA_features_by_Cohen.pdf', height = 20, width = 6)
dev.off()

#And we do the same for the wilcoxon as well, however unhelpful they are. 
wilcoxRes <- apply(logcounts(UCADatNewNames), 1, function(x) 
    wilcox.test(x[which(UCADatNewNames$Specific_UCA == "TRUE")],
                x[which(UCADatNewNames$Specific_UCA == "FALSE")])$p.value)
write.csv(wilcoxRes, "Results/Marker_gene_detection/Wilcoxon_UCA_specific_vs_non-specific.csv")

wilcoxOrdered <- wilcoxRes[order(wilcoxRes, decreasing=FALSE)]

plotExpression(UCADatNewNames, features=names(wilcoxOrdered)[1:20], 
               x="Specific_UCA", colour_by="Specific_UCA")
dev.copy(pdf,'Results/Marker_gene_detection/Top_20_UCA_features_by_Wilcox_p.pdf', height = 20, width = 6)
dev.off()

