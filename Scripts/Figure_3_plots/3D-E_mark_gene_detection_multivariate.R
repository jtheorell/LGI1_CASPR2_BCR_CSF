library(scran)
library(scater)
library(SingleCellExperiment)
library(ggplot2)
library(edgeR)
library(pheatmap)
#Here, we are going to try to model both cell type (i.e. B vs ASC) and specificity
#at the same time. We rely heavily on: 
#http://bioconductor.org/books/3.13/OSCA.multisample/multi-sample-comparisons.html
aeSce <- readRDS("Data/SingleCellExpFiles/4_all_spec_with_LLPC_info.rds")

#Now, we will remove all without a HGNC symbol and then change the names to them
dupHgnc <- rowData(aeSce)$hgnc_symbol[which(duplicated(rowData(aeSce)$hgnc_symbol))]

table(dupHgnc)
#      ABCF2    AHRR   ATXN7    GGT1   PINX1 POLR2J3    TBCE 
#606       1       1       1       1       1       1       1
aeSceHgnc <- aeSce[-which(rowData(aeSce)$hgnc_symbol %in% dupHgnc),]
row.names(aeSceHgnc) <- rowData(aeSceHgnc)$hgnc_symbol

##################
#EdgeR analysis
##################
current <- aggregateAcrossCells(SingleCellExperiment(assays = list("counts" = counts(aeSceHgnc))), 
                                ids=colData(aeSceHgnc)[,c("Specific", "cellType", "donor")])
current
# Creating up a DGEList object for use in edgeR:
y <- DGEList(counts(current), samples=colData(current))
y$samples$ncells
#14  3  3  2  5 72 11 31  2  2

#Here, we are extremely pragmatic. We have ten "samples" of which five have two
#or three cells in them, and we want to include all cells, so we will run this
#extremely low bulk as it is, and the results will be thereafter. 

design <- model.matrix(~factor(cellType)+factor(Specific), y$samples)

keep <- filterByExpr(y, design = design)
y <- y[keep,]
summary(keep)

#   Mode   FALSE    TRUE 
#logical   10454    8766

y <- calcNormFactors(y)

#design
y <- estimateDisp(y, design)
#summary(y$trended.dispersion)
fit <- glmQLFit(y, design, robust=TRUE)
#summary(fit$var.prior)
res <- glmQLFTest(fit, coef=ncol(design))

resSig <- res$table[which(abs(res$table$logFC) > 3 & res$table$PValue < 0.01),]

#So these 25 genes are plotted. First as a volcano. 
plotDat <- res$table
plotDat$significant <- FALSE
plotDat$significant[which(row.names(plotDat) %in% row.names(resSig))] <- TRUE
dir.create("Results/Figure_3_plots/3D-E_transcriptome")
ggplot(plotDat, aes(x = logFC, y = -log10(PValue), color = significant)) + 
    geom_point(size = 4) + theme_bw() + theme(aspect.ratio=1) + 
    scale_color_manual(values = c("grey", "red")) + 
    scale_x_continuous(limits = c(-15, 15), expand = c(0,0))
ggsave(paste0("Results/Figure_3_plots/3D-E_transcriptome/Multivariate_volcano_hit_colors.pdf"), plot = p,
       width = 5, height = 6)

#And now, a heatmap. 

aeSceHgncSig <- aeSceHgnc[which(row.names(aeSceHgnc) %in% row.names(resSig)),]
aeSceHgncSig$donSpecCell <- paste0(aeSceHgncSig$cellType, "_", 
                                   aeSceHgncSig$Specific, "_", 
                                   aeSceHgncSig$donor)
aveLogCountMat <- do.call("cbind", lapply(split(data.frame(t(logcounts(aeSceHgncSig))), 
                                                          aeSceHgncSig$donSpecCell), colMeans))

aveLogCountDfNorm <- data.frame(t(apply(aveLogCountMat, 1, function(x) x/max(x))))

#And finally, we reorder the columns. 
aveLogCountDfNormOrd <- aveLogCountDfNorm[,c(grep("B", colnames(aveLogCountDfNorm)),
                                            grep("ASC", colnames(aveLogCountDfNorm)))]

pdf("Results/Figure_3_plots/3D-E_transcriptome/Heatmap_of_diff_genes.pdf", height = 7, width = 5)
pheatmap(aveLogCountDfNormOrd,
         cluster_cols = FALSE,
         color = gray.colors(50, start =1, end = 0))
dev.off()

#And finally, we run a GO enrichment analysis here. 
library(clusterProfiler)
library(org.Hs.eg.db)
gene_list <- res$table$logFC
names(gene_list) <- row.names(res)

gene_list <- sort(gene_list, decreasing = TRUE)
set.seed(101)
gse <- gseGO(geneList=gene_list, 
             ont ="ALL", 
             keyType = "SYMBOL", 
             minGSSize = 10, 
             maxGSSize = 800, 
             pvalueCutoff = 0.05, 
             verbose = TRUE, 
             OrgDb = org.Hs.eg.db, 
             eps = 0)

dotplot(gse, split=".sign") + facet_grid(.~.sign)
ggsave("Results/Figure_3_plots/3D-E_transcriptome/Gene_set_enrichment_dotplot.pdf", width = 5, height = 5)


#Now, we also want to investigate the over-expressed genes in some more detail
resSigOrd <- resSig[order(resSig$PValue),]
