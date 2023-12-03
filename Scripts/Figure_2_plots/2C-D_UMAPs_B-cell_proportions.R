library(DepecheR)
library(scater)
library(ggplot2)
library(SingleCellExperiment)
#Here, we are going to make a plot showing the B-cell compartment in all its
#complexity and overlay the specificity information. 
csfSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")

#Now, we will generate an UMAP specifically to separate the ASC from the 
#B cells 

dir.create("Results/Figure_2_plots/Flow_specific", recursive = TRUE)
flowDatLoc <- as.data.frame(t(normcounts(altExp(csfSce, "flowData")))[,c("CD38", "CD138", "CD27", "IgD", "CD20")])
#As CD138 and CD20 have more influence on the B/ASC division we increase their weight them here: 
flowDatLoc$CD138 <- flowDatLoc$CD138*1.5
flowDatLoc$CD20 <- flowDatLoc$CD20*1.5
set.seed(11)
bUmap <- uwot::umap(flowDatLoc)
dColorPlot(flowDatLoc, xYData = bUmap, plotDir = "Results/Figure_2_plots/Flow_specific/Markers")

flowDatLoc$Specific <- csfSce$Specific
flowDatLoc$donor <- csfSce$donor
flowDatLoc$cellType <- csfSce$cellType
notTestRows <- which(flowDatLoc$Specific == "Not_tested")
flowDatTested <- flowDatLoc[-notTestRows,]
bUmapTested <- bUmap[-notTestRows,]
flowDatTested$Specific <- factor(flowDatTested$Specific, 
                                 levels = c("FALSE", "LGI1", "CASPR2"))

ggDat <- cbind(flowDatTested, as.data.frame(bUmapTested))
ggplot(ggDat, aes(x = V1, y = V2, 
             fill = Specific, 
             color = Specific)) +
    geom_point(size = 4) + scale_fill_manual(values = c("black", "orange", "#FF6633"))+
    scale_color_manual(values = c("black", "orange", "#FF6633"))+
    theme_bw() + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)  
ggsave("Results/Figure_2_plots/Flow_specific/Specific_spec_only_C_and_L_for_main.pdf", width = 6, height = 4)

ggplot(ggDat, aes(x = V1, y = V2, fill = Specific, shape = donor)) +
    geom_point(size = 4) + scale_fill_manual(values = c("black", "orange", "#FF6633"))+
    theme_bw() + scale_shape_manual(values = c(24, 22, 23)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)
ggsave("Results/Figure_2_plots/Flow_specific/Specific_spec_only_C_and_L_for_supp.pdf", width = 6, height = 4)


#Now, we plot cell type vs specificity

ggplot(ggDat, aes(x = V1, y = V2, fill = cellType, color = cellType)) +
    geom_point(size = 4) + scale_fill_manual(values = c("#2388DD","#BB3322"))+
    scale_color_manual(values = c("#2388DD","#BB3322"))+
    theme_bw() + theme(panel.grid.major = element_blank(), 
                       panel.grid.minor = element_blank(),
                       aspect.ratio=1)
ggsave("Results/Figure_2_plots/Flow_specific/Cell_type_spec_for_main.pdf", 
       width = 6, height = 4)
ggplot(ggDat, aes(x = V1, y = V2, fill = cellType, shape = donor)) +
    geom_point(size = 4) + scale_fill_manual(values = c("#2388DD","#BB3322"))+
    theme_bw() + scale_shape_manual(values = c(24, 22, 23)) + 
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(),
          aspect.ratio=1)
ggsave("Results/Figure_2_plots/Flow_specific/Cell_type_spec_for_supp.pdf", 
       width = 6, height = 4)

#Now, we are going to make bar graphs with this informaiton. 

spec_cell_df <- as.data.frame(table(as.character(flowDatTested$Specific), flowDatTested$cellType))
colnames(spec_cell_df) <- c("Specific", "Cell_type", "Freq")
spec_cell_df$B_subtype <- factor(spec_cell_df$Cell_type, c("ASC", "B"))
spec_cell_df$Specific <- factor(spec_cell_df$Specific, c("FALSE",
                                                         "LGI1",
                                                         "CASPR2"))
ggplot(spec_cell_df, aes(x=Specific, y=Freq, fill=Cell_type)) +
    geom_bar(stat="identity") + theme_bw() + 
    scale_y_continuous(limits = c(0,100),
                       expand = expansion(mult=c(0,0))) +
    scale_fill_manual(values = c("#2388DD","#BB3322"))
ggsave("Results/Figure_2_plots/Specificity_vs_cell_type_C_and_L.pdf", width = 5, height = 6)

spec_cell_df
#  Specific Cell_type Freq B_subtype
#1     CASPR2       ASC   31       ASC
#2      FALSE       ASC   17       ASC
#3       LGI1       ASC   83       ASC
#4     CASPR2         B    2         B
#5      FALSE         B   10         B
#6       LGI1         B    2         B

#We also make a per-donor version
spec_cell_df <- as.data.frame(table(paste0(flowDatTested$donor, "_", flowDatTested$Specific), flowDatTested$cellType))
colnames(spec_cell_df) <- c("Spec_don", "Cell_type", "Freq")
spec_cell_df$B_subtype <- factor(spec_cell_df$Cell_type, c("ASC", "B"))
spec_cell_df$Spec_don <- factor(spec_cell_df$Spec_don, c("1166_FALSE",
                                                         "1166_LGI1",
                                                         "1227_FALSE",
                                                         "1227_LGI1",
                                                         "1284_FALSE",
                                                         "1284_CASPR2"))
ggplot(spec_cell_df, aes(x=Spec_don, y=Freq, fill=Cell_type)) +
    geom_bar(stat="identity") + theme_bw() + 
    scale_y_continuous(limits = c(0,100),
                       expand = expansion(mult=c(0,0))) +
    scale_fill_manual(values = c("#2388DD","#BB3322"))
ggsave("Results/Figure_2_plots/Specificity_vs_cell_type_C_and_L_donors_separated.pdf", width = 6, height = 6)

spec_cell_df
#      Spec_don Cell_type Freq B_subtype
#1   1166_FALSE       ASC   14       ASC
#2    1166_LGI1       ASC   72       ASC
#3   1227_FALSE       ASC    0       ASC
#4    1227_LGI1       ASC   11       ASC
#5  1284_CASPR2       ASC   31       ASC
#6   1284_FALSE       ASC    3       ASC
#7   1166_FALSE         B    3         B
#8    1166_LGI1         B    2         B
#9   1227_FALSE         B    2         B
#10   1227_LGI1         B    0         B
#11 1284_CASPR2         B    2         B
#12  1284_FALSE         B    5         B

#So only a tiny fraction in the specific compatment are non-ASC, and the donor-
#to donor pattern is very similar, with the smallest sample, from 1227, being the 
#most exgtreme with no non-specifics being ASC and no specifics being B. 

#We will also make two a fisher tests here, one for LGI1 and one for CASPR2. 
fisherDf <- data.frame("LGI1" = c(83,2),
                       "Negative" = c(14,5))
row.names(fisherDf) <- c("ASC", "B")

fisherDf
#    LGI1 Negative
#ASC   83       14
#B      2        5

fisher.test(fisherDf)
#p-value = 0.002065, i.e. the B-cell frequency is "significantly" lower in the specific group. 
#THis is of course treating each antibody as a spearate observation, which is a violation
#of statistical rules. 

#We will also make two a fisher tests here, one for LGI1 and one for CASPR2. 
fisherDf <- data.frame("CASPR2" = c(31,2),
                       "Negative" = c(3,5))
row.names(fisherDf) <- c("ASC", "B")

fisherDf
#    CASPR2 Negative
#ASC     31        3
#B        2        5

fisher.test(fisherDf)
#p-value = 0.001357, so also here, there is a difference with the same caution as above. 

#And now, finally, separation of the LGI1 donors, for the supplement. 
#1166:
fisherDf <- data.frame("LGI1" = c(72,2),
                       "Negative" = c(14,3))
row.names(fisherDf) <- c("ASC", "B")

fisherDf
#    LGI1 Negative
#ASC   72       14
#B      2        3

fisher.test(fisherDf)
#p-value = 0.04341

#1227:
fisherDf <- data.frame("LGI1" = c(11,0),
                       "Negative" = c(0,2))
row.names(fisherDf) <- c("ASC", "B")

fisherDf
#    LGI1 Negative
#ASC   11       0
#B      0       2

fisher.test(fisherDf)
#p-value = 0.01282

