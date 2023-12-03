library(ggplot2)
epdDat <- read.csv("Data/BCR_auxiliaries/End-point_dilutions.csv")
BCR_all <- read.csv("Data/BCR_database_versions/7_Subclones_included.csv")
BCR_foc <- BCR_all[which(BCR_all$CELL %in% epdDat$CELL),]
BCR_H <- BCR_foc[which(BCR_foc$LOCUS == "H"),]

BCR_H$All_mutations_H_L <- sapply(BCR_H$CELL, function(x){
    sum(BCR_foc$All_mutations[which(BCR_foc$CELL == x)])
})

BCR_match <- BCR_H[match(epdDat$CELL, BCR_H$CELL),]

identical(epdDat$CELL, BCR_match$CELL) #TRUE
#So then we can create the figure. 
ggDat <- data.frame("HL_Mutations" = BCR_match$All_mutations_H_L,
                    "EPC" = epdDat$EPD,
                    "Group" = BCR_match$Specific_antigen)
ggplot(ggDat, aes(x = HL_Mutations, y = log10(EPC), color = Group)) +
    geom_point(size = 3) + theme_bw() +
    scale_color_manual(values = c("#FF6633", "orange")) + theme(aspect.ratio=1)
ggsave("Results/Figure_4_plots/4C_mutations_vs_EPC.pdf", 
       width = 6, height = 5)

#And some statistics
cor.test(ggDat$HL_Mutations, log10(ggDat$EPC), method = "spearman")
#p-value = 5.507e-05
#Rho -0.3847808
cor.test(ggDat$HL_Mutations, log10(ggDat$EPC), method = "pearson")
#p-value = 0.004687
#cor -0.2752349

#And if separated: 
cor.test(ggDat$HL_Mutations[which(ggDat$Group == "LGI1")], 
         log10(ggDat$EPC)[which(ggDat$Group == "LGI1")], method = "spearman")
#p-value = 0.0008112
#Ro -0.3525432
cor.test(ggDat$HL_Mutations[-which(ggDat$Group == "LGI1")], 
         log10(ggDat$EPC)[-which(ggDat$Group == "LGI1")], method = "spearman")
#p-value = 0.2919
#Ro 0.2714783 
#So this is not significant at all, but it goes if anything in the opposite direction. 

#What about if the LGI1 donors are separated?
cor.test(ggDat$HL_Mutations[which(BCR_match$Sample == "1166_1")], 
         log10(ggDat$EPC)[which(BCR_match$Sample == "1166_1")], method = "spearman")
#Rho -0.209233
#p-value 0.07161
cor.test(ggDat$HL_Mutations[which(BCR_match$Sample == "1227_1")], 
         log10(ggDat$EPC)[which(BCR_match$Sample == "1227_1")], method = "spearman")
#Rho -0.5228876
#p-value 0.0811



