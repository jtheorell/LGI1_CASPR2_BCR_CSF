library(ggplot2)
library(ggforce)
EPDs <- read.csv("Data/BCR_auxiliaries/End-point_dilutions.csv")

allUCAs <- read.csv("Data/BCR_auxiliaries/All_with_UCA.csv")

EPDs$UCA_binder <- sapply(EPDs$CELL, function(x){
    if(x %in% allUCAs$CELL){
        locResRaw <- unique(allUCAs$Specific_UCA_BS[which(allUCAs$CELL == x)])
    } else {
        NA
    }
    
})

EPDs <- EPDs[-which(is.na(EPDs$UCA_binder)),]

EPDs$donor <- substr(EPDs$CELL, 1, 4)
EPDs$antigen <- factor(sapply(EPDs$donor, switch, 
                       "1166" = "orange", 
                       "1227" = "orange", 
                       "1284" = "#FF6633"), 
                       levels = c("orange", "#FF6633"))

#Here, we order them so that the CASPR2 events are visible. 
EPDs <- EPDs[order(EPDs$UCA_binder, EPDs$antigen),]

set.seed(10)
ggplot(EPDs, aes(x = UCA_binder,
                 y = log10(EPD))) +
    geom_sina(jitter_y = FALSE, 
              scale = "width", 
              maxwidth = 0.4, 
              size = 4, color = EPDs$antigen) +
    geom_violin(alpha = 0, scale = "width", width = 0.4, linewidth = 1, color = "black") +
    theme_bw()

ggsave("Results/Figure_4_plots/4D_UCA_binders_vs_EPC.pdf", 
       width = 4, height = 5)

wilcox.test(EPDs$EPD[which(EPDs$UCA_binder == TRUE &
                               EPDs$antigen == "orange")],
            EPDs$EPD[which(EPDs$UCA_binder == FALSE &
                               EPDs$antigen == "orange")])
#p-value 0.0001007
wilcox.test(EPDs$EPD[which(EPDs$UCA_binder == TRUE &
                               EPDs$antigen == "#FF6633")],
            EPDs$EPD[which(EPDs$UCA_binder == FALSE &
                               EPDs$antigen == "#FF6633")])
#p-value 0.1869

#And the combination
wilcox.test(EPDs$EPD[which(EPDs$UCA_binder == TRUE)],
            EPDs$EPD[which(EPDs$UCA_binder == FALSE)])
#0.001663


