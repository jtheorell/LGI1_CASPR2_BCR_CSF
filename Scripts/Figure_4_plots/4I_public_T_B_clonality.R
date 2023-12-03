#Here, we have help from a pair of great Mayo Bioinformaticians,
#Ying Li and Rohini Mopuri. They have used the Gate paper for this analysis, see
#materials and methods. 

library(ggplot2)
library(ggforce)

t_clonal_dat <- read.csv("Data/TCR_data/HC_percent_clonality.csv")

t_clonal_dat$Group <- factor(t_clonal_dat$Group, levels = c("HC", "LGI1", "CASPR2"))

wilcox.test(t_clonal_dat$percent_clonality_sample[which(t_clonal_dat$Group == "HC")],
            t_clonal_dat$percent_clonality_sample[-which(t_clonal_dat$Group == "HC")],
            exact = FALSE)
#p-value = 0.01605

dir.create("Results/Figure_4_plots/TCR_BCR_clonality")
set.seed(22)
ggplot(t_clonal_dat, aes(x = Group, y = percent_clonality_sample, 
                         shape = Group, fill = Group)) +
    geom_sina(jitter_y = FALSE, scale = "width", size = 5) +
    geom_violin(alpha = 0, scale = "width", linewidth = 1) +
    theme_bw() + scale_fill_manual(values = c("grey", "black",  "white")) +
    scale_shape_manual(values = c(21,23,23)) +
    scale_y_continuous(limits = c(0,50), expand = c(0,0))
set.seed(22)
ggsave("Results/Figure_4_plots/TCR_BCR_clonality/Healthy_vs_AE_TCR_clonality.pdf",
       height = 5, width = 4)

median_HC_TCR_Clonality <- 
    median(t_clonal_dat$percent_clonality_sample[which(t_clonal_dat$Group == "HC")])

#The B cell clonality information is not shown straight on anywhere, so that will
#Have to be retrieved. 
aeSce <- readRDS("Data/SingleCellExpFiles/csfSce_2_norm.rds")
clonTab <- table(aeSce$donor, aeSce$Clonal)
clonTab/rowSums(clonTab)
#           FALSE      TRUE
#  1166 0.6453202 0.3546798
#  1227 0.8076923 0.1923077
#  1284 0.5833333 0.4166667
library(ggplot2)
#So, those are our datapoints then: 
ggDat <- data.frame("Donor" = c("#1", "#2", "#3"),
                    "T_clonal" = 100*c(0.02479339, 0.02072539, 0.1078717),
                    "B_clonal" = 100*c(0.3546798, 0.1923077, 0.4166667),
                    "Disease" = c("LGI1", "LGI1", "CASPR2"))

ggplot(ggDat, aes(x = B_clonal, y = T_clonal, fill = Disease)) + 
    geom_point(shape = 23, size = 5) + 
    scale_x_continuous(limits = c(-1,50), expand = c(0,0)) + 
    scale_y_continuous(limits = c(0,50), expand = c(0,0)) + 
    theme_bw() + geom_vline(xintercept = 0, linetype="dotted", 
                            linewidth=1) +
    geom_hline(yintercept = median_HC_TCR_Clonality, linetype="dotted", 
                 size=1) + scale_fill_manual(values = c("white", "black"))
ggsave("Results/Figure_4_plots/TCR_BCR_clonality/T_vs_B_clonality.pdf", height = 5, width = 5)
