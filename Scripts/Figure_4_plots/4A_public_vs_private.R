#This data and analysis comers from Matthew Rayboult. Here, only the figure
#is reproduced in a suitable R format. Before this, I have transposed the data. 

covidDat <- read.csv("Data/BCR_auxiliaries/Public_vs_private/Public_SC2_ASCs_TopMatch.csv",
                     header = FALSE)
medianCovid <- median(covidDat[,1])
patDat1166 <- read.csv("Data/BCR_auxiliaries/Public_vs_private/All_1166_TopMatch.csv",
                       header = FALSE)
median1166 <- median(patDat1166[,1])
patDat1227 <- read.csv("Data/BCR_auxiliaries/Public_vs_private/All_1227_TopMatch.csv",
                       header = FALSE)
median1227 <- median(patDat1227[,1])
patDat1284 <- read.csv("Data/BCR_auxiliaries/Public_vs_private/All_1284_TopMatch.csv",
                       header = FALSE)
median1284 <- median(patDat1284[,1])

ggDat <- data.frame("Group" = c(rep("Covid", nrow(covidDat)),
                                rep(1166, nrow(patDat1166)),
                                rep(1227, nrow(patDat1227)),
                                rep(1284, nrow(patDat1284))),
                    "Vals" = c(covidDat[,1], patDat1166[,1], patDat1227[,1], patDat1284[,1]))

ggplot(ggDat, aes(x=Vals, fill=Group)) +
    geom_histogram(aes(y=after_stat(density)), 
                   color="black", position="dodge",
                   binwidth = 0.05) + theme_bw() +
    scale_fill_manual(values = c("#FAA51A", "#FFA84A", "#FF6633", "black")) +
    scale_y_continuous(expand = c(0,0,0.05,0)) +
    geom_vline(xintercept = median1166, color = "#FAA51A", 
               linetype="dotted", linewidth = 2) +
    geom_vline(xintercept = median1227, color = "#FFA84A", 
               linewidth = 2) +
    geom_vline(xintercept = median1284, color = "#FF6633", 
               linetype="dotted", linewidth = 2) +
    geom_vline(xintercept = medianCovid, color = "black", 
               linetype="dotted", linewidth = 2)
ggsave("Results/Figure_4_plots/4A_public_vs_private.pdf", width = 8, height = 5)
