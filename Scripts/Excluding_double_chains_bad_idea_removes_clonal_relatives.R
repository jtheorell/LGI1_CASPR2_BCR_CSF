#Now, as an extra precaution against chains unlikely in use, we
#check the cells with two light chains, and if one of them is among the 10% 
#shortest chains of that kind, it is excluded. If both are, then the longer
#one is chosen. 
shortK <- quantile(nchar(bracerFile1H$SEQUENCE_VDJ_UNGAPPED[which(bracerFile1H$LOCUS == "K")]), 0.1)
shortK
#277.9
hist(nchar(bracerFile1H$SEQUENCE_VDJ_UNGAPPED[which(bracerFile1H$LOCUS == "K")]), breaks = 100)
abline(v = shortK)

shortL <- quantile(nchar(bracerFile1H$SEQUENCE_VDJ_UNGAPPED[which(bracerFile1H$LOCUS == "L")]), 0.1)
shortL
#311.4
hist(nchar(bracerFile1H$SEQUENCE_VDJ_UNGAPPED[which(bracerFile1H$LOCUS == "L")]), breaks = 100)
abline(v = shortL)

kappaCells <- bracerFile1H$CELL[which(bracerFile1H$LOCUS == "K")]
lambdaCells <- bracerFile1H$CELL[which(bracerFile1H$LOCUS == "L")]
twoLight <- kappaCells[which(kappaCells %in% lambdaCells)]

bracerFile1HDouble <- bracerFile1H[which(bracerFile1H$CELL %in% twoLight),]
bracerFile1HDouble <- 
    do.call("rbind", lapply(twoLight, function(x){
        locCell <- bracerFile1HDouble[which(bracerFile1HDouble$CELL == x),]
        lengthK <- nchar(locCell$SEQUENCE_VDJ_UNGAPPED[which(locCell$LOCUS == "K")])
        lengthL <- nchar(locCell$SEQUENCE_VDJ_UNGAPPED[which(locCell$LOCUS == "L")])
        anyShort <- c(FALSE, FALSE)
        if(lengthK < shortK){anyShort[1] <- TRUE}
        if(lengthL < shortL){anyShort[1] <- TRUE}
        if(any(anyShort)){
            shortest <- c("K", "L")[which.min(c(lengthK, lengthL))]
            locCellLong <- locCell[-which(locCell$LOCUS == shortest),]
        } else {
            locCell
        }
        
    }))

#This makes us lose 15 chains. 


#This was used in the selection sheet, but not in the original, so it cannot be used:; 
#Are there any cells that have both a too short K and a too shot L chain?
length(which(short_chains$K %in% short_chains$L))
#2
#So we will need to do this separately for the double cells. 
doubleShort <- short_chains[which(short_chains$K %in% short_chains$L)]

shortAll <- unique(unlist(short_chains))
shortNonDouble <- 
    shortAll[-which(shortAll %in% BCR_data$CELL[which(BCR_data$light_type == "double")])]

realShort <- c(doubleShort, shortNonDouble)
