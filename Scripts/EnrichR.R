#Check: https://cran.r-project.org/web/packages/enrichR/vignettes/enrichR.html

library(enrichR)
setEnrichrSite("Enrichr")

websiteLive <- TRUE
dbs <- listEnrichrDbs()
if (is.null(dbs)) websiteLive <- FALSE
if (websiteLive) head(dbs)

dbs <- c("GO_Molecular_Function_2015", "GO_Cellular_Component_2015", "GO_Biological_Process_2015")
if (websiteLive) {
    enriched <- enrichr(denHelGeneNames, dbs)
}

test <- enrich