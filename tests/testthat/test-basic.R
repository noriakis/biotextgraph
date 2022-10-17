test_that("test basic functions produce no errors", {
    geneList <- c("CXCL10","CXCL9")

    library(ggraph)
    expect_error( wcGeneSummary(geneList), NA)
    expect_error( wcGeneSummary(geneList, ora=TRUE), NA)
    expect_error( wcGeneSummary(geneList, plotType="network"), NA)
    expect_error( wcGeneSummary(geneList,
        plotType="network", genePlot=TRUE), NA)
    
    ## Test plotEigengeneNetworksWithWords
    library(igraph)
    mod <- returnExample()
    expect_error( plotEigengeneNetworksWithWords(mod$MEs, mod$colors), NA)
    

})