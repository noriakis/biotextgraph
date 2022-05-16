test_that("test basic functions produce no errors", {
    geneList <- c("2067","2068","2071","2072")

    library(ggraph)
    expect_error( wcGeneSummary(geneList), NA)
    expect_error( wcGeneSummary(geneList, ora=TRUE), NA)
    expect_error( wcGeneSummary(geneList, plotType="network"), NA)

    ## Test plotEigengeneNetworksWithWords
    library(igraph)
    mod <- returnExample()
    expect_error( plotEigengeneNetworksWithWords(mod$MEs, mod$colors), NA)
    

})