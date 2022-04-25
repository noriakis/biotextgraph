# wcGeneSummary

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/noriakis/wcGeneSummary/workflows/R-CMD-check/badge.svg)](https://github.com/noriakis/wcGeneSummary/actions)
  <!-- badges: end -->

Make wordcloud or a plot of correlation network of gene set from RefSeq description using R libraries [`GeneSummary`](https://bioconductor.org/packages/release/data/annotation/html/GeneSummary.html), [`tm`](https://www.jstatsoft.org/article/view/v025i05) and `wordcloud`. Input is gene list with the type `ENTREZID`. I think it is useful when the enrichment analysis returns no significant results for a set of genes. The idea for using RefSeq description in [`simplyfyEnrichment`](https://github.com/jokergoo/simplifyEnrichment) is already raised [here](https://github.com/jokergoo/simplifyEnrichment/issues/56).

### Installation
```R
devtools::install_github("noriakis/wcGeneSummary")
library(wcGeneSummary)
```

### Example of ERCC
```R
erccs <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")
entrezID = AnnotationDbi::select(org.Hs.eg.db, keys = erccs, columns = c("ENTREZID"), keytype = "SYMBOL")$ENTREZID
gwc <- wcGeneSummary(entrezID, excludeFreq = 5000, max.words=200, random.order=FALSE,
                     colors=palettetown::pokepal(150), shape="circle", rot.per=0.4)
# ggsave("erccWc.png", gwc$wc, width=8, height=8)
```
<img src="https://github.com/noriakis/software/blob/main/images/erccWc.png?raw=true" width="800px">

### Example of CXCL
```R
cxcls <- c()
for (i in c(1,2,3,5,6,8,9,10,11,12,13,14,16)){
    cxcls <- c(cxcls, paste0("CXCL",i))
}
entrezID = AnnotationDbi::select(org.Hs.eg.db, keys=cxcls, columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID
gwc <- wcGeneSummary(entrezID, excludeFreq=14000,
                     madeUpper=c("dna","rna",tolower(keys(org.Hs.eg.db, keytype="SYMBOL"))),
                     max.words=200, random.order=FALSE,
                     colors=palettetown::pokepal(151), shape="circle", rot.per=0.4)
# ggsave("cxclWc.png", gwc$wc, width=8, height=8)
```
<img src="https://github.com/noriakis/software/blob/main/images/cxclWc.png?raw=true" width="800px">

### Example of CCL (correlation network)
```R
library(wcGeneSummary)
library(org.Hs.eg.db)
library(ggraph)
ccls <- c()
for (i in c(1,2,3,4,5,6,7,8,9)){
    ccls <- c(ccls, paste0("CCL",i))
}
entrezID = AnnotationDbi::select(org.Hs.eg.db, keys=ccls, columns=c("ENTREZID"), keytype="SYMBOL")$ENTREZID
cclNet <- wcGeneSummary(entrezID, plotType="network",
                        layout="nicely",
                        madeUpper=c("dna","rna",tolower(keys(org.Hs.eg.db, keytype="SYMBOL"))),
                        numWords = 15, excludeFreq = 5000, edgeLink=FALSE, showLegend=FALSE)
cclNetTrans <- cclNet$net + theme(plot.background = element_rect(fill = "transparent",colour = NA))
# ggsave(file="cclNet.png", cclNetTrans, width=7, height=7, bg="transparent")
```

<img src="https://github.com/noriakis/software/blob/main/images/cclCorNetNicely.png?raw=true" width="800px">

### Example of a corelation network of words in the pathway
```R
library(org.Hs.eg.db)
keggPathways <- org.Hs.egPATH2EG
mappedKeys <- mappedkeys(keggPathways)
keggList <- as.list(keggPathways[mappedKeys])
## Hepatitis C
hCNet <- wcGeneSummary(keggList$`05160`, plotType="network",
                        layout="nicely", corThresh = 0.2,
                        madeUpper=c("dna","rna",tolower(keys(org.Hs.eg.db, keytype="SYMBOL"))),
                        numWords = 30, excludeFreq = 5000, colorText=TRUE,
                        edgeLink=FALSE, showLegend=FALSE)
hCNetTrans <- hCNet$net + theme(plot.background = element_rect(fill = "transparent",colour = NA))
# ggsave(file="hCNet.png", hCNetTrans, width=12, height=12, bg="transparent")
```
<img src="https://github.com/noriakis/software/blob/main/images/hCNet.png?raw=true" width="800px">

### Example of annotating dendrogram of gene cluster by words

The example of annotating dendrogram by pyramid plots of word counts is shown. In this example, `WGCNA` was used to cluster the gene expression values. Module eigengenes are further clustered by `pvclust`. `getWordsOnDendro` can be used to obtain the `patchworkGrob` list. Grobs can be plotted on dendrogram plot using the `annotation_custom`. Plotting is not functionalized due to positioning requirements.

```R
library(wcGeneSummary)
library(org.Hs.eg.db)
library(dendextend)
library(pvclust)
library(ggplot2)
library(grid)
library(gridExtra)
library(magrittr)

## An example to annotate module eigengenes in WGCNA
## Perform pvclust on ME data.frame
## block-wise module identification was performed beforehand
result <- pvclust(bwmod$MEs[,1:7], method.dist="cor", method.hclust="average", nboot=10)

## Make the dendrogram
dhc <- result %>%
    as.dendrogram() %>%
    hang.dendrogram()

## Make named vector
geneVec <- paste0("ME", bwmod$colors)
names(geneVec) <- names(bwmod$colors)

## Get pyramid plot list using the function.
## It takes time when geneNumLimit is large.
grobList <- getWordsOnDendro(dhc, geneVec, numberOfWords = 10, geneNumLimit = 1000,
                             geneVecType = "ENSEMBL", excludeFreq = 10000)

## Plot dendrogram ggplot, using the pvclust p-values
dendroPlot <- dhc %>% pvclust_show_signif_gradient(result) %>% ggplot() 

## Plot the grob on dendrogram using annotation_custom.
## If border is TRUE, border line is drawn using grid.rect.
border <- TRUE
for (gr in grobList){
    if (border){
        addPlot <- ggplotify::as.grob(function(){
            grid.arrange(gr$plot)
            grid.rect(width = .98, height = .98, gp = gpar(lwd = 1, col = "black", fill = NA))
        })
    } else {
        addPlot <- gr$plot
    }
    dendroPlot <- dendroPlot +
        annotation_custom(addPlot, xmin=gr$xmin, xmax=gr$xmax, ymin=gr$height+0.005, ymax=gr$heightup-0.005)
}

dendroPlot
```
<img src="https://github.com/noriakis/software/blob/main/images/plotDendro.png?raw=true" width="800px">

### The other examples

- [Example annotating gene clusters using WGCNA](https://noriakis.github.io/software/wcGeneSummary/)
- [Example in Bayesian network analysis](https://github.com/noriakis/compare_sign)

### References
- [Zuguang Gu. GeneSummary](https://doi.org/doi:10.18129/B9.bioc.GeneSummary)
- Feinerer I, Hornik K, Meyer D. Text Mining Infrastructure in R. J Stat Softw 2008;25:1–54. doi:10.18637/jss.v025.i05
- [Ian Fellows. wordcloud: Word Clouds. R package version 2.6. 2018.](https://cran.r-project.org/web/packages/wordcloud/index.html)
- [Lucas T. palettetown: Pokemon themed colour schemes for R.](https://github.com/timcdlucas/palettetown)