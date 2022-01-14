# wcGeneSummary

Make wordcloud or a plot of correlation network of gene set from RefSeq description using R libraries [`GeneSummary`](https://bioconductor.org/packages/release/data/annotation/html/GeneSummary.html), `tm` and `wordcloud`. Input is gene list with the type `ENTREZID`. I think it is useful when the enrichment analysis returns no significant results.

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
<img src="https://github.com/noriakis/wcGeneSummary/blob/main/images/erccWc.png?raw=true" width="800px">

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
<img src="https://github.com/noriakis/wcGeneSummary/blob/main/images/cxclWc.png?raw=true" width="800px">


### The other examples

- [Example annotating gene clusters using WGCNA](https://noriakis.github.io/software/wcGeneSummary/)
- [Example image of CCL genes (correlation network)](https://github.com/noriakis/software/blob/main/images/cclNet.png)
- [Example in Bayesian network analysis](https://github.com/noriakis/compare_sign)

### References
[Lucas T. palettetown: Pokemon themed colour schemes for R.](https://github.com/timcdlucas/palettetown)