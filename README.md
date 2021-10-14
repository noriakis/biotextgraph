# wcGeneSummary
Make word cloud of gene set using R libraries `GeneSummary`, `tm` and `wordcloud`. Input is gene list with the type `ENTREZID`. I think it is useful when GSEA or ORA returned no results.

```R
ercc <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")
entrezID = AnnotationDbi::select(org.Hs.eg.db, keys = ercc, columns = c("ENTREZID"), keytype = "SYMBOL")$ENTREZID
gwc <- wcGeneSummary(entrezID, excludeFreq = 5000, max.words=200, random.order=FALSE,
                     colors=palettetown::pokepal(150), shape="circle", rot.per=0.4)
# ggsave("erccWc.png", gwc$wc, width=8, height=8)
```

<img src="https://github.com/noriakis/wcGeneSummary/blob/main/images/erccWc.png?raw=true" width="800px">