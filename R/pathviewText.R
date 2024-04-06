# #'
# #' pathviewText
# #' 
# #' colorize pathview map, and link to word frequency
# #' 
# #' @param geneList gene ID list
# #' @param keyType keytype for queried gene
# #' @param pid pathway id
# #' @param pal palette to be passed to RColorBrewer
# #' for coloring of nodes
# #' @param org organism ID in KEGG, default to hsa
# #' @param orgDb organism database to convert symbol
# #' @param target target DB for word, "refseq" or "abstract"
# #' @param numWords number of words in barplot
# #' @param node.types passed to node.map e.g. "genes" and "ortholog"
# #' @param searchTerms search terms to be used in pubmed,
# #' if different than query geneList
# #' @param termMap data.frame consisting of columns
# #' "query" and "description", which shows query as "geneList",
# #' and description as "searchTerms".
# #' e.g. "K08097" and "phosphosulfolactate synthase"
# #' @param areas used in patchwork
# #' @param trans transpose the barplot
# #' @param argList passed to wc functions
# #' @param textSize text size in barplot
# #' @return biotext object, list of plots of pathway, and barplot and concatentaed image
# #' @import grid
# #' @importFrom utils head
# #' @importFrom RColorBrewer brewer.pal
# #' @examples
# #' query <- c("TP53","CDC45","CDC6")
# #' \dontrun{pathviewText(query, keyType = "SYMBOL", pid = "04110", org = "hsa")}
# #' @export
# #'
# #' 
# pathviewText <- function(geneList, keyType, pid, org="hsa",
#                          pal="RdBu",target="refseq",
#                          searchTerms=NULL, node.types="gene",
#                          termMap=NULL, orgDb=org.Hs.eg.db, textSize=12,
#                          numWords=20, trans=FALSE, areas=NULL, argList=list()) {
#     returnList <- list()
#     if (!keyType %in% c("KO","ENTREZID")) {
#         qqcat("Converting to ENTREZID\n")
#         rawGenes <- geneList
#         geneList <- AnnotationDbi::select(orgDb,
#                                           keys = geneList, columns = c("ENTREZID"),
#                                           keytype = keyType)$ENTREZID
#         geneList <- geneList[!is.na(geneList)]
#         qqcat("  Converted input genes: @{length(geneList)}\n")
#     } else {
#         rawGenes <- geneList
#     }
    
    
    
#     vec <- rep(1, length(geneList))
#     names(vec) <- geneList
    
#     pathview::download.kegg(pathway.id=pid,
#                   species=org)
#     xml.file <- paste0(org, pid,".xml")
#     node.data <- pathview::node.info(xml.file)
#     plot.data.gene <- pathview::node.map(mol.data=vec,
#                             node.data,
#                             node.types=node.types)
    
    
    
#     # pvParams <- pathview(geneList,pathway.id=pid,limit=list(gene=1))
#     # exportedPV <- pvParams$plot.data.gene
#     colorNum <- length(plot.data.gene$mol.data[!is.na(plot.data.gene$mol.data)])
#     allCol <- RColorBrewer::brewer.pal(colorNum, pal)
#     if (length(allCol)!=length(colorNum)) {
#         allCol <- colorRampPalette(allCol)(colorNum)
#     }
#     if (length(allCol) > colorNum) {
#         allCol <- allCol[1:colorNum]
#     }
    
    
#     molCol <- NULL
#     nm <- 1
#     colGeneMap <- NULL
#     for (i in row.names(plot.data.gene)){
#         if (!is.na(plot.data.gene[i,]$mol.data)) {
#             if (plot.data.gene[i,]$labels %in% colGeneMap[,2]) {
#                 tmpCol <- unique(colGeneMap[colGeneMap[,2] %in% plot.data.gene[i,]$labels, 1])
#                 molCol <- c(molCol, tmpCol)
#                 colGeneMap <- rbind(colGeneMap, c(tmpCol, plot.data.gene[i,]$labels))
#             } else {
#                 molCol <- c(molCol, allCol[nm])
#                 colGeneMap <- rbind(colGeneMap, c(allCol[nm], plot.data.gene[i,]$labels))
#             }
#             nm <- nm + 1
#         } else {
#             molCol <- c(molCol, "#FFFFFF")
#         }
#     }

#     pv.pars <- pathview::keggview.native(plot.data.gene=plot.data.gene,
#                              cols.ts.gene=molCol,
#                              pathway.name=paste0(org,pid),
#                              same.layer=TRUE,
#                              plot.col.key=FALSE,
#                              out.suffix = "custom.cols")
    
    
#     img <- png::readPNG(paste0(org,pid,".custom.cols.png"))
#     g <- rasterGrob(img, interpolate=FALSE)
#     if (target=="abstract") {
#         if (!is.null(searchTerms)) {
#             rawGenes <- searchTerms
#         }
#         argList[["queries"]] <- rawGenes
#         # argList[["numWords"]] <- numWords
#         argList[["plotType"]] <- "network"
#         argList[["genePlot"]] <- TRUE
#         # argList[["preserve"]] <- TRUE
#         # argList[["genePlotNum"]] <- length(geneList)
#         barp <- do.call("pubmed", argList)
#         barp@geneMap[,2] <- gsub(" \\(Q\\)", "", barp@geneMap[,2])
#     } else {
#         argList[["geneList"]] <- geneList
#         argList[["keyType"]] <- "ENTREZID"
#         # argList[["numWords"]] <- numWords
#         argList[["plotType"]] <- "network"
#         argList[["genePlot"]] <- TRUE
#         # argList[["preserve"]] <- TRUE
#         argList[["genePlotNum"]] <- length(geneList)
#         barp <- do.call("refseq", argList)
#     }
    
#     ## Obtain words related to genes listed in map
#     gmap <- barp@geneMap ## This not set to all lower
#     candWords <- NULL
#     if (is.null(termMap)) {
#         for (gn in unique(plot.data.gene[,2])) {
#             if (tolower(gn) %in% tolower(gmap[, 2])) {
#                 candWords <- rbind(candWords, gmap[tolower(gmap[,2]) %in% tolower(gn),])
#             }
#         }
#         candWords <- data.frame(candWords) |> `colnames<-`(c("word","query"))
#     } else {
#         for (gn in unique(plot.data.gene[,2])) {
#             if (gn %in% termMap$query) {
#                 mappedDesc <- unique(subset(termMap, termMap$query==gn)$description)
#                 for (dc in mappedDesc) {
#                     if (tolower(dc) %in% tolower(gmap[, 2])) {
#                         tmpgmap <- data.frame(gmap[tolower(gmap[,2]) %in% tolower(dc),]) |>
#                             `colnames<-`(c("word","description"))
#                         tmpgmap$query <- rep(gn, dim(tmpgmap)[1])
#                         candWords <- rbind(candWords, tmpgmap)
#                     }
#                 }
#             }
#         }    
#     }

#     frDf <- barp@freqDf
#     inBar <- candWords[tolower(candWords$word) %in% tolower(frDf$word),]
#     rePlot <- frDf#[1:numWords,]
#     barCols <- NULL
#     rePlotColor <- NULL

#     for (rn in row.names(rePlot)) {
#         wn <- rePlot[rn,]$word
#         if (tolower(wn) %in% tolower(inBar$word)) {
#             takeGene <- inBar[tolower(inBar$word) %in% tolower(wn),]
#             # if (is.null(dim(takeGene))) {
#                 # gn <- takeGene$query
#             # } else {
#             gn <- unique(takeGene$query)
#             # }
#             for (eg in gn) {
#                 freq <- rePlot[rn,]$freq
#                 col <- unique(subset(colGeneMap, colGeneMap[,2] %in% eg)[,1])
#                 rePlotColor <- rbind(rePlotColor, c(wn, eg, freq, col))
#             }
#         } else {
#             gn <- NA
#             freq <- rePlot[rn,]$freq
#             col <- "grey"    
#             rePlotColor <- rbind(rePlotColor, c(wn, gn, freq, col))
#         }
#     }
    
    
#     rePlotColor <- data.frame(rePlotColor) |> `colnames<-`(c("word","query","freq","color"))
#     rePlotColor$freq <- as.numeric(rePlotColor$freq)
#     occuNum <- rePlotColor |> dplyr::group_by(rePlotColor$word) |> 
#     dplyr::summarise(n=dplyr::n()) |> `colnames<-`(c("word","n"))
#     occuNumVec <- occuNum$n
#     names(occuNumVec) <- occuNum$word

#     relFreqs <- NULL
#     for (i in row.names(rePlotColor)) {
#         relFreq <- as.numeric(rePlotColor[i,]$freq / occuNumVec[rePlotColor[i,]$word])
#         relFreqs <- c(relFreqs, relFreq)
#     }
#     rePlotColor$stack <- relFreqs

#     colMat <- colGeneMap[!duplicated(colGeneMap[,2]),]
#     if (is.null(dim(colMat))) {
#         colVec <- colMat[1]
#         names(colVec) <- colMat[2]
#     } else {
#         colVec <- colMat[,1]
#         names(colVec) <- colMat[,2]
#     }
    
#     changedWord <- NULL
#     for (cn in tolower(rePlotColor$word)) {
#         if (cn %in% names(barp@dic)) {
#             changedWord <- c(changedWord, barp@dic[cn])
#         } else {
#             changedWord <- c(changedWord, cn)
#         }
#     }
#     rePlotColor$word <- changedWord

#     if (trans) {
    
#         replot <- rePlotColor |> head(n=numWords) |> 
#         ggplot(aes(x=reorder(.data$word,freq),
#             y=.data$stack,fill=.data$query))+
#             geom_bar(position="stack", stat="identity")+
#             scale_fill_manual(values=colVec, name="Gene")+
#             xlab("Words")+ylab("Frequency")+theme_minimal()+
#             theme(axis.text =element_text(angle=90, size=textSize))
        
#         if (is.null(areas)) {
#             areas <- "
#           AAAAAA
#           AAAAAA
#           AAAAAA
#           AAAAAA
#           #BBBB#
#         "
#         }
#         plt <- patchwork::wrap_plots(g, replot, ncol=1)+
#             plot_layout(design=areas)
#     } else {
#         replot <- rePlotColor |> head(n=numWords) |> 
#         ggplot(aes(y=reorder(.data$word,freq),
#             x=.data$stack,fill=.data$query))+
#             geom_bar(position="stack", stat="identity")+
#             scale_fill_manual(values=colVec, name="Gene")+
#             ylab("Words")+xlab("Frequency")+theme_minimal()+
#             theme(axis.text =element_text(size=textSize))
#         if (is.null(areas)) {
#             areas <- "
#             AAAAAABB
#             AAAAAABB
#             AAAAAABB
#             AAAAAABB
#             "
#         }
#         plt <- patchwork::wrap_plots(g, replot, nrow=1)+
#             plot_layout(design=areas)
#     }
    
#     returnList[["text"]] <- barp
#     returnList[["bar"]] <- replot
#     returnList[["pngGrob"]] <- g
#     returnList[["concat"]] <- plt
#     return(returnList)
# }