#'
#' pathviewText
#' 
#' colorize pathview map, and link to word frequency
#' 
#' @param geneList gene ID list
#' @param keyType keytype for queried gene
#' @param pid pathway id
#' @param pal palette to be passed to RColorBrewer
#' for coloring of nodes
#' @param org organism ID in KEGG, default to hsa
#' @param target target DB for word
#' @param numWords number of words in barplot
#' @param node.types passed to node.map
#' @param ... passed to wc functions
#' @import png grid pathview
#' @export
#'
#' 
pathviewText <- function(geneList, keyType, pid, org="hsa",
  pal="RdBu",target="refseq", searchTerms=NULL, node.types="gene",
  numWords=20, ...) {
  returnList <- list()
  if (!keyType %in% c("KO","ENTREZID")) {
    qqcat("converting to ENTREZID\n")
    rawGenes <- geneList
    geneList <- AnnotationDbi::select(orgDb,
                                      keys = geneList, columns = c("ENTREZID"),
                                      keytype = keyType)$ENTREZID
    geneList <- geneList[!is.na(geneList)]
    qqcat("converted input genes: @{length(geneList)}\n")
  } else {
    rawGenes <- geneList
  }



  vec <- rep(1, length(geneList))
  names(vec) <- geneList
  
  download.kegg(pathway.id=pid,
    species=org)
  xml.file=paste0(org, pid,".xml")
  node.data=node.info(xml.file)
  plot.data.gene=node.map(mol.data=vec,
                          node.data,
                          node.types=node.types)



  # pvParams <- pathview(geneList,pathway.id=pid,limit=list(gene=1))
  # exportedPV <- pvParams$plot.data.gene
  colorNum <- length(plot.data.gene$mol.data[!is.na(plot.data.gene$mol.data)])
  allCol <- RColorBrewer::brewer.pal(colorNum, pal)
  if (length(allCol) > colorNum) {
    allCol <- allCol[1:colorNum]
  }


  molCol <- NULL
  nm <- 1
  colGeneMap <- NULL
  for (i in row.names(plot.data.gene)){
    if (!is.na(plot.data.gene[i,]$mol.data)) {
      if (plot.data.gene[i,]$labels %in% colGeneMap[,2]) {
        tmpCol <- unique(colGeneMap[colGeneMap[,2] %in% plot.data.gene[i,]$labels, 1])
        molCol <- c(molCol, tmpCol)
        colGeneMap <- rbind(colGeneMap, c(tmpCol, plot.data.gene[i,]$labels))
      } else {
        molCol <- c(molCol, allCol[nm])
        colGeneMap <- rbind(colGeneMap, c(allCol[nm], plot.data.gene[i,]$labels))
      }
      nm <- nm + 1
    } else {
      molCol <- c(molCol, "#FFFFFF")
    }
  }

  pv.pars= keggview.native(plot.data.gene=plot.data.gene,
                           cols.ts.gene=molCol,
                           pathway.name=paste0(org,pid),
                           same.layer=TRUE,
                           plot.col.key=FALSE,
                           out.suffix = "custom.cols")


  img <- readPNG(paste0(org,pid,".custom.cols.png"))
  g <- rasterGrob(img, interpolate=FALSE)
  if (target=="abstract") {
    if (!is.null(searchTerms)) {
      rawGenes <- searchTerms
    }
    barp <- wcAbst(rawGenes, keyType="ENTREZID",
                   numWords=numWords,
                    plotType="network",
                    genePlot=TRUE,
                    genePlotNum=length(geneList), ...)
  } else {
    barp <- wcGeneSummary(geneList, keyType="ENTREZID",
                          numWords=numWords,
                          plotType="network",
                          genePlot=TRUE,
                          genePlotNum=length(geneList), ...)
  }

  ## Obtain words related to genes listed in map
  gmap <- barp@geneMap
  candWords <- NULL
  for (gn in plot.data.gene[,2]) {
    if (gn %in% barp@geneMap[, 2]) {
      candWords <- rbind(candWords, gmap[gmap[,2] %in% gn,])
    }
  }

  inBar <- candWords[candWords[,1] %in% barp@freqDf$word,]
  rePlot <- barp@freqDf[1:numWords,]
  barCols <- NULL
  rePlotColor <- NULL
  for (rn in row.names(rePlot)) {
    wn <- rePlot[rn,]$word
    if (rn %in% inBar[,1]) {
      takeGene <- inBar[inBar[,1] %in% rn,]
      if (is.null(dim(takeGene))) {
        gn <- takeGene[2]
      } else {
        gn <- unique(takeGene[,2])
      }
      for (eg in gn) {
        freq <- rePlot[rn,]$freq
        col <- unique(subset(colGeneMap, colGeneMap[,2] %in% eg)[,1])
        rePlotColor <- rbind(rePlotColor, c(wn, eg, freq, col))
      }
    } else {
      gn <- NA
      freq <- rePlot[rn,]$freq
      col <- "grey"    
      rePlotColor <- rbind(rePlotColor, c(wn, gn, freq, col))
    }
  }


  rePlotColor <- data.frame(rePlotColor)
  rePlotColor$X3 <- as.numeric(rePlotColor$X3)
  occuNum <- rePlotColor |> group_by(X1) |> summarise(n=n())
  occuNumVec <- occuNum$n
  names(occuNumVec) <- occuNum$X1
  
  newX3 <- NULL
  for (i in row.names(rePlotColor)) {
    newNum <- as.numeric(rePlotColor[i,]$X3 / occuNumVec[rePlotColor[i,]$X1])
    newX3 <- c(newX3, newNum)
  }
  rePlotColor$stack <- newX3 
  colMat <- colGeneMap[!duplicated(colGeneMap[,2]),]
  if (is.null(dim(colMat))) {
    colVec <- colMat[1]
    names(colVec) <- colMat[2]
  } else {
    colVec <- colMat[,1]
    names(colVec) <- colMat[,2]
  }
  

  replot <- rePlotColor |> ggplot(aes(x=reorder(X1,X3),y=stack,fill=X2))+
    geom_bar(position="stack", stat="identity")+
    scale_fill_manual(values=colVec, name="Gene")+
    xlab("Words")+ylab("Frequency")+theme_minimal()+
    theme(axis.text =element_text(angle=90))


  areas <- "
    AAAAAA
    AAAAAA
    AAAAAA
    AAAAAA
    #BBBB#
  "
  plt <- patchwork::wrap_plots(g, replot, ncol=1)+
    plot_layout(design=areas)
  returnList[["bar"]] <- replot
  returnList[["pngGrob"]] <- g
  returnList[["concat"]] <- plt
  return(returnList)
}

