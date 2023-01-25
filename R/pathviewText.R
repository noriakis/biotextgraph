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
#' @param orgDb organism database to convert symbol
#' @param target target DB for word, "refseq" or "abstract"
#' @param numWords number of words in barplot
#' @param node.types passed to node.map e.g. "genes" and "ortholog"
#' @param searchTerms search terms to be used in wcAbst,
#' if different than query geneList
#' @param termMap data.frame consisting of columns
#' "query" and "description", which shows query as "geneList",
#' and description as "searchTerms".
#' e.g. "K08097" and "phosphosulfolactate synthase"
#' @param areas used in patchwork
#' @param trans transpose the barplot
#' @param argList passed to wc functions
#' @import grid
#' @export
#'
#' 
pathviewText <- function(geneList, keyType, pid, org="hsa",
                         pal="RdBu",target="refseq",
                         searchTerms=NULL, node.types="gene",
                         termMap=NULL, orgDb=org.Hs.eg.db,
                         numWords=20, trans=FALSE, areas=NULL, argList=list()) {
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
    if (length(allCol)!=length(colorNum)) {
        allCol <- colorRampPalette(allCol)(colorNum)
    }
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
        argList[["queries"]] <- rawGenes
        argList[["keyType"]] <- "ENTREZID"
        argList[["numWords"]] <- numWords
        argList[["plotType"]] <- "network"
        argList[["genePlot"]] <- TRUE
        argList[["preserve"]] <- TRUE
        argList[["genePlotNum"]] <- length(geneList)
        barp <- do.call("wcAbst", argList)
        barp@geneMap[,2] <- gsub(" \\(Q\\)", "", barp@geneMap[,2])
    } else {
        argList[["geneList"]] <- geneList
        argList[["keyType"]] <- "ENTREZID"
        argList[["numWords"]] <- numWords
        argList[["plotType"]] <- "network"
        argList[["genePlot"]] <- TRUE
        argList[["preserve"]] <- TRUE
        argList[["genePlotNum"]] <- length(geneList)
        barp <- do.call("wcGeneSummary", argList)
    }
    
    ## Obtain words related to genes listed in map
    gmap <- barp@geneMap
    candWords <- NULL
    if (is.null(termMap)) {
        for (gn in plot.data.gene[,2]) {
            if (gn %in% barp@geneMap[, 2]) {
                candWords <- rbind(candWords, gmap[gmap[,2] %in% gn,])
            }
        }
        candWords <- data.frame(candWords) |> `colnames<-`(c("word","query"))
    } else {
        for (gn in plot.data.gene[,2]) {
            if (gn %in% termMap$query) {
                mappedDesc <- unique(subset(termMap, query==gn)$description)
                for (dc in mappedDesc) {
                    if (dc %in% gmap[, 2]) {
                        tmpgmap <- data.frame(gmap[gmap[,2] %in% dc,]) |>
                            `colnames<-`(c("word","description"))
                        tmpgmap$query <- rep(gn, dim(tmpgmap)[1])
                        candWords <- rbind(candWords, tmpgmap)
                    }
                }
            }
        }    
    }
    
    
    inBar <- candWords[candWords$word %in% barp@freqDf$word,]
    rePlot <- barp@freqDf[1:numWords,]
    barCols <- NULL
    rePlotColor <- NULL
    for (rn in row.names(rePlot)) {
        wn <- rePlot[rn,]$word
        if (rn %in% inBar$word) {
            takeGene <- inBar[inBar$word %in% rn,]
            if (is.null(dim(takeGene))) {
                gn <- takeGene$query
            } else {
                gn <- unique(takeGene$query)
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
    
    changeX1 <- NULL
    for (cn in tolower(rePlotColor$X1)) {
        if (cn %in% names(barp@dic)) {
            changeX1 <- c(changeX1, barp@dic[cn])
        } else {
            changeX1 <- c(changeX1, cn)
        }
    }
    rePlotColor$X1 <- changeX1


    if (trans) {
    
        replot <- rePlotColor |> ggplot(aes(x=reorder(X1,X3),y=stack,fill=X2))+
            geom_bar(position="stack", stat="identity")+
            scale_fill_manual(values=colVec, name="Gene")+
            xlab("Words")+ylab("Frequency")+theme_minimal()+
            theme(axis.text =element_text(angle=90))
        
        if (is.null(areas)) {
            areas <- "
          AAAAAA
          AAAAAA
          AAAAAA
          AAAAAA
          #BBBB#
        "
        }
        plt <- patchwork::wrap_plots(g, replot, ncol=1)+
            plot_layout(design=areas)
    } else {
        replot <- rePlotColor |> ggplot(aes(y=reorder(X1,X3),x=stack,fill=X2))+
            geom_bar(position="stack", stat="identity")+
            scale_fill_manual(values=colVec, name="Gene")+
            ylab("Words")+xlab("Frequency")+theme_minimal()
        if (is.null(areas)) {
            areas <- "
            AAAAAABB
            AAAAAABB
            AAAAAABB
            AAAAAABB
            "
        }
        plt <- patchwork::wrap_plots(g, replot, nrow=1)+
            plot_layout(design=areas)
    }
    
    returnList[["text"]] <- barp
    returnList[["bar"]] <- replot
    returnList[["pngGrob"]] <- g
    returnList[["concat"]] <- plt
    return(returnList)
}