#' @rdname generalf
#' @import tm
#' @import bugsigdbr
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @importFrom XML xmlParse xpathSApply xmlTreeParse xmlElementsByTagName xmlValue
#' @importFrom dplyr filter
#' @importFrom stats dist
#' @importFrom grDevices palette
#' @importFrom tidygraph tbl_graph
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency graph_from_edgelist
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' 
#' @examples
#' mbList <- c("Veillonella dispar")
#' \dontrun{
#'     bugsigdb(mbList, plotType="wc")
#' }
#' @export
#' 
bugsigdb <- function (mbList,
    excludeFreq=1000, exclude="frequency",
    excludeType=">", normalize=FALSE, takeMean=FALSE,
    additionalRemove=NA, tfidf=FALSE,
    target="title", apiKey=NULL, takeMax=FALSE,
    pre=FALSE, pvclAlpha=0.95, numOnly=TRUE,
    madeUpper=c("dna","rna"), redo=NULL, fontFamily="sans",
    pal=c("blue","red"), numWords=15, preserve=FALSE,
    metab=NULL, metThresh=0.2, curate=TRUE,
    abstArg=list(), tagPalette=NULL, metCol=NULL,
    scaleRange=c(5,10), showLegend=FALSE, ecPlot=FALSE,
    edgeLabel=FALSE, mbPlot=FALSE, onlyTDM=FALSE,
    ecFile=NULL, upTaxFile=NULL, filterMax=FALSE, mbColor="grey",
    useUdpipe=FALSE, colorize=FALSE, cooccurrence=FALSE,
    udpipeModel="english-ewt-ud-2.5-191206.udpipe", scaleFreq=NULL,
    ngram=1, plotType="network", disPlot=FALSE, onWholeDTM=FALSE,
    naEdgeColor="grey50", useggwordcloud=TRUE, wcScale=10,addFreqToMB=FALSE,
    catColors=NULL, useSeed=42,discreteColorWord=FALSE,
    colorText=FALSE, corThresh=0.2, tag="none", tagWhole=FALSE, stem=FALSE,
    layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, cl=FALSE,
    autoThresh=TRUE, argList=list(), docsum=FALSE,  absolute=TRUE,
    corOption=list(), cache=TRUE) {
    

    if (useUdpipe) {
        qqcat("Using udpipe mode\n")
        plotType="network"
        ngram <- 1

        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)
    }
    ret <- new("biotext")
    ret@query <- mbList
    ret@type <- paste0("BSDB_",target)
    addNet <- list()

    if (pre) {additionalRemove <- c(additionalRemove, "microbiota",
        "microbiome","relative","abundance","abundances",
        "including","samples","sample","otu","otus","investigated",
        "taxa","taxon")}
    
    if (is.null(redo)) {
        qqcat("Input microbes: @{length(mbList)}\n")
        tb <- importBugSigDB(cache=cache)
        subTb <- c()        
        for (m in mbList) {
            tmp <- tb[grepl(m, tb$`MetaPhlAn taxon names`,
                     fixed = TRUE),]
            if (dim(tmp)[1]>0) {
                qqcat("  Found @{length(unique(tmp$PMID))} entries for @{m}\n")
                tmp$query <- m
                subTb <- rbind(subTb, tmp)
            }
        }
        qqcat("Including @{dim(subTb)[1]} entries\n")

        titles <- unique(subTb$Title)
        titles <- titles[!is.na(titles)]
        fil <- c()
        for (title in titles){
            tmp <- subset(subTb, subTb$Title==title)
            if (dim(tmp)[1]>1){
                ## Duplicated entry is deleted based on query and paper title
                tmp$title2 <- paste0(tmp$Title,"_",tmp$query)
                tmp2 <- tmp[!duplicated(tmp$title2),]
                fil <- rbind(fil, c(paste(tmp2$query, collapse=","),
                    unique(tmp2$Title), unique(tmp2$PMID)))
            } else {
                fil <- rbind(fil, c(tmp$query, tmp$Title, tmp$PMID))
            }
        }
        fil <- data.frame(fil)
        colnames(fil) <- c("query","text","ID")
        fil <- fil[!is.na(fil$ID),] # Some PMIDs have NA
        # returnList[["filterWords"]] <- filterWords
        ret@pmids <- fil$ID
        ret@rawTextBSDB <- subTb

    } else {
        qqcat("Redoing abstract query for microbes\n")
        ret <- redo
        target <- "abstract"
    }
    if (curate) {
        if (target=="abstract"){
            qqcat("Target is abstract\n")
            if (is.null(redo)) {
                pmids <- unique(fil$ID)
                # pmids <- pmids[!is.na(pmids)]
                qqcat("  Querying PubMed for @{length(pmids)} pmids\n")
                queryUrl <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                                   "efetch.fcgi?db=pubmed&retmode=XML&id=",
                                       paste(pmids, collapse=","))
                if (!is.null(apiKey)){
                    queryUrl <- paste0(queryUrl, "&api_key=", apiKey)
                } else {
                    qqcat("  Querying without API key\n")
                }
                onequery <- url(queryUrl, open = "rb", encoding = "UTF8")
                xmlString <- readLines(onequery, encoding="UTF8")
                parsedXML <- xmlParse(xmlString)
                obtainText <- function(pmid) {xpathSApply(parsedXML,
                    paste('//PubmedArticle/MedlineCitation[PMID=',pmid,
                        ']/Article/Abstract/AbstractText'),
                    xmlValue)}
                PMIDs <- as.numeric(xpathSApply(parsedXML,
                    "//PubmedArticle/MedlineCitation/PMID",
                    xmlValue))
                abstDf <- NULL
                for (pmid in PMIDs) {
                    if (length(obtainText(pmid))!=0) {
                        for (text in obtainText(pmid)) {
                            tax <- unique(subset(fil, fil$ID==pmid)$query)
                            abstDf <- rbind(abstDf, c(pmid, text, tax))
                        }
                    }
                }
                abstDf <- data.frame(abstDf) |> `colnames<-`(c("ID","text","query"))
                # abstset <- xmlElementsByTagName(parsedXML$doc$children$PubmedArticleSet,
                #                                 "Abstract",
                #                                 recursive = TRUE)
                # absttext <- as.character(xmlValue(abstset))
                ret@rawText <- abstDf
            } else {
                abstDf <- ret@rawText
                filterWords <- ret@filtered
                subTb <- ret@rawTextBSDB
            }
            ## TODO: ask taking unique text or not?
            docs <- VCorpus(VectorSource(abstDf$text))
        } else {
            docs <- VCorpus(VectorSource(fil$text))
        }
    } else {
        abstArg[["queries"]] <- mbList
        abstArg[["quote"]] <- TRUE
        abstArg[["target"]] <- target
        abstArg[["onlyDf"]] <- TRUE

        abstDf <- do.call(pubmed, abstArg)
        ret@rawText <- abstDf
        ret@type <- paste0("BSDB_PubMed_",target)
        docs <- VCorpus(VectorSource(abstDf$text))
    }
    ## Make corpus
    ## Filter high frequency words if needed
    if (exclude=="frequency") {
        pref = "BSDB_Freq"
    } else {
        pref = "BSDB_TfIdf"
    }
    if (filterMax) {
        pref <- paste0(pref, "_Max")
    }
    filterWords <- retFiltWords(useFil=pref, filType=excludeType, filNum=excludeFreq)

    qqcat("Filtering @{length(filterWords)} words (frequency and/or tfidf)\n")
    if (preserve) {
        pdic <- preserveDict(docs, ngram, numOnly, stem)
        ret@dic <- pdic
    }
    docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly, stem)
    if (length(filterWords)!=0 | length(additionalRemove)!=0){
        allfils <- c(filterWords, additionalRemove)
        allfils <- allfils[!is.na(allfils)]
        if (length(allfils)!=0) {
            ret@filtered <- allfils
        }
    }
    ret@corpus <- docs
    ret@tag <- tag

    ## Make additional network
    if (disPlot & plotType=="network") {## This does not need to be deduplicated
        dis <- NULL
        for (i in seq_len(nrow(subTb))){
            dis <- rbind(dis,
            c(subTb[i, "Condition"], subTb[i, "query"]))
        }
        dis <- dis[!is.na(dis[,1]),]
        dis <- dis[!is.na(dis[,2]),]
        vtx <- data.frame(cbind(c(dis[,1], dis[,2]), c(rep("Diseases",length(dis[,1])),
            rep("Microbes",length(dis[,2]))))) |> `colnames<-`(c("name","type"))
        vtx <- vtx[!duplicated(vtx),]
        dmg <- tbl_graph(nodes=vtx, edges=data.frame(dis), directed=FALSE)
        addNet[["Diseases"]] <- dmg
    }

    ## Add metab if present
    ## TODO: somehow show edge weights other than
    ## correlation between words
    if (!is.null(metab) & plotType=="network") {
        if (is.null(metCol)) {
            stop("No column names specified")
        }
        qqcat("Checking metabolites\n")
        metabGraph <- NULL
        for (sp in mbList) {
            tmp <- metab[grepl(sp,metab[[metCol[1]]]),]
            tmp <- tmp[abs(tmp[[metCol[3]]])>metThresh,]
            if (dim(tmp)[1]!=0) {
                for (met in tmp[[metCol[2]]]) {
                    metabGraph <- rbind(metabGraph, c(sp, met))
                }
            } else {
                qqcat("  Found no metabolites for @{sp}\n")
            }
        }
        if (!is.null(metabGraph)) {
            vtx <- data.frame(
                cbind(c(metabGraph[,1], metabGraph[,2]),
                      c(rep("Microbes",length(metabGraph[,1])),
                        rep("Metabolites",length(metabGraph[,2]))))) |> `colnames<-`(c("name","type"))
            vtx <- vtx[!duplicated(vtx),]
            metabGraph <- tbl_graph(nodes=vtx, edges=data.frame(metabGraph), directed=FALSE)
            addNet[["Metabolites"]] <- metabGraph
        }
    }

    ## Add EC if present
    if (ecPlot & plotType=="network") {
        mbPlot <- TRUE
        if (is.null(ecFile)) {stop("Please provide EC file")}
        if (is.null(upTaxFile)) {stop("Please provide UniProt taxonomy file")}
        ecDf <- enzyme(file=ecFile, ecnum="all", taxec=TRUE,
            taxFile=upTaxFile, candTax=mbList)
        if (!is.null(ecDf)) {
            ecDf <- ecDf[,c("desc","query")]
            vtx <- data.frame(
                cbind(c(ecDf[,1], ecDf[,2]),
                      c(rep("Enzymes",length(ecDf[,1])),
                        rep("Microbes",length(ecDf[,2]))))) |> 
            `colnames<-`(c("name","type"))

            vtx <- vtx[!duplicated(vtx),]
            ecGraph <- tbl_graph(nodes=vtx,
                edges=data.frame(ecDf),
                directed=FALSE)
            addNet[["Enzymes"]] <- ecGraph
           # ecg <- simplify(graph_from_data_frame(, 
           #      directed=FALSE))
           # addNet[["Enzymes"]] <- as_tbl_graph(ecg)
        }
    }


    if (useUdpipe) {
        if (curate & target!="abstract") {
            abstDf <- fil
        }
        ret <- retUdpipeNet(ret=ret,texts=abstDf,udmodel_english=udmodel_english,
          orgDb=NULL, filterWords=filterWords, additionalRemove=additionalRemove,
          colorText=colorText,edgeLink=edgeLink,queryPlot=mbPlot, layout=layout,
          pal=pal,showNeighbors=NULL, showFreq=NULL, nodePal=tagPalette,addNet=addNet,
          queryName="Microbes")
        return(ret)
    }


    ## Set parameters for correlation network
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    if (ngram!=1){
        NgramTokenizer <- function(x)
            unlist(lapply(ngrams(words(x), ngram),
                paste, collapse = " "),
                use.names = FALSE)
        if (tfidf) {
            docs <- TermDocumentMatrix(docs,
                control = list(tokenize = NgramTokenizer,
                    weighting = weightTfIdf))
        } else {
            docs <- TermDocumentMatrix(docs,
                control = list(tokenize = NgramTokenizer))
        }
    } else {
        if (tfidf) {
            docs <- TermDocumentMatrix(docs,
                control = list(weighting = weightTfIdf))
        } else {
            docs <- TermDocumentMatrix(docs)
        }
    }

    mat <- as.matrix(docs)
    if (docsum) {
        mat <- apply(mat, 2, function(x) ifelse(x>0, 1, 0))
    }    
    if (normalize) {
        mat <- sweep(mat, 2, colSums(mat), `/`)
    }

    if (takeMax & takeMean) {stop("Should either of specify takeMax or takeMean")}
    if (takeMax) {
        perterm <- apply(mat, 1, max, na.rm=TRUE)
    } else {
        if (takeMean) {
            perterm <- apply(mat,1,mean)
        } else {
            perterm <- rowSums(mat)
        }
    }

    matSorted <- sort(perterm, decreasing=TRUE)
    ret@wholeFreq <- matSorted

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }
    ret@numWords <- numWords
    if (onlyTDM) {
        return(docs)
    }
    # returnList[["rawfrequency"]] <- matSorted
    ret@TDM <- docs

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        freqWords <- names(matSorted)
        DTM <- t(as.matrix(docs))
        returnDf <- data.frame(word = names(matSorted),freq=matSorted) |>
            na.omit()
        ret@freqDf <- returnDf

        if (tag=="tdm") {# Needs rework
            if (!is.null(redo)){
                if (!is.null(redo@pvpick)) {
                    pvcl <- ret@pvpick
                } else {
                    if (tagWhole){
                        pvc <- pvclust(as.matrix(dist(t(DTM))),parallel=cl, method.dist=tag)
                    } else {
                        pvc <- pvclust(as.matrix(dist(
                            t(
                                DTM[,colnames(DTM) %in% freqWords]
                                )
                            )), parallel=cl, method.dist=tag)
                    }
                    pvcl <- pvpick(pvc, alpha=pvclAlpha)
                    ret@pvclust <- pvc
                    ret@pvpick <- pvcl
                }
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(t(DTM))),parallel=cl, method.dist=tag)
                } else {
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            DTM[,colnames(DTM) %in% freqWords]
                            )
                        )),parallel=cl, method.dist=tag)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
        }

        ## mbPlot: plot associated mbs
        if (disPlot) {mbPlot <- TRUE}
        if (!is.null(metab)) {mbPlot <- TRUE}
        if (ecPlot) {mbPlot <- TRUE}

        if (mbPlot) {
            if (target=="abstract") {
                row.names(DTM) <- abstDf$query
            } else {
                if (curate) {
                    row.names(DTM) <- fil$query
                } else {
                    row.names(DTM) <- abstDf$query
                }
            }
        }

        matrixs <- obtainMatrix(ret, FALSE, NULL, DTM, freqWords,
            corThresh, cooccurrence, onWholeDTM, numWords, autoThresh, absolute, corOption)
        
        coGraph <- matrixs$coGraph
        ret <- matrixs$ret
        if (tag=="cor") {
			ret <- tag_words(ret, cl,
				pvclAlpha, whole=tagWhole,
				num_words=ret@numWords)
            pvc <- ret@pvclust
            pvcl <- ret@pvpick
        }
        ret@igraphRaw <- coGraph

        ## Subset to frequent-words
        coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }

        nodeName <- V(coGraph)$name
        # dtmCol <- colnames(DTM)
        # for (i in madeUpper) {
        #     dtmCol[dtmCol == i] <- toupper(i)
        #     nodeName[nodeName == i] <- toupper(i)
        # }
        # V(coGraph)$name <- nodeName
        # colnames(DTM) <- dtmCol

        if (mbPlot) {
            mbmap <- c()
            for (rn in nodeName){
                tmp <- DTM[ ,rn]
                for (nm in names(tmp[tmp!=0])){
                    if (grepl(",", nm, fixed = TRUE)) {
                        for (microbe in unlist(strsplit(nm, ","))){
                            mbmap <- rbind(mbmap, c(rn, microbe))
                        }
                    } else {
                        mbmap <- rbind(mbmap, c(rn, nm))
                    }
                }
            }
            mbmap <- mbmap[mbmap[,2]!="",]
            gcnt <- table(mbmap[,2])
            if (length(mbList)!=1) {
                gcnt <- gcnt[order(gcnt, decreasing=TRUE)]
            }
            if (is.table(gcnt)) {
                ret@geneCount <- gcnt
            }
            
            incGene <- names(gcnt)[1:length(mbList)]
            mbmap <- mbmap[mbmap[,2] %in% incGene,]
            ret@geneMap <- mbmap

            vtx <- data.frame(cbind(c(mbmap[,1], mbmap[,2]), c(rep("Words",length(mbmap[,1])),
                rep("Microbes",length(mbmap[,2]))))) |> `colnames<-`(c("name","type"))
            vtx <- vtx[!duplicated(vtx),]
            
            ####
            ## It is possible not to distinguish between query microbe name 
            ## and words as in `pubmed()`, as the words are lowercases here, 
            ## and altered to titlecase after
            ####

            mbmap <- tbl_graph(nodes=vtx, edges=data.frame(mbmap), directed=FALSE)
            V(coGraph)$type <- "Words"
            coGraph <- graph_join(as_tbl_graph(coGraph),
                as_tbl_graph(mbmap))
            coGraph <- coGraph |> activate("nodes") |>
                mutate(type=ifelse(is.na(.data$Freq),"Microbes","Words"))

            ## If present, add additional graphs
            if (length(addNet)!=0) {
                for (netName in names(addNet)) {
                    tmpAdd <- addNet[[netName]]
                    coGraph <- graph_join(as_tbl_graph(coGraph),
                        as_tbl_graph(tmpAdd))
                    coGraph <- coGraph |> activate("nodes") |>
                        mutate(type=ifelse(is.na(.data$type),netName,.data$type))
                }
            }
            ## Set edge weight
            ## Probably set to NA would be better.
            E(coGraph)$edgeColor <- E(coGraph)$weight
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshMbPlot <- 0.01} else {
                corThreshMbPlot <- corThresh - 0.1
            }
            tmpW[is.na(tmpW)] <- corThreshMbPlot
            E(coGraph)$weight <- tmpW

        } else {
            coGraph <- as_tbl_graph(coGraph) |> activate("nodes") |>
                mutate(type=ifelse(is.na(.data$Freq),"Microbes","Words"))
            E(coGraph)$edgeColor <- E(coGraph)$weight
        }
        nodeN <- (coGraph |> activate("nodes") |> data.frame())$type
        V(coGraph)$nodeCat <- nodeN
        names(nodeN) <- names(V(coGraph))

        if (tag!="none") {
            netCol <- tolower(names(V(coGraph)))
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    netCol[netCol==j] <- paste0("cluster",i)
            }
            netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
            V(coGraph)$tag <- netCol

            ## Add disease and other labs
            if (!is.null(nodeN)) {
                addC <- V(coGraph)$tag
                for (nn in seq_along(names(V(coGraph)))) {
                    if (names(V(coGraph))[nn] %in% names(nodeN)) {
                        addC[nn] <- nodeN[names(V(coGraph))[nn]]
                    } else {
                        next
                    }
                }
                V(coGraph)$tag <- addC
            }
        }

        if (preserve) {
            nodeDf <- coGraph |> as_tbl_graph() |>activate("nodes") |> data.frame()
            V(coGraph)$name <- apply(nodeDf,
                  1,
                  function(x) {ifelse(x["type"]=="Words",
                  	ifelse(is.na(pdic[x["name"]]), x["name"], pdic[x["name"]]),
                    x["name"])})
            # coGraph <- set.vertex.attribute(coGraph, "name", value=newGname)
        }

        nodeName <- V(coGraph)$name
        for (i in madeUpper) {
          nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName

        if (colorize) {addFreqToMB <- TRUE}
        if (addFreqToMB) {
            ## Set pseudo freq based on min value of freq
            fre <- V(coGraph)$Freq
            fre[is.na(fre)] <- min(fre, na.rm=TRUE)
            V(coGraph)$Freq <- fre
        }


          ## Define colors

        if (tag!="none") { ## Obtain tag coloring
            if (is.null(tagPalette)) {
              cols <- V(coGraph)$tag |> unique()
              if (length(cols)>2) {
                  tagPalette <- RColorBrewer::brewer.pal(8, "Dark2")
                  tagPalette <- colorRampPalette(tagPalette)(length(cols))
              } else {
                  tagPalette <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(length(cols))]
              }
              names(tagPalette) <- cols
              tagPalette["Microbes"] <- mbColor
            }
        }

        if (is.null(catColors)) { ## Obtain category coloring
            catLen <- length(unique(V(coGraph)$nodeCat))
            if (catLen>2) {
                catColors <- RColorBrewer::brewer.pal(catLen, "Dark2")
            } else {
                catColors <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(catLen)]
            }
            names(catColors) <- unique(V(coGraph)$nodeCat)
            catColors["Microbes"] <- mbColor
        }

          if (!is.tbl_graph(coGraph)) {
              ret@igraph <- coGraph
          } else {
              ret@igraph <- as.igraph(coGraph)
          }

        ## Main plot
        netPlot <- ggraph(coGraph, layout=layout)

        netPlot <- appendEdges(netPlot, FALSE, edgeLink,
            edgeLabel, showLegend, fontFamily)

        netPlot <- appendNodesAndTexts(netPlot,tag,colorize,tagPalette,
                          showLegend,catColors,pal,fontFamily,colorText,scaleRange,
                          useSeed, ret, tagColors=tagPalette, discreteColorWord=discreteColorWord)
        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=c(1,3), name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation", na.value=naEdgeColor)+
            theme_graph()

        ret@net <- netPlot
    } else {
        ## WC
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word = names(matSorted),freq=matSorted) |>
            na.omit()

        if (tag!="none") {
            freqWords <- names(matSorted)
            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            if (!is.null(redo)){
                if (!is.null(redo@pvpick)) {
                    pvcl <- ret@pvpick
                } else {
                    if (tagWhole){
                        pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl, method.dist=tag)
                    } else {
                        pvc <- pvclust(as.matrix(dist(
                            t(
                                freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                            )
                        )), parallel=cl, method.dist=tag)
                    }
                    pvcl <- pvpick(pvc, alpha=pvclAlpha)
                    ret@pvclust <- pvc
                    ret@pvpick <- pvcl
                }
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl, method.dist=tag)
                } else {
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                        )
                    )),parallel=cl, method.dist=tag)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
            if (is.null(tagPalette)) {
              tagPalette <- RColorBrewer::brewer.pal(length(unique(pvcl$clusters)), "Dark2")
            }
            wcCol <- returnDf$word
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- tagPalette[i]
            }
            wcCol[!wcCol %in% tagPalette] <- "grey"
        }
        
        for (i in madeUpper) {
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        if (preserve) {
            for (nm in unique(returnDf$word)) {
                if (nm %in% names(pdic)) {
                    returnDf[returnDf$word == nm, "word"] <- pdic[nm]
                }
            }
        }
        if (!is.null(scaleFreq)) {
            showFreq <- returnDf$freq*scaleFreq
        } else {
            showFreq <- returnDf$freq
        }

        if (tag!="none"){
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- showFreq
            argList[["family"]] <- fontFamily
            argList[["colors"]] <- wcCol
            argList[["random.order"]] <- FALSE
            argList[["ordered.colors"]] <- TRUE

            if (useggwordcloud) {
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill="transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        } else {
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- showFreq
            argList[["family"]] <- fontFamily
            if (useggwordcloud) {
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill = "transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        }
        ret@freqDf <- returnDf
        ret@wc <- wc
    }
    return(ret)
}
