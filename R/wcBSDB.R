#' wcBSDB
#' 
#' Visualize BugSigDB
#' 
#' @param mbList microbe list
#' @param additionalRemove specific words to be excluded
#' @param exclude "frequency" or "tfidf",
#' @param excludeFreq default to 1000 (no words are excluded)
#' @param excludeType ">" or "<"
#' @param mbPlot plot microbe names
#' @param disPlot plot diseases
#' @param madeUpper make the words uppercase in resulting plot
#' @param pal palette for color gradient in correlation network and tag coloring
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param corThresh the correlation threshold
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @param colorText color text label based on frequency in correlation network
#' @param ngram default to NA (1)
#' @param tag perform pvclust on words and colorlize them in wordcloud
#' @param tagWhole tag whole set or subset
#' @param target "title" or "abstract"
#' @param cl cluster object passed to pvclust
#' @param apiKey api key for eutilities
#' @param redo if "abstract" is chosen in target, one can provide resulting object again
#' @param pre predefined filter words
#' @param tfidf use TfIdf when making TDM
#' @param pvclAlpha alpha for pvpick()
#' @param numOnly delete number only
#' @param preserve try to preserve the original characters
#' @param onlyTDM return only TDM
#' @param onWholeDTM calculate correlation network
#'                   on whole dataset or top-words specified by numWords
#' @param metab tibble of metabolite - taxon association
#' @param metThresh threshold of association
#' @param metCol metabolite data frame column name in the order of
#' "candidate taxon", "metabolite", "quantitative values for thresholding"
#' @param stem whether to use stemming
#' @param curate include articles in bugsigdb
#' @param abstArg passed to PubMed function when using curate=FALSE
#' @param nodePal node palette when tag is TRUE
#' @param ecPlot plot link between enzyme and microbes
#' this option requires two files to be passed to wcEC() and getUPTax().
#' @param ecFile enzyme database file
#' @param upTaxFile UniProt taxonomy file
#' @param takeMax when summarizing term-document matrix, take max.
#' @param filterMax use pre-calculated filter based on max-values when excluding TfIdf
#' Otherwise take sum.
#' @param useUdpipe use udpipe to make a network
#' @param udpipeModel udpipe model file name
#' @param argList parameters to pass to wordcloud()
#' @param normalize sum normalize the term frequency document-wise
#' @param takeMean take mean values for each term in term-document matrix
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import bugsigdbr
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @import XML
#' @importFrom dplyr filter
#' @importFrom stats dist
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency graph_from_edgelist
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' 
#' @examples
#' mbList <- c("Veillonella dispar")
#' wcBSDB(mbList)
#' @export
#' 
wcBSDB <- function (mbList,
                    excludeFreq=1000, exclude="frequency",
                    excludeType=">", normalize=FALSE, takeMean=FALSE,
                    additionalRemove=NA, tfidf=FALSE,
                    target="title", apiKey=NULL, takeMax=FALSE,
                    pre=FALSE, pvclAlpha=0.95, numOnly=TRUE,
                    madeUpper=c("dna","rna"), redo=NULL,
                    pal=c("blue","red"), numWords=15, preserve=TRUE,
                    metab=NULL, metThresh=0.2, curate=TRUE,
                    abstArg=list(), nodePal=palette(), metCol=NULL,
                    scaleRange=c(5,10), showLegend=FALSE, ecPlot=FALSE,
                    edgeLabel=FALSE, mbPlot=FALSE, onlyTDM=FALSE,
                    ecFile=NULL, upTaxFile=NULL, filterMax=FALSE,
                    useUdpipe=FALSE,
                    udpipeModel="english-ewt-ud-2.5-191206.udpipe",
                    ngram=NA, plotType="wc", disPlot=FALSE, onWholeDTM=FALSE,
                    colorText=FALSE, corThresh=0.2, tag=FALSE, tagWhole=FALSE, stem=FALSE,
                    layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, cl=FALSE, argList=list()) {

    if (useUdpipe) {
        qqcat("Using udpipe mode\n")
        plotType="network"
        ngram <- NA

        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)
    }
    ret <- new("osplot")
    ret@query <- mbList
    ret@type <- paste0("BSDB_",target)
    addNet <- list()

    if (pre) {additionalRemove <- c(additionalRemove, "microbiota",
        "microbiome","relative","abundance","abundances",
        "including","samples","sample","otu","otus","investigated",
        "taxa","taxon")}
    if (is.null(redo)) {
        qqcat("Input microbes: @{length(mbList)}\n")
        tb <- importBugSigDB()
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
        # returnList[["subsetBSDB"]] <- subTb

        titles <- unique(subTb$Title)
        titles <- titles[!is.na(titles)]
        fil <- c()
        for (title in titles){
            tmp <- subset(subTb, Title==title)
            if (dim(tmp)[1]>1){
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
        ret@rawText <- subTb

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
                            tax <- unique(subset(fil, ID==pmid)$query)
                            abstDf <- rbind(abstDf, c(pmid, text, tax))
                        }
                    }
                }
                abstDf <- data.frame(abstDf) |> `colnames<-`(c("ID","text","query"))
                # abstset <- xmlElementsByTagName(parsedXML$doc$children$PubmedArticleSet,
                #                                 "Abstract",
                #                                 recursive = TRUE)
                # absttext <- as.character(xmlValue(abstset))
                ret@rawTextBSDB <- abstDf
            } else {
                abstDf <- ret@rawTextBSDB
                filterWords <- ret@filtered
                subTb <- ret@rawText
            }
            docs <- VCorpus(VectorSource(abstDf$text))
        } else {
            docs <- VCorpus(VectorSource(fil$text))
        }
    } else {
        abstDf <- do.call(wcAbst, c(list(queries=mbList,
            quote=TRUE, target=target, onlyDf=TRUE), abstArg))
        ret@rawTextBSDB <- abstDf
        docs <- VCorpus(VectorSource(abstDf$text))
    }
    ## Make corpus
    ## Filter high frequency words
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

    if (disPlot & plotType=="network") {## This does not need to be deduplicated
        dis <- NULL
        for (i in seq_len(nrow(subTb))){
            dis <- rbind(dis,
            c(subTb[i, "Condition"], subTb[i, "query"]))
        }
        dis <- dis[!is.na(dis[,1]),]
        dis <- dis[!is.na(dis[,2]),]
        dmg <- simplify(graph_from_data_frame(dis, directed=FALSE))
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
            metabGraph <- simplify(graph_from_edgelist(metabGraph,
                directed=FALSE))
            addNet[["Metabolites"]] <- metabGraph
        }
    }

    ## Add EC if present
    if (ecPlot & plotType=="network") {
        mbPlot <- TRUE
        if (is.null(ecFile)) {stop("Please provide EC file")}
        if (is.null(upTaxFile)) {stop("Please provide UniProt taxonomy file")}
        ecDf <- wcEC(file=ecFile, ecnum="all", taxec=TRUE,
            taxFile=upTaxFile, candTax=mbList)
        if (!is.null(ecDf)) {
           ecg <- simplify(graph_from_data_frame(ecDf[,c("desc","query")], 
                directed=FALSE))
           addNet[["Enzymes"]] <- ecg
        }
    }


    if (useUdpipe) {
        if (curate & target!="abstract") {
            abstDf <- fil
        }
        ret <- retUdpipeNet(ret=ret,texts=abstDf,udmodel_english=udmodel_english,
          orgDb=orgDb, filterWords=filterWords, additionalRemove=additionalRemove,
          colorText=colorText,edgeLink=edgeLink,queryPlot=mbPlot, layout=layout,
          pal=pal,showNeighbors=NULL, showFreq=NULL, nodePal=nodePal,addNet=addNet)
        ## TODO: add disPlot and the other options
        return(ret)
    }


    ## Set parameters for correlation network
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    if (!is.na(ngram)){
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
        # freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        ## TODO: before or after?
        freqWordsDTM <- t(as.matrix(docs))
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
        ret@freqDf <- returnDf

        if (tag) {# Needs rework
            if (!is.null(redo)){
                if (!is.null(redo@pvpick)) {
                    pvcl <- ret@pvpick
                } else {
                    if (tagWhole){
                        pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl)
                    } else {
                        pvc <- pvclust(as.matrix(dist(
                            t(
                                freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                                )
                            )), parallel=cl)
                    }
                    pvcl <- pvpick(pvc, alpha=pvclAlpha)
                    ret@pvclust <- pvc
                    ret@pvpick <- pvcl
                }
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl)
                } else {
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                            )
                        )),parallel=cl)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
        }

        ## genePlot: plot associated genes
        if (disPlot) {mbPlot <- TRUE}
        if (!is.null(metab)) {mbPlot <- TRUE}
        if (ecPlot) {mbPlot <- TRUE}

        # if (target=="abstract" & mbPlot) {stop("mbPlot is not supported in abstract as target.")}
        if (mbPlot) {
            if (target=="abstract") {
                row.names(freqWordsDTM) <- abstDf$query
            } else {
                if (curate) {
                    row.names(freqWordsDTM) <- fil$query
                } else {
                    row.names(freqWordsDTM) <- abstDf$query
                }
            }
        }

        ## Check correlation
        if (onWholeDTM){
            corData <- cor(freqWordsDTM)
        } else {
            corData <- cor(freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords])
        }
        ret@corMat <- corData
        ret@corThresh <- corThresh

        ## Set correlation below threshold to zero
        corData[corData<corThresh] <- 0
        coGraph <- graph.adjacency(corData, weighted=TRUE,
                    mode="undirected", diag = FALSE)
        ## before or after?
        coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }

        nodeName <- V(coGraph)$name
        dtmCol <- colnames(freqWordsDTM)
        for (i in madeUpper) {
            dtmCol[dtmCol == i] <- toupper(i)
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        colnames(freqWordsDTM) <- dtmCol

        if (mbPlot) {
            mbmap <- c()
            for (rn in nodeName){
                tmp <- freqWordsDTM[ ,rn]
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

            candMb <- unique(mbmap[,2])
            nodeN <- rep("Microbes",length(candMb))
            names(nodeN) <- candMb

            mbmap <- simplify(graph_from_edgelist(mbmap,
                directed = FALSE))
            coGraph <- igraph::union(coGraph, mbmap)
            
            ## If present, add additional edges
            if (length(addNet)!=0) {
                for (netName in names(addNet)) {
                    tmpAdd <- addNet[[netName]]
                    tmpNN <- names(V(tmpAdd))
                    tmpNN <- tmpNN[!tmpNN %in% names(nodeN)]

                    newNN <- rep(netName, length(tmpNN))
                    names(newNN) <- tmpNN
                    nodeN <- c(nodeN, newNN)

                    coGraph <- igraph::union(coGraph, tmpAdd)
                }
            }

            ## Set edge weight
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshMbPlot <- 0.01} else {
                corThreshMbPlot <- corThresh - 0.1
            }
            tmpW[is.na(tmpW)] <- corThreshMbPlot
            E(coGraph)$weight <- tmpW

        }



        if (tag) {
            netCol <- tolower(names(V(coGraph)))
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    netCol[netCol==j] <- paste0("cluster",i)
            }
            netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
            V(coGraph)$tag <- netCol

            ## Add disease and other labs
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

        if (preserve) {
            newGname <- NULL
            for (nm in names(V(coGraph))) {
              if (nm %in% names(pdic)) {
                newGname <- c(newGname, pdic[nm])
              } else {
                newGname <- c(newGname, nm)
                  }
            }
            coGraph <- set.vertex.attribute(coGraph, "name", value=newGname)
        }


        ## Add pseudo-freq
        naRm <- V(coGraph)$Freq
        naRm[is.na(naRm)] <- min(naRm, na.rm=TRUE)
        V(coGraph)$Freq <- naRm

        ret@igraph <- coGraph

        ## Main plot
        netPlot <- ggraph(coGraph, layout=layout)


        if (edgeLink){
            if (edgeLabel){
                netPlot <- netPlot +
                            geom_edge_link(
                                aes(width=weight,
                                color=weight,
                                label=round(weight,3)),
                                angle_calc = 'along',
                                label_dodge = unit(2.5, 'mm'),
                                alpha=0.5,
                                show.legend = showLegend)
            } else {
                netPlot <- netPlot +
                            geom_edge_link(aes(width=weight, color=weight),
                                alpha=0.5, show.legend = showLegend)
            }
        } else {
            if (edgeLabel){
                netPlot <- netPlot +
                            geom_edge_diagonal(
                                aes(width=weight,
                                color=weight,
                                label=round(weight,3)),
                                angle_calc = 'along',
                                label_dodge = unit(2.5, 'mm'),
                                alpha=0.5,
                                show.legend = showLegend)
            } else {
                netPlot <- netPlot +
                            geom_edge_diagonal(aes(width=weight, color=weight),
                                alpha=0.5, show.legend = showLegend)                
            }
        }

        if (tag) {
            netPlot <- netPlot + geom_node_point(aes(size=Freq, color=tag),
                                                show.legend = showLegend) +
                                 scale_color_manual(values=nodePal)
        } else { 
            netPlot <- netPlot + geom_node_point(aes(size=Freq, color=Freq),
                                                show.legend = showLegend)+
                                 scale_color_gradient(low=pal[1],high=pal[2],
                                                      name = "Frequency")
        }

        if (colorText){
            if (tag) {
                netPlot <- netPlot + 
                    geom_node_text(aes(label=name, size=Freq, color=tag),
                        check_overlap=TRUE, repel=TRUE,# size = labelSize,
                        bg.color = "white", segment.color="black",
                        bg.r = .15, show.legend=showLegend)
            } else {
                netPlot <- netPlot + 
                    geom_node_text(aes(label=name, size=Freq, color=Freq),
                        check_overlap=TRUE, repel=TRUE,# size = labelSize,
                        bg.color = "white", segment.color="black",
                        bg.r = .15, show.legend=showLegend)
            }
        } else {
            netPlot <- netPlot +
                        geom_node_text(aes(label=name, size=Freq),
                            check_overlap=TRUE, repel=TRUE,# size = labelSize,
                            color = "black",
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend) 
        }
        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=c(1,3), name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation")+
            theme_graph()

        ret@net <- netPlot
    } else {
        ## WC
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)

        if (tag) {
            freqWords <- names(matSorted)
            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            if (!is.null(redo)){
                if (!is.null(redo@pvpick)) {
                    pvcl <- ret@pvpick
                } else {
                    if (tagWhole){
                        pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl)
                    } else {
                        pvc <- pvclust(as.matrix(dist(
                            t(
                                freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                            )
                        )), parallel=cl)
                    }
                    pvcl <- pvpick(pvc, alpha=pvclAlpha)
                    ret@pvclust <- pvc
                    ret@pvpick <- pvcl
                }
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),parallel=cl)
                } else {
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]
                        )
                    )),parallel=cl)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
            wcCol <- returnDf$word
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- pal[i]
            }
            wcCol[!wcCol %in% pal] <- "grey"
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
        if (tag){
            wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, 
                                               freq = returnDf$freq,
                                               colors = wcCol,
                                               random.order=FALSE,
                                               ordered.colors=TRUE)))
        } else {
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- returnDf$freq
            wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
        }
        ret@freqDf <- returnDf
        ret@wc <- wc
    }
    return(ret)
}
