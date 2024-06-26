#' @rdname generalf
#' @export
#' @examples \dontrun{pubmed("DDX41")}
#' @import tm
#' @import GeneSummary
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @importFrom GetoptLong qqcat
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
pubmed <- function(queries, useRawQuery=FALSE,
	redo=NULL, madeUpper=c("dna","rna"),
    target="abstract", useFil=NA, filType="above",
    filNum=0, sortOrder="relevance", fontFamily="sans",
    pvclAlpha=0.95, numOnly=TRUE, delim="OR",
    apiKey=NULL, tfidf=FALSE, cl=FALSE, autoThresh=TRUE,
    pal=c("blue","red"), numWords=30, scaleRange=c(5,10),
    showLegend=FALSE, plotType="network", colorText=FALSE, quote=FALSE,
    corThresh=0.2, layout="nicely", tag="none", tagWhole=FALSE,
    onlyCorpus=FALSE, onlyTDM=FALSE, retMax=10,
    edgeLabel=FALSE, edgeLink=TRUE, ngram=1, genePlot=FALSE, scaleFreq=NULL,
    onlyDf=FALSE, tagPalette=NULL, preserve=FALSE, takeMax=FALSE,
    catColors=NULL, perQuery=FALSE,
    discreteColorWord=FALSE,
    useUdpipe=FALSE, udpipeOnlyFreq=FALSE, udpipeOnlyFreqNB=FALSE,
    addFreqToQuery=FALSE,
    naEdgeColor="grey50", cooccurrence=FALSE, colorize=FALSE,
    queryColor="grey",
    useggwordcloud=TRUE, wcScale=10, distinguish_query=TRUE, useSeed=42,
    udpipeModel="english-ewt-ud-2.5-191206.udpipe", normalize=FALSE,
    takeMean=FALSE,  absolute=TRUE,
    corOption=list(),
    deleteZeroDeg=TRUE, additionalRemove=NA, orgDb=org.Hs.eg.db,
    onlyGene=FALSE, filterByGO=FALSE, docsum=FALSE,
    pre=FALSE, onWholeDTM=FALSE, madeUpperGenes=TRUE, stem=FALSE,
    argList=list(), dateRange=NULL, cc0=FALSE)
{
    if (cc0) {
        if (requireNamespace("pubmedMini", quietly = TRUE)) {
            ## Installation by 
            ## devtools::install_github("noriakis/pubmedMini")
            pmc <- pubmedMini::loadpubmedMini()
            commons <- intersect(unique(pmc$query), queries)
            if (length(commons)<1) {
                stop("No gene query found in `pubmedMini` data.")
            }
            pmc <- pmc[pmc$query %in% commons, ]
        } else {
            stop("Please install `pubmedMini` library, by devtools::install_github('noriakis/pubmedMini')")
        }
    }

    if (useUdpipe) {
        qqcat("Using udpipe mode\n")
        plotType="network"
        ngram <- NA
        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)
    }
    if (madeUpperGenes){
        madeUpper <- c(madeUpper, tolower(keys(orgDb, "SYMBOL")))
    }
    if (pre) {
        additionalRemove <- c(additionalRemove,"genes","gene","patients","hub",
            "analysis","cells","cell","expression","doi",
            "deg","degs","author","authors","elsevier",
            "oxford","wiley")
    }
    if (is.null(redo)) {

        ret <- new("biotext")
        ret@date <- Sys.time()
        ret@type <- paste0("pubmed_",target)

        if (cc0) { ## PMC_CC0 article
            allDataDf <- pmc
            allDataDf <- allDataDf[toupper(allDataDf$section) %in% toupper(target), ]
            if (dim(allDataDf)[1]==0) {stop("No text available")}
            ret@rawText <- allDataDf
            ret@type <- paste0("PMC_CC0_",target)
        } else {
            if (useRawQuery) {
                query <- queries
            } else {
                if (perQuery) {
                    ## Disabled the limit of the number of query
                    # if (length(queries)>limit){
                    #     stop("Number of queries exceeded specified limit number")
                    # }
                    query <- queries    
                } else {
                    if ((length(queries)>10) & delim=="OR") {
                        message("Warning: major genes could dominate the search",
                        " results when the number of queries is large")
                    }
                    if (quote) {
                        query <- paste(dQuote(queries,options(useFancyQuotes = FALSE)),
                            collapse=paste0(" ",delim," "))
                    } else {
                        query <- paste(queries, collapse=paste0(" ",delim," "))
                    }
                }           
            }
            ret@query <- query
            clearQuery <- gsub('\"', '', queries)
            ret <- getPubMed(ret, query, clearQuery, type=target, apiKey=apiKey,
               retMax=retMax, sortOrder=sortOrder, perQuery=perQuery, dateRange=dateRange)
            ret@retMax <- retMax
            allDataDf <- ret@rawText            
        }
    } else {
        if (!is(redo, 'biotext')) {
            stop("Please provide biotext class object")
        }
        qqcat("Resuming from the previous results\n")
        ret <- redo
        allDataDf <- ret@rawText
    }

    if (onlyDf) {
        return(allDataDf)
    }

    ## Filter the words based on specified criteria
    if (!is.na(useFil)){
        filterWords <- retFiltWords(useFil, filType, filNum)
    } else {
        filterWords <- NA
    }
    if (length(filterWords)!=0 | length(additionalRemove)!=0){
        allfils <- c(filterWords, additionalRemove)
        allfils <- allfils[!is.na(allfils)]
        if (length(allfils)!=0) {
            ret@filtered <- allfils
        }
    }

    if (target=="abstract"){
        docs <- VCorpus(VectorSource(allDataDf$text))
    } else if (target=="title"){
        docs <- VCorpus(VectorSource(allDataDf$text))
    } else {
        stop("specify 'title' or 'abstract' to target argument")
    }
    if (preserve) {
        pdic <- preserveDict(docs, ngram, numOnly, stem)
        ret@dic <- pdic
    }
    
    docs <- makeCorpus(docs, filterWords,
        additionalRemove, numOnly, stem)
    ret@tag <- tag
    ret@corpus <- docs
    if (onlyCorpus){
        return(docs)
    }
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
  
    if (onlyTDM){
        return(docs)
    }
  
    mat <- as.matrix(docs)
    if (docsum) {
        mat <- apply(mat, 2, function(x) ifelse(x>0, 1, 0))
    }
    if (normalize) {
        mat <- sweep(mat, 2, colSums(mat), `/`)
    }

    if (takeMax & takeMean) {
        stop("Should either of specify takeMax or takeMean")
    }
    if (takeMax) {
        perterm <- apply(mat, 1, max, na.rm=TRUE)
    }
    if (takeMean) {
        perterm <- apply(mat, 1, mean, na.rm=TRUE)
    }
    if ((takeMax + takeMean)==0) {
        perterm <- rowSums(mat)
    }

    matSorted <- sort(perterm, decreasing=TRUE)
    ret@wholeFreq <- matSorted

    ret@TDM <- docs

    ## If filter by GO terms
    if (filterByGO) {
        qqcat("`filterByGO` option enabled. Filtering by GO terms ...\n")
        data_env <- new.env(parent = emptyenv())
        load(system.file("extdata", "sysdata.rda", package = "biotextgraph"),
            envir=data_env)
        goWords <- data_env[["goWords"]]
        goWords2gram <- data_env[["goWords2gram"]]

        if (ngram==1) {
            filtered_by_GO <- names(matSorted)[tolower(names(matSorted)) %in% goWords]
            matSorted <- matSorted[filtered_by_GO]
        } else if (ngram==2) {
            filtered_by_GO <- names(matSorted)[tolower(names(matSorted)) %in% goWords2gram]
            matSorted <- matSorted[filtered_by_GO]          
        } else {# Do nothing
        }
    }

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }
    ret@numWords <- numWords
    

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word=names(matSorted),freq=matSorted)  |>
            na.omit()
        for (i in madeUpper) {
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        ret@freqDf <- returnDf
    
        freqWords <- names(matSorted)

        if (useUdpipe) {
            if (udpipeOnlyFreq & udpipeOnlyFreqNB) {
                stop("Cannot specify both of these options")
            }
            if (udpipeOnlyFreq) {
                showNeighbors <- NULL
                showFreq <- freqWords
            }
            if (udpipeOnlyFreqNB) {
                showNeighbors <- freqWords
                showFreq <- NULL
            }
            ret <- retUdpipeNet(ret=ret,texts=allDataDf,
                udmodel_english=udmodel_english,
                orgDb=orgDb, filterWords=filterWords,
                additionalRemove=additionalRemove,
                colorText=colorText,edgeLink=edgeLink,
                queryPlot=genePlot, layout=layout,
                pal=pal,showNeighbors=showNeighbors,
                showFreq=showFreq, nodePal=tagPalette)
            ret@type <- "udpipe"
            return(ret)
        }


        DTM <- t(as.matrix(docs))
        row.names(DTM) <- allDataDf$query

        if (tag!="none") {
            if (!is.null(ret) & length(ret@pvpick)!=0){
                qqcat("Using previous pvclust results")
                pvcl <- ret@pvpick
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(
                            dist(as.matrix(docs))
                        ), parallel=cl, method.dist=tag)
                } else {
          
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            DTM[, colnames(DTM) %in% freqWords]
                        ))), parallel=cl, method.dist=tag)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
        }
    
        matrixs <- obtainMatrix(ret, FALSE, NULL, DTM, freqWords,
            corThresh, cooccurrence, onWholeDTM, numWords, autoThresh,
            absolute, corOption)
    

        coGraph <- matrixs$coGraph
        ret <- matrixs$ret
        
        # if (tag=="cor") {
		#     ret <- tag_words(ret, cl,
		# 	    pvclAlpha, whole=tagWhole,
		# 	    num_words=ret@numWords,
		# 	    corMat=TRUE, mat=matrixs$ret@corMat)
        #     pvc <- ret@pvclust
        #     pvcl <- ret@pvpick
        # }

        ret@igraphRaw <- coGraph
        coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]



        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }
    
        nodeName <- V(coGraph)$name

        incQuery <- NULL
        distinguish_query <- FALSE
        if (genePlot) {
            genemap <- c()
            for (rn in nodeName){
                tmp <- DTM[ ,rn]
                for (nm in names(tmp[tmp!=0])){
                    if (nm!=""){
                        if (grepl(",",nm,fixed=TRUE)){
                            for (nm2 in unlist(strsplit(nm, ","))){
                                if (distinguish_query) {
                                    genemap <- rbind(genemap,
                                        c(rn, paste(nm2, "(Q)")))
                                } else {
                                    genemap <- rbind(genemap, c(rn, nm2))
                                }
                            }
                        } else {
                            if (distinguish_query) {
                                genemap <- rbind(genemap,
                                    c(rn, paste(nm, "(Q)")))
                            } else {
                                genemap <- rbind(genemap, c(rn, nm))
                            }
                        }
                    }
                }
            }
            incQuery <- unique(genemap[,2])
            retGenemap <- genemap
            ret@geneMap <- retGenemap

            vtx <- data.frame(cbind(c(genemap[,2], genemap[,1]),
                c(rep("Queries",length(genemap[,2])),
                rep("Words",length(genemap[,1]))))) |> 
                `colnames<-`(c("name","type"))
            vtx <- vtx[!duplicated(vtx),]


            vtx <- vtx |> `rownames<-`(1:nrow(vtx))
            eds <- data.frame(genemap)
            words <- vtx |> subset(vtx$type=="Words")
            queriesDf <- vtx |> subset(vtx$type=="Queries")

            row.names(words)[which(words$name %in% eds[,1])]
            row.names(queriesDf)[which(queriesDf$name %in% eds[,2])]

            eds[,1] <- sapply(eds[,1], function(x) {
                as.integer(row.names(words)[which(words$name %in% x)])
            })

            eds[,2] <- sapply(eds[,2], function(x) {
                as.integer(row.names(queriesDf)[which(queriesDf$name %in% x)])
            })

            genemap <- tbl_graph(nodes=vtx,edges=eds,directed=FALSE)
            V(coGraph)$type <- "Words"
            coGraph <- graph_join(as_tbl_graph(coGraph), genemap)

            E(coGraph)$edgeColor <- E(coGraph)$weight
            tmpW <- E(coGraph)$weight

            ## Replace NA edges with slightly thin edge
            if (corThresh < 0.1) {
                corThreshGenePlot <- 0.01
            } else {
                corThreshGenePlot <- corThresh - 0.1
            }
            tmpW[is.na(tmpW)] <- corThreshGenePlot
            E(coGraph)$weight <- tmpW
        } else { ## If not genePlot
            coGraph <- as_tbl_graph(coGraph) |> activate("nodes") |>
                mutate(type=ifelse(is.na(.data$Freq),"Queries","Words"))
            E(coGraph)$edgeColor <- E(coGraph)$weight
        }

        ## Assign node category
        nodeN <- V(coGraph)$type
        V(coGraph)$nodeCat <- nodeN
        names(nodeN) <- V(coGraph)$name
    
        if (tag!="none") {
            netCol <- tolower(names(V(coGraph)))
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    netCol[netCol==j] <- paste0("cluster",i)
            }
            netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
            V(coGraph)$tag <- netCol
            if (!is.null(nodeN)) {
                addC <- V(coGraph)$tag
                for (nn in seq_along(names(V(coGraph)))) {
                    if (V(coGraph)$tag[nn]!="not_assigned"){next}
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
            nodeDf <- coGraph |>
                as_tbl_graph() |>
                activate("nodes") |>
                data.frame()
            V(coGraph)$name <- apply(nodeDf, 1,
                function(x) {
                    ifelse(x["type"]=="Words",
              	        ifelse(is.na(pdic[x["name"]]),
                            x["name"], pdic[x["name"]]),
                        x["name"])
                }
            )
        }

        nodeName <- V(coGraph)$name
        for (i in madeUpper) {
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName

        ## Set pseudo freq as min value of freq

        if (addFreqToQuery) {
            fre <- V(coGraph)$Freq
            fre[is.na(fre)] <- min(fre, na.rm=TRUE)
            V(coGraph)$Freq <- fre
        }

        if (colorize) {
            fre <- V(coGraph)$Freq
            fre[is.na(fre)] <- min(fre, na.rm=TRUE)
            V(coGraph)$Freq <- fre
        }
        
        ## Main plot

        if (onlyGene) {
            qqcat("Subsetting to the gene symbol in orgDb\n")
            included <- names(V(coGraph))[
                tolower(names(V(coGraph))) %in% tolower(keys(orgDb, "SYMBOL"))
            ]
            coGraph <- induced_subgraph(coGraph, included)
        }


        if (!is.tbl_graph(coGraph)) {
            ret@igraph <- coGraph
        } else {
            ret@igraph <- as.igraph(coGraph)
        }

        netPlot <- ggraph(coGraph, layout=layout)
        netPlot <- appendEdges(netPlot, FALSE, edgeLink,
            edgeLabel, showLegend, fontFamily)

        if (tag!="none") { ## Obtain tag coloring
            if (is.null(tagPalette)) {
                cols <- V(coGraph)$tag |> unique()
                if (length(cols)>2) {
                    tagPalette <- RColorBrewer::brewer.pal(8, "Dark2")
                    tagPalette <- colorRampPalette(tagPalette)(length(cols))
                } else {
                    tagPalette <- RColorBrewer::brewer.pal(3,"Dark2")[
                        seq_len(length(cols))
                    ]
                }
                names(tagPalette) <- cols
                tagPalette["Queries"] <- queryColor
            }
        }

        if (is.null(catColors)) { ## Obtain category coloring
            catLen <- length(unique(V(coGraph)$nodeCat))
            if (catLen>2) {
                catColors <- RColorBrewer::brewer.pal(catLen, "Dark2")
            } else {
                catColors <- RColorBrewer::brewer.pal(3,"Dark2")[
                    seq_len(catLen)
                ]
            }
            names(catColors) <- unique(V(coGraph)$nodeCat)
            catColors["Queries"] <- queryColor
        }
        
        netPlot <- appendNodesAndTexts(netPlot,tag,colorize,tagPalette,
            showLegend,catColors,pal,fontFamily,colorText,scaleRange,
            useSeed,ret,tagColors=tagPalette,
            discreteColorWord=discreteColorWord)


        netPlot <- netPlot+
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=c(1,3), name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                na.value=naEdgeColor,
                name = "Correlation")+
            theme_graph()
        ret@net <- netPlot
    } else {
        ## WC part
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word = names(matSorted),freq=matSorted) |>
            na.omit()
        freqWords <- names(matSorted)
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ])) 
        if (tag!="none") {
            if (!is.null(redo) & length(redo@pvpick)!=0) {
                qqcat("Using previous pvclust results")
                pvcl <- redo@pvpick
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(as.matrix(docs))),
                        parallel=cl, method.dist=tag)
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
            wcCol <- returnDf$word
	        if (is.null(tagPalette)) {
	        	tagPalette <- colorRampPalette(brewer.pal(11, "RdBu"))(length(pvcl$clusters |> unique()))
	        	names(tagPalette) <- pvcl$clusters |> unique()
	        }
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- tagPalette[i]
            }
            wcCol[!wcCol %in% tagPalette] <- "grey"
        }
        for (i in madeUpper) {
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