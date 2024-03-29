#' pubmed
#' 
#' make word cloud or correlation network from PubMed
#' 
#' @param queries gene symbols
#' @param redo if plot in other parameters, input the previous list
#' @param madeUpper make the words uppercase in resulting plot
#' @param madeUpperGenes make genes upper case automatically (default to TRUE)
#' @param pal palette for color gradient in correlation network
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param cooccurrence default to FALSE, if TRUE, use cooccurrence instead of correlation
#' @param corThresh the correlation threshold
#' @param autoThresh automatically determine the threshold value to show `numWords`
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @param colorText color text label based on frequency in correlation network
#' @param ngram default to NA (1)
#' @param additionalRemove specific words to be excluded
#' @param target "abstract" or "title"
#' @param tag perform pvclust on words and colorlize them in wordcloud or network
#' argument of "cor" or "tdm". Default to "none", which performs no tagging.
#' If wordcloud, tagging will be performed on TDM.
#' @param tagWhole tag based on whole data or subset
#' @param genePlot plot associated genes (default: FALSE)
#' Query gene name is shown with (Q)
#' @param useFil filter based on "GS_TfIdf" (whole gene summary tf-idf)
#'  or "BSDB_TfIdf" (whole bugsigdb tf-idf)
#' @param filNum specify filter tfidf
#' @param filType "above" or "below"
#' @param apiKey api key for eutilities
#' @param tfidf use TfIdf when making TDM
#' @param pvclAlpha alpha for pvpick()
#' @param onlyCorpus return only corpus
#' @param onlyTDM return only TDM
#' @param pre filter preset words like publisher's names
#' @param preserve try to preserve the original characters
#' @param numOnly delete number only
#' @param cl cluster to pass to pvclust (snow::makeCluster(n))
#' @param bn perform bootstrap-based Bayesian network inference 
#' instead of correlation using bnlearn
#' @param R how many bootstrap when bn is stated
#' @param delim delimiter for queries
#' @param retMax how many items are to be retlieved?
#' @param orgDb org database, default to org.Hs.eg.db
#' @param quote whether to quote the queries
#' @param onWholeDTM calculate correlation network
#'                   on whole dataset or top-words specified by numWords
#' @param limit limit number for query count, default to 10
#' @param sortOrder sort order, passed to rentrez function
#' @param stem whether to use stemming
#' @param onlyDf return only the raw data.frame of searching PubMed
#' @param tagPalette tag palette when tag is TRUE
#' @param takeMax when summarizing term-document matrix, take max.
#' Otherwise take sum.
#' @param onlyGene plot only the gene symbol
#' (orgDb with SYMBOL key can be used)
#' @param argList parameters to pass to wordcloud()
#' @param useUdpipe use udpipe to make a network
#' @param udpipeModel udpipe model file name
#' @param udpipeOnlyFreq when using udpipe, include only high-frequent words
#' @param udpipeOnlyFreqNB when using udpipe, include only the neighbors of
#' high-frequent words
#' @param normalize sum normalize the term frequency document-wise
#' @param takeMean take mean values for each term in term-document matrix
#' @param naEdgeColor edge color linking query with the other category than text
#' @param useggwordcloud default to TRUE
#' @param wcScale scaling size for ggwordcloud
#' @param fontFamily font family to use, default to "sans"
#' @param distinguish_query if TRUE, distinguish query and returned texts
#' by appending (Q) on query
#' @param colorize color the nodes and texts based on their category
#' @param queryColor color for associated queries with words
#' @param useSeed use seed
#' @param scaleFreq scale the frequency
#' @param addFreqToQuery add pseudo-frequency to query node
#' @param discreteColorWord colorize words by "Words" category, not frequency.
#' @param catColors colors for words ant texts when colorize=TRUE and discreteColorWord is TRUE
#' @export
#' @examples \dontrun{pubmed("DDX41")}
#' @return object consisting of data frame and ggplot2 object
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
pubmed <- function(queries, redo=NULL, madeUpper=c("dna","rna"),
   target="abstract", useFil=NA, filType="above",
   filNum=0, sortOrder="relevance", fontFamily="sans",
   pvclAlpha=0.95, numOnly=TRUE, delim="OR", limit=10,
   apiKey=NULL, tfidf=FALSE, cl=FALSE, autoThresh=TRUE,
   pal=c("blue","red"), numWords=30, scaleRange=c(5,10),
   showLegend=FALSE, plotType="network", colorText=FALSE, quote=FALSE,
   corThresh=0.2, layout="nicely", tag="none", tagWhole=FALSE,
   onlyCorpus=FALSE, onlyTDM=FALSE, bn=FALSE, R=20, retMax=10,
   edgeLabel=FALSE, edgeLink=TRUE, ngram=1, genePlot=FALSE, scaleFreq=NULL,
   onlyDf=FALSE, tagPalette=NULL, preserve=TRUE, takeMax=FALSE,
   catColors=NULL,
   discreteColorWord=FALSE,
   useUdpipe=FALSE, udpipeOnlyFreq=FALSE, udpipeOnlyFreqNB=FALSE,
   addFreqToQuery=FALSE,
   naEdgeColor="grey50", cooccurrence=FALSE, colorize=FALSE,
   queryColor="grey",
   useggwordcloud=TRUE, wcScale=10, distinguish_query=TRUE, useSeed=42,
   udpipeModel="english-ewt-ud-2.5-191206.udpipe", normalize=FALSE,
   takeMean=FALSE,
   deleteZeroDeg=TRUE, additionalRemove=NA, orgDb=org.Hs.eg.db,
   onlyGene=FALSE,
   pre=FALSE, onWholeDTM=FALSE, madeUpperGenes=TRUE, stem=FALSE,
   argList=list())
{
	if (!tag %in% c("none","tdm","cor")) {
		stop("tag should be none, tdm, or cor.")
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
        ret@type <- paste0("pubmed_",target)
        ## Disabled the limit of the number of query
        # if (length(queries)>limit){
        #     stop("Number of queries exceeded specified limit number")
        # }
        # ret@query <- queries
        # ret@delim <- delim
        if (quote) {
            query <- paste(dQuote(queries,options(useFancyQuotes = FALSE)),
                collapse=paste0(" ",delim," "))
        } else {
            query <- paste(queries, collapse=paste0(" ",delim," "))
        }
        ret@query <- query
        clearQuery <- gsub('\"', '', queries)
        ret <- getPubMed(ret, query, clearQuery, type=target, apiKey=apiKey,
           retMax=retMax, sortOrder=sortOrder)
        ret@retMax <- retMax
        allDataDf <- ret@rawText
    } else {
        if (class(redo)!="biotext") {
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
        stop("specify target or abstract")
    }
    if (preserve) {
        pdic <- preserveDict(docs, ngram, numOnly, stem)
        ret@dic <- pdic
    }
    
    docs <- makeCorpus(docs, filterWords,
        additionalRemove, numOnly, stem)
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

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }
    ret@numWords <- numWords  
  
    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word=names(matSorted),freq=matSorted)
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
            return(ret)
        }


        DTM <- t(as.matrix(docs))
        row.names(DTM) <- allDataDf$query

        if (tag=="tdm") {
            if (!is.null(ret) & length(ret@pvpick)!=0){
                qqcat("Using previous pvclust results")
                pvcl <- ret@pvpick
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(
                            dist(as.matrix(docs))
                        ), parallel=cl)
                } else {
          
                    pvc <- pvclust(as.matrix(dist(
                        t(
                            DTM[, colnames(DTM) %in% freqWords]
                        ))), parallel=cl)
                }
                pvcl <- pvpick(pvc, alpha=pvclAlpha)
                ret@pvclust <- pvc
                ret@pvpick <- pvcl
            }
        }
    
        matrixs <- obtainMatrix(ret, FALSE, NULL, DTM, freqWords,
            corThresh, cooccurrence, onWholeDTM, numWords, autoThresh)
    

        coGraph <- matrixs$coGraph
        ret <- matrixs$ret
        
        if (tag=="cor") {
		    ret <- tag_words(ret, cl,
			    pvclAlpha, whole=tagWhole,
			    num_words=ret@numWords,
			    corMat=TRUE, mat=matrixs$ret@corMat)
            pvc <- ret@pvclust
            pvcl <- ret@pvpick
        }

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
        netPlot <- appendEdges(netPlot, bn, edgeLink,
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
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
        freqWords <- names(matSorted)
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ])) 
        if (tag!="none") {
            if (!is.null(redo) & length(redo@pvpick)!=0) {
                qqcat("Using previous pvclust results")
                pvcl <- redo@pvpick
            } else {
                if (tagWhole){
                    pvc <- pvclust(as.matrix(dist(as.matrix(docs))),
                        parallel=cl)
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
            wcCol <- returnDf$word
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