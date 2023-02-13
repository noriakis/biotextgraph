#' wcAbst
#' 
#' make word cloud or correlation network from PubMed
#' 
#' 
#' @param queries gene symbols
#' @param redo if plot in other parameters, input the previous list
#' @param madeUpper make the words uppercase in resulting plot
#' @param madeUpperGenes make genes upper case automatically (default to TRUE)
#' @param pal palette for color gradient in correlation network
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
#' @param additionalRemove specific words to be excluded
#' @param target "abstract" or "title"
#' @param tag cluster the words based on text using pvclust
#' @param tagWhole tag based on whole data or subset
#' @param genePlot plot associated genes (default: FALSE)
#' Query gene name is shown with (Q)
#' @param useFil filter based on "GS_TfIdf" (whole gene summary tf-idf)
#'  or "BSDB_TfIdf" (whole bugsigdb tf-idf)
#' @param filNum specify filter tfidf
#' @param filType "above" or "below"
#' @param geneUpper make queries uppercase
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
#' @param nodePal node palette when tag is TRUE
#' @param takeMax when summarizing term-document matrix, take max.
#' Otherwise take sum.
#' @param onlyGene plot only the gene symbol
#' (orgDb with SYMBOL key can be used)
#' @param argList parameters to pass to wordcloud()
#' @param useUdpipe use udpipe to make a network
#' @param udpipeModel udpipe model file name
#' @param udpipeOnlyFreq when using udpipe, include only high-frequent words
#' @param udpipeOnlyFreqN when using udpipe, include only the neighbors of
#' high-frequent words
#' @param normalize sum normalize the term frequency document-wise
#' @param takeMean take mean values for each term in term-document matrix
#' @param naEdgeColor edge color linking query with the other category than text
#' @export
#' @examples \donttest{wcAbst("DDX41")}
#' @return object consisting of data frame and ggplot2 object
#' @import tm
#' @import GeneSummary
#' @import wordcloud
#' @import igraph
#' @import stringr
#' @import ggraph ggplot2
#' @importFrom GetoptLong qqcat
#' @importFrom AnnotationDbi keys
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
wcAbst <- function(queries, redo=NULL, madeUpper=c("dna","rna"),
                   target="abstract", useFil=NA, filType="above",
                   filNum=0, sortOrder="relevance",
                   pvclAlpha=0.95, numOnly=TRUE, delim="OR", limit=10,
                   geneUpper=TRUE, apiKey=NULL, tfidf=FALSE, cl=FALSE,
                   pal=c("blue","red"), numWords=30, scaleRange=c(5,10),
                   showLegend=FALSE, plotType="wc", colorText=FALSE, quote=FALSE,
                   corThresh=0.2, layout="nicely", tag=FALSE, tagWhole=FALSE,
                   onlyCorpus=FALSE, onlyTDM=FALSE, bn=FALSE, R=20, retMax=10,
                   edgeLabel=FALSE, edgeLink=TRUE, ngram=NA, genePlot=FALSE,
                   onlyDf=FALSE, nodePal=palette(), preserve=TRUE, takeMax=FALSE,
                   useUdpipe=FALSE, udpipeOnlyFreq=FALSE, udpipeOnlyFreqN=FALSE,
                   naEdgeColor="grey50",
                   udpipeModel="english-ewt-ud-2.5-191206.udpipe", normalize=FALSE, takeMean=FALSE,
                   deleteZeroDeg=TRUE, additionalRemove=NA, orgDb=org.Hs.eg.db, onlyGene=FALSE,
                   pre=FALSE, onWholeDTM=FALSE, madeUpperGenes=TRUE, stem=FALSE, argList=list())
{
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
    ret <- new("osplot")
    ret@type <- paste0("pubmed_",target)
    if (length(queries)>limit){
      stop("Number of queries exceeded specified limit number")}
    # ret@query <- queries
    # ret@delim <- delim
    if (quote) {
      query <- paste(dQuote(queries,options(useFancyQuotes = FALSE)), collapse=paste0(" ",delim," "))
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
    qqcat("Resuming from the previous results\n")
    ret <- redo
    allDataDf <- ret@rawText
  }

  if (onlyDf) {
    return(allDataDf)
  }

  
  ## Probably set default filnum?
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

  if (geneUpper){
    ## Maybe duplicate to madeUpperGenes
    aq <- allDataDf$query
    aq <- aq[aq!=""]
    madeUpper <- c(madeUpper, tolower(unique(unlist(strsplit(aq, ",")))))
    # print(madeUpper)
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

  docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly, stem)
  ret@corpus <- docs
  if (onlyCorpus){
    return(docs)
  }
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
  
  if (onlyTDM){
    return(docs)
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

  # fetched[["rawfrequency"]] <- matSorted
  ret@TDM <- docs

  if (numWords > length(matSorted)){
        numWords <- length(matSorted)
  }
  ret@numWords <- numWords  
  
  if (plotType=="network"){
    matSorted <- matSorted[1:numWords]
    returnDf <- data.frame(word=names(matSorted),freq=matSorted)
    for (i in madeUpper) {
      # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
      returnDf[returnDf$word == i,"word"] <- toupper(i)
    }
    ret@freqDf <- returnDf
    
    freqWords <- names(matSorted)

    if (useUdpipe) {
      if (udpipeOnlyFreq & udpipeOnlyFreqN) {stop("Cannot specify both of these options")}
      if (udpipeOnlyFreq) {
        showNeighbors <- NULL
        showFreq <- freqWords
      }
      if (udpipeOnlyFreqN) {
        showNeighbors <- freqWords
        showFreq <- NULL
      }
      ret <- retUdpipeNet(ret=ret,texts=allDataDf,udmodel_english=udmodel_english,
          orgDb=orgDb, filterWords=filterWords, additionalRemove=additionalRemove,
          colorText=colorText,edgeLink=edgeLink,queryPlot=genePlot, layout=layout,
          pal=pal,showNeighbors=showNeighbors, showFreq=showFreq, nodePal=nodePal)
      return(ret)
    }


    # freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
    ## TODO: before or after?
    freqWordsDTM <- t(as.matrix(docs))
    row.names(freqWordsDTM) <- allDataDf$query
    if (tag) {
      if (!is.null(ret) & length(ret@pvpick)!=0){
        qqcat("Using previous pvclust results")
        pvcl <- ret@pvpick
      } else {
        if (tagWhole){
          pvc <- pvclust(as.matrix(dist(as.matrix(docs))), parallel=cl)
        } else {
          # pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
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
    }
    
    ## Check correlation
    if (bn) {
      qqcat("bn specified, R=@{R}\n")
      # To avoid computaitonal time, subset to numWords
      bnboot <- bnlearn::boot.strength(
        data.frame(freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]),
        algorithm = "hc", R=R)
      ret@strength <- bnboot
      av <- bnlearn::averaged.network(bnboot)
      avig <- bnlearn::as.igraph(av)
      el <- data.frame(as_edgelist(avig))
      colnames(el) <- c("from","to")
      mgd <- merge(el, bnboot, by=c("from","to"))
      colnames(mgd) <- c("from","to","weight","direction")
      coGraph <- graph_from_data_frame(mgd, directed=TRUE)
    } else {
      ## Check correlation
      if (onWholeDTM){
        corData <- cor(freqWordsDTM)
      } else {
        corData <- cor(freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords])
      }
      ret@corMat <- corData
      
      ## Set correlation below threshold to zero
      corData[corData<corThresh] <- 0
      coGraph <- graph.adjacency(corData, weighted=TRUE,
                                 mode="undirected", diag = FALSE)
    }
    ## before or after?
    coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
    V(coGraph)$Freq <- matSorted[V(coGraph)$name]
    ## Set pseudo freq as min value of freq
    # fre <- V(coGraph)$Freq
    # fre[is.na(fre)] <- min(fre, na.rm=TRUE)
    # V(coGraph)$Freq <- fre

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
    
    if (genePlot) {
      genemap <- c()
      for (rn in nodeName){
        tmp <- freqWordsDTM[ ,rn]
        for (nm in names(tmp[tmp!=0])){
          if (nm!=""){
            if (grepl(",",nm,fixed=TRUE)){
              for (nm2 in unlist(strsplit(nm, ","))){
                genemap <- rbind(genemap, c(rn, paste(nm2, "(Q)")))
              }
            } else {
              genemap <- rbind(genemap, c(rn, paste(nm, "(Q)")))
            }
          }
        }
      }
      # if (preserve) {
      #   retGenemap <- genemap
      #   gmnew <- NULL
      #   for (q in retGenemap[,1]) {
      #     if (q %in% names(pdic)){
      #       gmnew <- c(gmnew, pdic[q])
      #     } else {
      #       gmnew <- c(gmnew, q)
      #     }
      #   }
      #   retGenemap[,1] <- gmnew
      # } else {
      #   retGenemap <- genemap
      # }
      retGenemap <- genemap
      ret@geneMap <- retGenemap
      genemap <- simplify(igraph::graph_from_edgelist(genemap, directed = FALSE))
      coGraph <- igraph::union(coGraph, genemap)

      E(coGraph)$edgeColor <- E(coGraph)$weight
      tmpW <- E(coGraph)$weight
      if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
        corThreshGenePlot <- corThresh - 0.1}
      tmpW[is.na(tmpW)] <- corThreshGenePlot
      E(coGraph)$weight <- tmpW
    } else {
      E(coGraph)$edgeColor <- E(coGraph)$weight
    }
    
    if (tag) {
      netCol <- tolower(names(V(coGraph)))
      for (i in seq_along(pvcl$clusters)){
        for (j in pvcl$clusters[[i]])
          netCol[netCol==j] <- paste0("cluster",i)
      }
      netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
      V(coGraph)$tag <- netCol
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

    ## Main plot

    if (onlyGene) {
      qqcat("Subsetting to the gene symbol in orgDb\n")
      included <- names(V(coGraph))[tolower(names(V(coGraph))) %in% tolower(keys(orgDb, "SYMBOL"))]
      coGraph <- induced_subgraph(coGraph, included)
    }

    ret@igraph <- coGraph
    netPlot <- ggraph(coGraph, layout=layout)
    
    if (bn){
      if (edgeLink){
        if (edgeLabel){
          netPlot <- netPlot +
            geom_edge_link(
              aes(width=.data$weight,
                  color=.data$edgeColor,
                  label=round(.data$weight,3)),
              angle_calc = 'along',
              label_dodge = unit(2.5, 'mm'),
              arrow = arrow(length = unit(4, 'mm')), 
              start_cap = circle(3, 'mm'),
              end_cap = circle(3, 'mm'),
              alpha=0.5,
              show.legend = showLegend)
        } else {
          netPlot <- netPlot +
            geom_edge_link(aes(width=.data$weight, color=.data$edgeColor),
                           arrow = arrow(length = unit(4, 'mm')), 
                           start_cap = circle(3, 'mm'),
                           end_cap = circle(3, 'mm'),
                           alpha=0.5, show.legend = showLegend)
        }
      } else {
        if (edgeLabel){
          netPlot <- netPlot +
            geom_edge_diagonal(
              aes(width=.data$weight,
                  color=.data$edgeColor,
                  label=round(.data$weight,3)),
              angle_calc = 'along',
              label_dodge = unit(2.5, 'mm'),
              arrow = arrow(length = unit(4, 'mm')), 
              start_cap = circle(3, 'mm'),
              end_cap = circle(3, 'mm'),
              alpha=0.5,
              show.legend = showLegend)
        } else {
          netPlot <- netPlot +
            geom_edge_diagonal(aes(width=.data$weight, color=.data$edgeColor),
                               arrow = arrow(length = unit(4, 'mm')), 
                               start_cap = circle(3, 'mm'),
                               end_cap = circle(3, 'mm'),                                    
                               alpha=0.5, show.legend = showLegend)                
        }
      }
    } else {
      if (edgeLink){
        if (edgeLabel){
          netPlot <- netPlot +
            geom_edge_link(
              aes(width=.data$weight,
                  color=.data$edgeColor,
                  label=round(.data$weight,3)),
              angle_calc = 'along',
              label_dodge = unit(2.5, 'mm'),
              alpha=0.5,
              show.legend = showLegend)
        } else {
          netPlot <- netPlot +
            geom_edge_link(aes(width=.data$weight, color=.data$edgeColor),
                           alpha=0.5, show.legend = showLegend)
        }
      } else {
        if (edgeLabel){
          netPlot <- netPlot +
            geom_edge_diagonal(
              aes(width=.data$weight,
                  color=.data$edgeColor,
                  label=round(.data$weight,3)),
              angle_calc = 'along',
              label_dodge = unit(2.5, 'mm'),
              alpha=0.5,
              show.legend = showLegend)
        } else {
          netPlot <- netPlot +
            geom_edge_diagonal(aes(width=.data$weight, color=.data$edgeColor),
                               alpha=0.5, show.legend = showLegend)                
        }
      }
    }
    if (tag) {
      netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$tag),
                                           show.legend = showLegend)+
      scale_color_manual(values=nodePal)
    } else { 
      netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$Freq),
                                           show.legend = showLegend)+
        scale_color_gradient(low=pal[1],high=pal[2],
                             name = "Frequency")
    }
    if (colorText){
      if (tag) {
        netPlot <- netPlot + 
          geom_node_text(aes(label=.data$name, size=.data$Freq, color=.data$tag),
                         check_overlap=TRUE, repel=TRUE,# size = labelSize,
                         bg.color = "white", segment.color="black",
                         bg.r = .15, show.legend=showLegend)
      } else {
        netPlot <- netPlot + 
          geom_node_text(aes(label=.data$name, size=.data$Freq, color=.data$Freq),
                         check_overlap=TRUE, repel=TRUE,# size = labelSize,
                         bg.color = "white", segment.color="black",
                         bg.r = .15, show.legend=showLegend)
      }
    } else {
      netPlot <- netPlot + 
        geom_node_text(aes(label=.data$name, size=.data$Freq),
                       check_overlap=TRUE, repel=TRUE,# size = labelSize,
                       color = "black",
                       bg.color = "white", segment.color="black",
                       bg.r = .15, show.legend=showLegend) 
    }
    netPlot <- netPlot+
      scale_size(range=scaleRange, name="Frequency")+
      scale_edge_width(range=c(1,3), name = "Correlation")+
      scale_edge_color_gradient(low=pal[1],high=pal[2],na.value=naEdgeColor,
                                name = "Correlation")+
      theme_graph()
    ret@net <- netPlot
  } else {
    ## WC part
    matSorted <- matSorted[1:numWords]
    returnDf <- data.frame(word = names(matSorted),freq=matSorted)
    freqWords <- names(matSorted)
    freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
    
    if (tag) {
      if (!is.null(redo) & length(redo@pvpick)!=0) {
        qqcat("Using previous pvclust results")
        pvcl <- redo@pvpick
      } else {
        if (tagWhole){
          pvc <- pvclust(as.matrix(dist(as.matrix(docs))), parallel=cl)
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
                                         ordered.colors = TRUE)))
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
