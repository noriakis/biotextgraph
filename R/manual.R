#' manual
#' 
#' Produce networks using manual input
#' 
#' @param df df
#' @param madeUpper make the words uppercase in resulting plot
#' @param pal palette for color gradient in correlation network
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param corThresh the correlation threshold
#' @param autoThresh automatically determine the threshold value to show `numWords`
#' @param cooccurrence default to FALSE, if TRUE, use cooccurrence instead of correlation
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @param colorText color text label based on frequency in correlation network
#' @param ngram default to NA (1)
#' @param additionalRemove specific words to be excluded
#' @param tag perform pvclust on words and colorlize them in wordcloud or network
#' argument of "cor" or "tdm". Default to "none", which performs no tagging.
#' If wordcloud, tagging will be performed on TDM.
#' @param tagWhole tag based on whole data or subset
#' @param useFil filter based on "GS_TfIdf" (whole gene summary tf-idf)
#'  or "BSDB_TfIdf" (whole bugsigdb tf-idf)
#' @param filNum specify filter tfidf
#' @param filType "above" or "below"
#' @param tfidf use TfIdf when making TDM
#' @param pvclAlpha alpha for pvpick()
#' @param onlyCorpus return only corpus
#' @param onlyTDM return only TDM, if quanteda, return DFM
#' @param preserve try to preserve the original characters
#' @param numOnly delete number only
#' @param cl cluster to pass to pvclust (snow::makeCluster(n))
#' @param bn perform bootstrap-based Bayesian network inference 
#' instead of correlation using bnlearn
#' @param R how many bootstrap when bn is stated
#' @param onWholeDTM calculate correlation network
#'                   on whole dataset or top-words specified by numWords
#' @param stem whether to use stemming
#' 
#' @param colorize color the nodes and texts based on their category
#' @param catColors named vector showing colors for each category
#' (e.g. `query`, `Words` and specified column names in data.frame)
#' @param discreteColorWord colorize words by "Words" category, not frequency.
#' @param tagPalette palette for each tag when tag is TRUE, default to NULL
#' (e.g. `cluster_1`, `cluster_2`, ...)
#' @param queryColor color for query in the network plot
#' if tagPalette or catColors is NULL, determined palettes with 
#' queryColor replaced with this color will be used.
#' 
#' @param takeMax when summarizing term-document matrix, take max.
#' Otherwise take sum.
#' @param argList parameters to pass to wordcloud()
#' @param normalize sum normalize the term frequency document-wise
#' @param takeMean take mean values for each term in term-document matrix
#' @param queryPlot plot query
#' @param useQuanteda use quanteda functions to generate
#' @param quantedaArgs list of arguments to be passed to tokens()
#' @param naEdgeColor edge color values for NA
#' @param collapse default to FALSE, collapse all the sentences
#' @param useUdpipe use udpipe to make a network
#' @param udpipeModel udpipe model file name
#' @param useggwordcloud default to TRUE
#' @param wcScale scaling size for ggwordcloud
#' @param fontFamily font family to use, default to "sans"
#' @param useSeed seed
#' @param scaleFreq scale the frequency
#' @param addFreqToNonWords add pseudo-frequency corresponding to minimum
#' frequency of the words to nodes other than words
#' @examples
#' ret <- refseq("DDX41", plotType="wc")
#' manual(ret@rawText$Gene_summary, plotType="wc")
#' @export
#' @return list of data frame and ggplot2 object
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
manual <- function(df, madeUpper=NULL,
   useFil=NA, filType="above", cooccurrence=FALSE,
   filNum=0, useQuanteda=FALSE, quantedaArgs=list(),
   pvclAlpha=0.95, numOnly=TRUE, tfidf=FALSE, cl=FALSE,
   pal=c("blue","red"), numWords=30, scaleRange=c(5,10), scaleFreq=NULL,
   showLegend=FALSE, plotType="network", colorText=FALSE,
   corThresh=0.2, layout="nicely", tag="none", tagWhole=FALSE,
   onlyCorpus=FALSE, onlyTDM=FALSE, bn=FALSE, R=20, queryColor="grey",
   edgeLabel=FALSE, edgeLink=TRUE, ngram=1, colorize=FALSE,
   tagPalette=NULL, preserve=TRUE, takeMax=FALSE, catColors=NULL,
   deleteZeroDeg=TRUE, additionalRemove=NA, naEdgeColor="grey50",
   normalize=FALSE, takeMean=FALSE, queryPlot=FALSE, collapse=FALSE,
   onWholeDTM=FALSE, stem=FALSE, argList=list(), useUdpipe=FALSE,
   discreteColorWord=FALSE, autoThresh=TRUE,
   useggwordcloud=TRUE, wcScale=10, fontFamily="sans",
   addFreqToNonWords=FALSE,
   udpipeModel="english-ewt-ud-2.5-191206.udpipe", useSeed=42)
{

	if (!tag %in% c("none","tdm","cor")) {
		stop("tag should be none, tdm, or cor.")
	}

    if (useUdpipe) {
      qqcat("Using udpipe mode\n")
      plotType="network"
      udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)
    }
    if (!is.data.frame(df)) {
      if (is.vector(df)) {
        df <- data.frame(df) |> `colnames<-`(c("text"))
      }
    }
    if (queryPlot) {
      if (!"query" %in% colnames(df)) {stop("There is no query specified")}
    }

    ret <- new("biotext")
    ret@type <- paste0("manual")
    ret@rawText <- df

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
    if (useQuanteda) {
      preserve <- FALSE
      ret <- returnQuanteda(ret, quantedaArgs,numWords,ngram,
                           filterWords,additionalRemove, tfidf, collapse)
      if (onlyCorpus) {return(ret@corpusQuanteda)}
      if (onlyTDM) {return(ret@dfm)}
      docs <- t(as.matrix(ret@dfm))
      mat <- t(as.matrix(ret@dfm))
    } else {
      if (collapse) {
        docs <- VCorpus(VectorSource(paste(df$text, collapse=" ")))
      } else {
        docs <- VCorpus(VectorSource(df$text))
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
      ret@TDM <- docs
      if (onlyTDM){
        return(docs)
      }
      mat <- as.matrix(docs)
    }

    if (useUdpipe) {
        if (!"query" %in% colnames(df)) {
          df$query <- seq_len(nrow(df))
          queryPlot <- FALSE
        }
        df$ID <- 1:nrow(df)#df$query
        ret <- retUdpipeNet(ret=ret, texts=df,udmodel_english=udmodel_english,
            orgDb=NULL, filterWords=filterWords, additionalRemove=additionalRemove,
            colorText=colorText,edgeLink=edgeLink,queryPlot=queryPlot, layout=layout,
            pal=pal, showNeighbors=NULL, showFreq=NULL, nodePal=tagPalette)
        return(ret)
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
  
    if (plotType=="network"){
      matSorted <- matSorted[1:numWords]
      returnDf <- data.frame(word=names(matSorted),freq=matSorted)
      for (i in madeUpper) {
        returnDf[returnDf$word == i,"word"] <- toupper(i)
      }
      ret@freqDf <- returnDf
    
      freqWords <- names(matSorted)
      DTM <- t(as.matrix(docs))
      row.names(DTM) <- df$query
      
      if (tag=="tdm") {
        if (!is.null(ret) & length(ret@pvpick)!=0){
          qqcat("Using previous pvclust results")
          pvcl <- ret@pvpick
        } else {
          if (tagWhole){
            pvc <- pvclust(as.matrix(dist(as.matrix(docs))), parallel=cl)
          } else {
            pvc <- pvclust(as.matrix(dist(
              t(
                DTM[, colnames(DTM) %in% freqWords]
              )
            )), parallel=cl)
          }
          pvcl <- pvpick(pvc, alpha=pvclAlpha)
          ret@pvclust <- pvc
          ret@pvpick <- pvcl
        }
      }
    
      matrixs <- obtainMatrix(ret, bn, R, DTM, freqWords,
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
      V(coGraph)$type <- "Words"
    
      if (deleteZeroDeg){
        coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
      }
    
      nodeName <- V(coGraph)$name
      
      incCols <- colnames(df)
      incCols <- incCols[!incCols %in% c("query","text")]
      
      if (length(incCols)!=0) {
        qqcat("Including columns @{paste(incCols, collapse=' and ')} to link with query\n")

        for (ic in incCols) {
          querymap <- df[,c(ic, "query")]
          vtx <- data.frame(cbind(c(querymap[,2], querymap[,1]),
            c(rep("Queries",length(querymap[,2])),
              rep(ic,length(querymap[,1]))))) |> 
          `colnames<-`(c("name","type"))
          vtx <- vtx[!duplicated(vtx),]
          vtx <- vtx |> `rownames<-`(1:nrow(vtx))
          eds <- data.frame(querymap)
          words <- vtx |> subset(vtx$type==ic)
          queriesDf <- vtx |> subset(vtx$type=="Queries")
          row.names(words)[which(words$name %in% eds[,1])]
          row.names(queriesDf)[which(queriesDf$name %in% eds[,2])]
          eds[,1] <- sapply(eds[,1], function(x) {
            as.integer(row.names(words)[which(words$name %in% x)])
          })
          eds[,2] <- sapply(eds[,2], function(x) {
            as.integer(row.names(queriesDf)[which(queriesDf$name %in% x)])
          })
          qmap <- tbl_graph(nodes=vtx,edges=eds,directed=FALSE)
          coGraph <- graph_join(as_tbl_graph(coGraph),
                                qmap)
        }
      }


      ## `nodeN` holds the named vector
      ## of which category each node in the network
      ## belongs to.
      nodeN <- NULL
      for (coln in c(incCols, "query")) {
        tmpn <- df[[coln]]
        tmpnn <- rep(coln, length(tmpn))
        names(tmpnn) <- tmpn
        nodeN <- c(nodeN, tmpnn)
      }


      if (queryPlot) {
        genemap <- c()
        for (rn in nodeName){
          tmp <- DTM[ ,rn]
          for (nm in names(tmp[tmp!=0])){
            if (nm!=""){
              if (grepl(",",nm,fixed=TRUE)){
                for (nm2 in unlist(strsplit(nm, ","))){
                  genemap <- rbind(genemap, c(rn, paste(nm2)))
                }
              } else {
                genemap <- rbind(genemap, c(rn, paste(nm)))
              }
            }
          }
        }
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
        qmap <- tbl_graph(nodes=vtx,edges=eds,directed=FALSE)
        coGraph <- graph_join(as_tbl_graph(coGraph),
                              qmap)
        # genemap <- simplify(igraph::graph_from_edgelist(genemap, directed = FALSE))
        # coGraph <- igraph::union(coGraph, genemap)
      }

      E(coGraph)$edgeColor <- E(coGraph)$weight
      tmpW <- E(coGraph)$weight
      if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
        corThreshGenePlot <- corThresh - 0.1}
      tmpW[is.na(tmpW)] <- corThreshGenePlot
      E(coGraph)$weight <- tmpW

      if (tag!="none") {
        netCol <- tolower(names(V(coGraph)))
        for (i in seq_along(pvcl$clusters)){
          for (j in pvcl$clusters[[i]])
            netCol[netCol==j] <- paste0("cluster",i)
        }
        netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
        V(coGraph)$tag <- netCol

        addC <- V(coGraph)$tag
        for (nn in seq_along(names(V(coGraph)))) {
            if (names(V(coGraph))[nn] %in% names(nodeN)) {
                if (addC[nn]=="not_assigned") {
                  addC[nn] <- nodeN[names(V(coGraph))[nn]]
                } else {
                  addC[nn] <- paste0(addC[nn], ";", nodeN[names(V(coGraph))[nn]])
                }
            } else {
                next
            }
        }
        V(coGraph)$tag <- addC
      }

      if (addFreqToNonWords) {
        fre <- V(coGraph)$Freq
        fre[is.na(fre)] <- min(fre, na.rm=TRUE)
        V(coGraph)$Freq <- fre
      }

      if (colorize) {
        ## Set pseudo freq based on min value of freq
        fre <- V(coGraph)$Freq
        fre[is.na(fre)] <- min(fre, na.rm=TRUE)
        V(coGraph)$Freq <- fre

        if (tag!="none") {
        	qqcat("Overriding tagged information by pvclust by colorize option\n")}
        if (!is.null(nodeN)) {
            addC <- NULL
            for (nn in seq_along(names(V(coGraph)))) {
                if (names(V(coGraph))[nn] %in% names(nodeN)) {
                    addC[nn] <- nodeN[names(V(coGraph))[nn]]
                } else {
                    addC[nn] <- "Words"
                }
            }
            V(coGraph)$tag <- addC
        }
      }

      if (!is.null(nodeN)) {
          nodeCat <- NULL
          for (nn in seq_along(names(V(coGraph)))) {
              if (names(V(coGraph))[nn] %in% names(nodeN)) {
                  nodeCat[nn] <- nodeN[names(V(coGraph))[nn]]
              } else {
                  nodeCat[nn] <- "Words"
              }
          }
      } else {
          nodeCat <- rep("Words",length(V(coGraph)))
      }
      V(coGraph)$nodeCat <- nodeCat


      if (preserve) {
        nodeDf <- coGraph |> as_tbl_graph() |> activate("nodes") |> data.frame()
        V(coGraph)$name <- apply(nodeDf,
              1,
              function(x) {ifelse(x["type"]=="Words",
              	ifelse(is.na(pdic[x["name"]]), x["name"], pdic[x["name"]]),
                x["name"])})
      }

      nodeName <- V(coGraph)$name
      for (i in madeUpper) {
        nodeName[nodeName == i] <- toupper(i)
      }
      V(coGraph)$name <- nodeName

    
      ## Main plot
      if (!is.tbl_graph(coGraph)) {
          ret@igraph <- coGraph
      } else {
          ret@igraph <- as.igraph(coGraph)
      }
      netPlot <- ggraph(coGraph, layout=layout)
      netPlot <- appendEdges(netPlot, bn, edgeLink,
            edgeLabel, showLegend, fontFamily)

      ## Define colors
      if (tag!="none") {
        if (is.null(tagPalette)) {
          cols <- V(coGraph)$tag |> unique()
          tagPalette <- RColorBrewer::brewer.pal(8, "Dark2")
          tagPalette <- colorRampPalette(tagPalette)(length(unique(V(coGraph)$tag)))
          names(tagPalette) <- cols
          tagPalette["query"] <- queryColor
        }
      }

      if (is.null(catColors)) {
        catColors <- RColorBrewer::brewer.pal(length(unique(V(coGraph)$nodeCat)), "Dark2")
        names(catColors) <- unique(V(coGraph)$nodeCat)
        catColors["query"] <- queryColor
      }

      tagColors <- tagPalette
      netPlot <- appendNodesAndTexts(netPlot,tag,colorize,tagPalette,
                          showLegend,catColors,pal,fontFamily,colorText,scaleRange,
                          useSeed, ret, tagColors, discreteColorWord=discreteColorWord)
      netPlot <- netPlot+
        scale_size(range=scaleRange, name="Frequency")+
        scale_edge_width(range=c(1,3), name = "Correlation")+
        scale_edge_color_gradient(low=pal[1],high=pal[2],
                                  name = "Correlation",
                                  na.value=naEdgeColor)+
        theme_graph()
      ret@net <- netPlot
  } else {
    ## WC part
    matSorted <- matSorted[1:numWords]
    returnDf <- data.frame(word = names(matSorted),freq=matSorted)
    freqWords <- names(matSorted)
    freqWordsDTM <- t(as.matrix(docs[row.names(docs) %in% freqWords, ]))
    
    if (tag!="none") {
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
