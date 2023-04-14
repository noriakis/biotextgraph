#'
#' retUdpipeNet
#'
#' return network using udpipe dependencies
#'
#' @param ret biotext object
#' @param texts data.frame
#' @param showNeighbors only show nodes with the neighbors of these words
#' @noRd
#' 
retUdpipeNet <- function(ret,texts,udmodel_english,orgDb,
                         filterWords,additionalRemove,colorText,
                         edgeLink, queryPlot, layout, pal,
                         showNeighbors, showFreq, nodePal,addNet=NULL) {

  ret@model <- "udpipe"
  ## Frequency
  freq <- udpipe::document_term_frequencies(texts$text) |>
    dplyr::group_by(.data$term) |> dplyr::summarise(sum=sum(freq))
  vfreq <- freq$sum
  names(vfreq) <- freq$term

  ## Queries
  if (!"query" %in% colnames(texts)) {
    revID <- AnnotationDbi::select(orgDb,
                                   keys = as.character(texts$Gene_ID), 
                                   columns = c("SYMBOL"),
                                   keytype = "ENTREZID")$SYMBOL
    texts$query <- revID
  }

  ## Annotate using udpipe
  alledges <- NULL
  allwordatt <- NULL
  allqueries <- NULL
  
  ## Not include these in graph
  notInc <- c("ADP","AUX","DET","PUNCT",
    "PART","PRON","SCONJ","CCONJ")

  for (gid in texts$ID) { ## To match with query

    edges <- NULL
    ges <- NULL
    tmp <- subset(texts, texts$ID==gid)
    gsym <- tmp$query
    if (grepl(",",gsym)) {
      gsym <- unlist(strsplit(gsym,","))
    }
    tmpm <- udpipe::udpipe_annotate(udmodel_english, tmp$text)
    x <- data.frame(tmpm)
    
    ## Use token, not lemma
    wordatt <- x$upos
    names(wordatt) <- x$token
    
    for (sent in x$sentence_id) {
      one <- subset(x, x$sentence_id==sent)
      for (tkid in one$token_id) {
        sampleTkID <- subset(one, one$token_id==tkid)
        if (!sampleTkID$upos %in% notInc){
          e1 <- sampleTkID$token
          e2n <- sampleTkID$head_token_id
          e2 <- subset(one, one$token_id==e2n)$token
          if (length(e1)!=0 & length(e2)!=0) {
            edges <- rbind(edges,c(e1, e2))
            for (cgsym in gsym) {
              ges <- rbind(ges, c(e1, cgsym))
              ges <- rbind(ges, c(e2, cgsym))
            }
          }
        }
      }
      
    }
    alledges <- rbind(alledges,edges)
    allqueries <- rbind(allqueries, ges)
    
    allwordatt <- c(allwordatt, wordatt)
  }

  allwordatt <- tapply(allwordatt,
    names(allwordatt),
    function(x) {
      if(length(unique(x))!=1){
        paste(unique(x), collapse=",")
      } else {
        unique(x)
      }})

  geGraph <- simplify(igraph::graph_from_data_frame(allqueries, directed=FALSE))
  udpGraph <- simplify(igraph::graph_from_data_frame(alledges,directed=FALSE))

  if (queryPlot) {
    udpGraph <- igraph::union(udpGraph, geGraph)
  }

  nodeN <- NULL  
  if (!is.null(addNet)) {
    for (netName in names(addNet)) {
        tmpAdd <- addNet[[netName]]
        tmpNN <- names(V(tmpAdd))
        tmpNN <- tmpNN[!tmpNN %in% names(nodeN)]

        newNN <- rep(netName, length(tmpNN))
        names(newNN) <- tmpNN
        nodeN <- c(nodeN, newNN)

        udpGraph <- igraph::union(udpGraph, tmpAdd)
    }
  }


  cat <- NULL
  fre <- NULL
  for ( i in names(V(udpGraph))) {
    fre <- c(fre, vfreq[i])
    if (i %in% names(allwordatt)) {
      cat <- c(cat, allwordatt[i])
    } else if (i %in% texts$query ){
      cat <- c(cat, "Query")
    } else if (i %in% names(nodeN)) {
      cat <- c(cat, nodeN[i])
    } else {
      cat <- c(cat, "Others")
    }
  }

  ## Set pseudo freq as min value of freq
  fre[is.na(fre)] <- min(fre, na.rm=TRUE)
  
  V(udpGraph)$cat <- cat
  V(udpGraph)$freq <- fre

  udpGraph <- induced_subgraph(udpGraph,
    !tolower(names(V(udpGraph))) %in% tolower(c(filterWords,
                                          additionalRemove)))

  if (!is.null(showFreq)) {
    qqcat("Subsetting to the specified @{length(showFreq)} words\n")
    showFreq <- showFreq[!tolower(showFreq) %in% tolower(unique(allqueries[,2]))]
    udpGraph <- induced_subgraph(udpGraph,
      tolower(names(V(udpGraph))) %in% tolower(showFreq))
  }

  if (!is.null(showNeighbors)){
    qqcat("Subsetting to the neighbors of specified @{length(showNeighbors)} words\n")
    inc <- NULL
    ## Exclude queries
    showNeighbors <- showNeighbors[!tolower(showNeighbors) %in% tolower(unique(allqueries[,2]))]
    for (nn in showNeighbors) {
      if (tolower(nn) %in% tolower(names(V(udpGraph))))
        nn2 <- names(V(udpGraph))[tolower(names(V(udpGraph)))==tolower(nn)]
      inc <- c(inc, 
        names(neighbors(udpGraph, nn2)))
    }
    udpGraph <- induced_subgraph(udpGraph,
      tolower(names(V(udpGraph))) %in% tolower(inc))
  }
  
  ret@igraph <- udpGraph

  net <- ggraph(udpGraph, layout=layout)
  if (edgeLink) {
    net <- net + geom_edge_link(color="grey")
  } else {
    net <- net + geom_edge_diagonal(color="grey")
  }
  

  if (colorText) {
    net <- net + geom_node_point(aes(color=freq, size=freq))
    net <- net + geom_node_text(aes(label=.data$name, color=freq, size=freq),
                                check_overlap=TRUE, repel=TRUE,
                                bg.color = "white", segment.color="black",
                                bg.r = .15, show.legend=FALSE)+
                scale_color_gradient(low=pal[1],high=pal[2], name="Frequency")
  } else {
    net <- net + geom_node_point(aes(color=cat, size=freq))
    net <- net + geom_node_text(aes(label=.data$name, size=freq),
                   check_overlap=TRUE, repel=TRUE,
                   bg.color = "white", segment.color="black",
                   bg.r = .15, show.legend=FALSE)+
    scale_color_manual(values=nodePal, name="Category")
  }
  net <- net + scale_size(range=c(3,9), name="Frequency")+ theme_graph()
  ret@net <- net
  return(ret)
}