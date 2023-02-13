#' compareWordNet
#' 
#' compare two gene clusters based on words 
#' 
#' @param listOfNets list consisting results of wc* functions (plotType="network")
#' @param titles title to be shown on plot
#' @param layout layout
#' @param hull show category by hull
#' @param size node size, freq, overlap or numeric number
#' @param tag show tag on plot
#' @param tagLevel words that tagged in this many networks will be included
#' @param conc concavity parameter
#' @param edgeLink whether to use link or diagonal
#' @param freqMean how to concatenate frequency (TRUE: mean, FALSE: sum)
#' @param scaleRange point size scaling
#' @param returnNet return only the network (ig)
#' @param colPal color palette to be used in RColorBrewer
#' @param colNum color number to be used in plot
#' @param colorText whether to color text based on category
#' @param ovlThresh show text with this number of overlap between graphs
#' 
#' @export
#' @examples
#' net1 <- wcGeneSummary(c("DDX41","IRF3"), plotType="network")
#' net2 <- wcGeneSummary(c("DDX41","PNKP"), plotType="network")
#' compare <- compareWordNet(list(net1, net2))
#' @return plot comparing gene clusters
#' @importFrom grDevices colorRampPalette
compareWordNet <- function(listOfNets, titles=NULL,
                           layout="nicely", hull=FALSE, size="freq", conc=1,
                           tag=FALSE, tagLevel=1, edgeLink=TRUE,
                           freqMean=FALSE, scaleRange=c(5,10), ovlThresh=0,
                           returnNet=FALSE, colPal="Pastel1", colNum=20,
                           colorText=FALSE) {
  listOfIGs <- list()
  listOfNodes <- list()

  if (is.null(titles)){
    titles <- c()
    for (e in seq_along(listOfNets)){
      titles <- c(titles, paste0("title",e))
    }
  }

  for (e in seq_along(listOfNets)) {
      listOfIGs[[e]] <- listOfNets[[e]]@igraph
      listOfNodes[[e]] <- names(V(listOfNets[[e]]@igraph))
  }

  names(listOfIGs) <- titles

  if (tag) {
    for (g in names(listOfIGs)) {
      tmpg <- listOfIGs[[g]]
      V(tmpg)$tag <- paste0(g,"_",V(tmpg)$tag)
      listOfIGs[[g]] <- tmpg
    }
  }

  commonNodes <- Reduce(intersect, listOfNodes)
  uig <- simplify(Reduce(igraph::union, listOfIGs))
  nodeAttr <- names(get.vertex.attribute(uig))
  tagName <- nodeAttr[grepl("tag", nodeAttr)]
  freqName <- nodeAttr[grepl("Freq", nodeAttr)]

  ## concatenate tags
  if (tag) {
      tags <- c()
      for (tn in tagName){
          tags <- cbind(tags, get.vertex.attribute(uig, tn))
      }
      contag <- apply(tags, 1, function(x){
          septag <- x[!is.na(x) & x!="not_assigned"]
          
        if (length(septag)>=tagLevel) {
            return(paste(septag, collapse="_"))
        } else {
            return(NA)
        }
      } )
      V(uig)$tag <- contag
  }


  frqs <- c()
  for (fq in freqName){
    tmpfrq <- get.vertex.attribute(uig, fq)
    frqs <- cbind(frqs, tmpfrq)
  }
  if (freqMean){
    V(uig)$Freqs <- apply(frqs, 1, function(x) mean(x, na.rm=TRUE))
  } else {
    V(uig)$Freqs <- apply(frqs, 1, function(x) sum(x, na.rm=TRUE))
  }


  col <- NULL
  ovl <- NULL
  for (node in names(V(uig))){
      if (node %in% commonNodes) {
          col <- c(col, "Common")
          howm <- length(listOfNets)
      } else {
          tmpcol <- c()
          for (e in seq_along(listOfIGs)){
              if (node %in% names(V(listOfIGs[[e]]))) {
                  if (is.null(titles)){
                    tmpcol <- c(tmpcol, e)
                  } else {
                    tmpcol <- c(tmpcol, titles[e])
                  }
              }
          }
          howm <- length(tmpcol)
          tmpcol <- paste(tmpcol, collapse="_")
          col <- c(col, tmpcol)
      }
      ovl <- c(ovl, howm)
  }
  V(uig)$col <- col
  V(uig)$ovl <- ovl
  
  ## concatenate tags
  # if (tag) {
  #   contag <- c()
  #   for (i in seq_along(names(V(uig)))){
  #     t1 <- str_replace(V(uig)$tag_1[i], "cluster", "")
  #     t2 <- str_replace(V(uig)$tag_2[i], "cluster", "")
  #     if (is.na(t1) | is.na(t2)) {
  #       contag <- c(contag, NA)
  #     } else if (t1=="not_assigned"|t2=="not_assigned"){
  #       contag <- c(contag, NA)
  #     } else {
  #       contag <- c(contag, paste0(t1,t2))
  #     }
  #   }

  #   V(uig)$tag <- contag
  # }
  # print(contag)
  
  # col <- c()
  # for (node in names(V(uig))){
  #   if (node %in% commonNodes) {
  #     col <- c(col, "Common")
  #   } else if (node %in% names(V(wc1$ig))){
  #     col <- c(col, titles[1])
  #   } else {
  #     col <- c(col, titles[2])
  #   }
  # }
  # V(uig)$col <- col

  if (size=="freq") {
    V(uig)$size <- V(uig)$Freqs
  } else if (size=="ovl") {
    V(uig)$size <- V(uig)$ovl
  } else {
    V(uig)$size <- rep(size, length(names(V(uig))))
  }

  if (returnNet){
    return(uig)
  }

  catNum <- length(unique(V(uig)$col))
  ## You can change it later
  cs <- RColorBrewer::brewer.pal(catNum, colPal)
  if (length(cs)<colNum) {
    cs <- colorRampPalette(cs)(colNum)
  }


  comNet <- ggraph(uig, layout=layout)

  if (edgeLink){
    comNet <- comNet+
    geom_edge_link(color="grey")
  } else {
    comNet <- comNet+
    geom_edge_diagonal(color="grey")
  }
  
  if (tag){
    ## Hull sometimes make groups inconsistent
    comNet <- comNet + 
              ggforce::geom_mark_hull(
                aes(.data$x,
                  .data$y,
                  group=.data$tag,
                  fill=.data$tag,
                filter=!is.na(.data$tag)),
                concavity=conc,
                alpha=0.25,
                na.rm=FALSE,
                show.legend=TRUE,
                inherit.aes=TRUE)
    ## TODO:
    ## specifying label produces an error,
    ## thus show.legend=TRUE is specified
    comNet <- comNet + 
      geom_node_point(aes(color=.data$col,
        size=.data$size))+
        scale_color_discrete(name="Group")
    # comNet <- comNet + geom_mark_hull(
    #   aes(comNet$data$x,
    #       comNet$data$y,
    #       group=tag,
    #       fill=tag,
    #       label=tag,
    #       filter = !is.na(tag)),
    #   concavity = conc,
    #   # expand = unit(2, "mm"),
    #   alpha = 0.25,
    #   na.rm = TRUE,
    #   # label.fill="transparent",
    #   show.legend = FALSE
    # )
  } else {
  
    if (hull) {
      comNet <- comNet + 
        geom_node_point(aes(color=col, size=size))+
        scale_color_discrete(name="Group")
      comNet <- comNet + 
        ggforce::geom_mark_hull(
        aes(comNet$data$x,
            comNet$data$y,
            group = col,
            label=col, fill=col,
            filter = !is.na(col)),
        concavity = conc,
        # expand = unit(2, "mm"),
        alpha = 0.25,
        na.rm = TRUE,
        # label.fill="transparent",
        show.legend = TRUE
      )
    } else {
      comNet <- comNet + 
        geom_node_point(aes(color=col, size=size)) +
        scale_color_manual(name="Group",values=cs)
    }
  }
  if (colorText) {
    comNet <- comNet +
    geom_node_text(
      aes(label=comNet$data$name, color=col, size=size, filter=!is.na(tag) & ovl > ovlThresh),
      check_overlap=TRUE, repel=TRUE,# size = labelSize,
      bg.color = "white", segment.color="black",
      bg.r = .15, show.legend=FALSE)+
    scale_size(range=scaleRange, name="Frequency")+
    theme_graph()
  } else {
    comNet <- comNet +
    geom_node_text(
      aes(label=comNet$data$name, size=size, filter=!is.na(tag) & ovl > ovlThresh),
      check_overlap=TRUE, repel=TRUE,# size = labelSize,
      bg.color = "white", segment.color="black",
      bg.r = .15, show.legend=FALSE)+
    scale_size(range=scaleRange, name="Frequency")+
    theme_graph()
  }
  comNet
}


#' plotDynamic
#' 
#' list network of words using graphlayouts::layout_as_dynamic
#' 
#' @param listOfNets list consisting results of wc* functions (plotType="network")
#' @param concat "union" or "intersection"
#' @param alpha pass to layout_as_dynamic
#' @param titles title to be shown on plot
#' @param tag show tag on plot
#' 
#' @export
#' @examples
#' library(igraph)
#' wc1 <- wcGeneSummary(c("DDX41","IRF3"), plotType="network")
#' wc2 <- wcGeneSummary(c("DDX41","PNKP"), plotType="network")
#' compare <- plotDynamic(list(wc1, wc2))
#' @return plot comparing gene clusters
#' @importFrom graphlayouts layout_as_dynamic
plotDynamic <- function(listOfNets,concat="union",alpha=0.5,titles=NULL,tag=FALSE){

  if (is.null(titles)){
    titles <- c()
    for (e in seq_along(listOfNets)){
      titles <- c(titles, paste0("title",e))
    }
  }

  if (concat=="union"){
    allNodes <- c()
    for (n in listOfNets){
        allNodes <- c(allNodes, V(n@igraph)$name)
    }
    igList <- list()
    for (e in seq_along(listOfNets)){
        tmpadd <- setdiff(allNodes, names(V(listOfNets[[e]]@igraph)))
        igList[[e]] <- add_vertices(listOfNets[[e]]@igraph,
                                      length(tmpadd),
                                      attr=list(name=tmpadd))
    }
  } else {
    uniqNodeNames <- list()
    for (e in seq_along(listOfNets)){
        uniqNodeNames[[e]] <- names(V(listOfNets[[e]]@igraph))
    }
    uniqNodeNames <- Reduce(intersect, uniqNodeNames)
    igList <- list()
    for (e in seq_along(listOfNets)){
        igList[[e]] <- induced_subgraph(listOfNets[[e]]@igraph,
          names(V(listOfNets[[e]]@igraph)) %in% uniqNodeNames)
    }
  }

  xy <- layout_as_dynamic(igList, alpha = alpha)
  pList <- vector("list",length(igList))

  for(i in seq_along(igList)){
    pList[[i]] <- ggraph(igList[[i]],layout="manual",x=xy[[i]][,1],y=xy[[i]][,2])+
        geom_edge_link0(edge_width=0.6,edge_colour="grey66")
    if (tag){
      pList[[i]] <- pList[[i]] + geom_node_point(aes(size=.data$Freq,fill=.data$tag),
        shape=21, show.legend=FALSE)
    } else {
      pList[[i]] <- pList[[i]] + geom_node_point(aes(size=.data$Freq,fill=.data$Freq),
        shape=21, show.legend=FALSE)+
                    scale_fill_gradient(low="blue",high="red",
                           name = "Frequency")
    }
      pList[[i]] <- pList[[i]] +
        geom_node_text(aes(label=.data$name),repel = TRUE,bg.color="white")+
        theme_graph()+
        theme(legend.position="bottom")+
        labs(title=titles[i])
  }
  Reduce("+",pList)
}