#' compareWordNet
#' 
#' @description Compare multiple networks based on words 
#' @details The function accepts list (named) of biotext object, and 
#' plot the merged network highlighting the intersection of the network
#' and identified clusters.
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
#' @param colorText whether to color text based on category
#' @param ovlThresh show text with this number of overlap between graphs
#' @param community compare based on community (igraph), override tag 
#' @param returnClass return biotext class object, default to TRUE
#' 
#' @export
#' @examples
#' net1 <- refseq(c("DDX41","IRF3","ERCC1","ERCC2","XRCC1"), plotType="network")
#' net2 <- refseq(c("DDX41","PNKP","ERCC3","IRF3","COPA"), plotType="network")
#' compare <- compareWordNet(list(net1, net2))
#' @return plot comparing gene clusters
#' @importFrom grDevices colorRampPalette
#' @importFrom ggforce geom_mark_hull
compareWordNet <- function(listOfNets, titles=NULL,
                           layout="nicely", hull=FALSE, size="freq", conc=1,
                           tag=FALSE, tagLevel=1, edgeLink=TRUE,
                           freqMean=FALSE, scaleRange=c(5,10), ovlThresh=0,
                           returnNet=FALSE, colPal="Pastel1",
                           colorText=FALSE, community=FALSE, returnClass=TRUE) {
  ret <- new("biotext")
  ret@type <- "combine"
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
  listOfTGs <- list()
  for (e in seq_along(listOfIGs)) {
    listOfTGs[[e]] <- as_tbl_graph(listOfIGs[[e]])
  }
  uig <- Reduce(function(x,y) graph_join(x,y, by="name"), listOfTGs)
  #simplify(Reduce(igraph::union, listOfIGs))
  nodeAttr <- names(get.vertex.attribute(uig))
  tagName <- nodeAttr[grepl("tag", nodeAttr)]
  communityName <- nodeAttr[grepl("community", nodeAttr)]
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

  if (community) {
      tags <- c()
      for (tn in communityName){
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
  
  ret@igraphRaw <- as.igraph(uig)

  catNum <- length(unique(V(uig)$col))
  ## You can change it later
  cs <- RColorBrewer::brewer.pal(catNum, colPal)
  if (length(cs)<catNum) {
    cs <- colorRampPalette(cs)(catNum)
  }


  comNet <- ggraph(uig, layout=layout)

  if (edgeLink){
    comNet <- comNet+
    geom_edge_link(color="grey")
  } else {
    comNet <- comNet+
    geom_edge_diagonal(color="grey")
  }
  
  if (tag | community){
    if (hull) {
      ## Hull sometimes make groups inconsistent
      comNet <- comNet + 
                geom_mark_hull(
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
    } else {
      comNet <- comNet + 
        geom_node_point(aes(color=.data$tag,
          size=.data$size))+
          scale_color_discrete(name="Tag")      
    }

  } else {
  
    if (hull) {
      comNet <- comNet + 
        geom_node_point(aes(color=col, size=size))+
        scale_color_discrete(name="Group")
      comNet <- comNet + 
        geom_mark_hull(
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
  ret@net <- comNet
  ret@numWords <- length(V(uig))
  if (returnClass) {return(ret)}
  comNet
}


#' plotDynamic
#' 
#' @description List network of words using graphlayouts::layout_as_dynamic
#' @details The function accepts the list of biotext object storing inferred networks.
#' The networks are aligned by the specific layout and plotted.
#' 
#' @param listOfNets list consisting results of wc* functions (plotType="network")
#' @param concat "union" or "intersection"
#' @param alpha pass to layout_as_dynamic
#' @param titles title to be shown on plot
#' @param tag show tag on plot
#' @param useDynamic use layout_as_dynamic
#' 
#' @export
#' @examples
#' library(igraph)
#' wc1 <- refseq(c("DDX41","IRF3","XRCC1","ERCC1","ERCC2","ERCC3"), plotType="network")
#' wc2 <- refseq(c("DDX41","PNKP","XRCC1","COPA","CD4","NLRP3"), plotType="network")
#' compare <- plotDynamic(list(wc1, wc2))
#' @return plot comparing gene clusters
#' @importFrom dplyr arrange
#' @importFrom tidygraph as_tbl_graph activate
#' @importFrom graphlayouts layout_as_dynamic
plotDynamic <- function(listOfNets,concat="union",alpha=0.5,titles=NULL,tag=FALSE,
  useDynamic=TRUE){

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
    allNodes <- allNodes |> unique()
    for (e in seq_along(listOfNets)){
        tmpadd <- setdiff(allNodes, names(V(listOfNets[[e]]@igraph)))
        tmpg <- add_vertices(listOfNets[[e]]@igraph,
                                      length(tmpadd),
                                      attr=list(name=tmpadd))
        tmpg <- tmpg |> 
        as_tbl_graph() |> 
        activate("nodes") |> 
        arrange(.data$name)
         igList[[e]] <- tmpg
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

  if (useDynamic) {
    xy <- layout_as_dynamic(igList, alpha = alpha)
  } else {
    xy1 <- layout_nicely(igList[[1]])
    xy <- list()
    for (i in seq_along(igList)) {
      xy[[i]] <- xy1
    }
  }

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