#' compareWordNet
#' 
#' compare two gene clusters based on words 
#' 
#' @param listOfNets list consisting results of wc* functions (plotType="network")
#' @param titles title to be shown on plot
#' @param layout layout
#' @param hull show category by hull
#' @param size node size
#' @param tag show tag on plot
#' @param tagLevel words that tagged in this many networks will be included
#' @param conc concavity parameter
#' @param edgeLink whether to use link or diagonal
#' @param freq concatenate frequency and show it on plot
#' @param freqMean how to concatenate frequency (TRUE: mean, FALSE: sum)
#' @param scaleRange point size scaling
#' @param returnNet return only the network (ig)
#' 
#' @export
#' @examples
#' net1 <- wcGeneSummary(c("DDX41","IRF3"), plotType="network")
#' net2 <- wcGeneSummary(c("DDX41","PNKP"), plotType="network")
#' compare <- compareWordNet(list(net1, net2))
#' @return plot comparing gene clusters
#' @import ggforce
#' @importFrom stringr str_replace
compareWordNet <- function(listOfNets, titles=NULL,
                           layout="nicely", hull=FALSE, size=4, conc=1,
                           tag=FALSE, tagLevel=1, edgeLink=TRUE,
                           freq=TRUE, freqMean=FALSE, scaleRange=c(5,10),
                           returnNet=FALSE) {
  listOfIGs <- list()
  listOfNodes <- list()

  if (is.null(titles)){
    titles <- c()
    for (e in seq_along(listOfNets)){
      titles <- c(titles, paste0("title",e))
    }
  }

  for (e in seq_along(listOfNets)) {
      listOfIGs[[e]] <- listOfNets[[e]]$ig
      listOfNodes[[e]] <- names(V(listOfNets[[e]]$ig))
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

  if (freq) {
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
  }

  col <- c()
  for (node in names(V(uig))){
      if (node %in% commonNodes) {
          col <- c(col, "Common")
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
          tmpcol <- paste(tmpcol, collapse="_")
          col <- c(col, tmpcol)
      }
  }
  V(uig)$col <- col
  if (returnNet){
    return(uig)
  }
  
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
  comNet <- ggraph(uig, layout=layout)

  if (edgeLink){
    comNet <- comNet+
    geom_edge_link()
  } else {
    comNet <- comNet+
    geom_edge_diagonal()
  }
  
  if (tag){
    comNet <- comNet + 
              geom_mark_hull(
                aes(x,y,group=tag,fill=tag,
                filter=!is.na(tag)),
                concavity=conc,
                alpha=0.25, na.rm=FALSE,
                show.legend=TRUE,
                inherit.aes=TRUE)
    ## TODO:
    ## specifying label produces an error,
    ## thus show.legend=TRUE is specified
    if (freq) {
      comNet <- comNet + 
        geom_node_point(aes(color=col, size=Freqs))+
          scale_color_discrete(name="Group")
    } else {
      comNet <- comNet + 
        geom_node_point(aes(color=col), size=size)+
          scale_color_discrete(name="Group")
    }
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
      if (freq) {
        comNet <- comNet + geom_node_point(aes(size=Freqs))
      } else {
        comNet <- comNet + geom_node_point(size=size)
      }
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
        show.legend = FALSE
      )
    } else {
      if (freq) {
        comNet <- comNet + geom_node_point(aes(color=col, size=Freqs))
      } else {
        comNet <- comNet + geom_node_point(aes(color=col), size=size)
      }
      ## You can change it later
      if (length(listOfNets)==2){
        cs <- c("tomato","steelblue")
      } else if (length(listOfNets)==3){
        cs <- c("tomato","gold","steelblue")
      } else {
        cs <- palette()
      }
      comNet <- comNet +
        scale_color_manual(name="Group",values=cs)
    }
  }

  comNet +
  geom_node_text(
    aes(label=name),
    check_overlap=TRUE, repel=TRUE,# size = labelSize,
    bg.color = "white", segment.color="black",
    bg.r = .15, show.legend=FALSE)+
  scale_size(range=scaleRange, name="Frequency")+
  theme_graph()
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
        allNodes <- c(allNodes, V(n$ig)$name)
    }
    igList <- list()
    for (e in seq_along(listOfNets)){
        tmpadd <- setdiff(allNodes, names(V(listOfNets[[e]]$ig)))
        igList[[e]] <- add_vertices(listOfNets[[e]]$ig,
                                      length(tmpadd),
                                      attr=list(name=tmpadd))
    }
  } else {
    uniqNodeNames <- list()
    for (e in seq_along(listOfNets)){
        uniqNodeNames[[e]] <- names(V(listOfNets[[e]]$ig))
    }
    uniqNodeNames <- Reduce(intersect, uniqNodeNames)
    igList <- list()
    for (e in seq_along(listOfNets)){
        igList[[e]] <- induced_subgraph(listOfNets[[e]]$ig,
          names(V(listOfNets[[e]]$ig)) %in% uniqNodeNames)
    }
  }

  xy <- layout_as_dynamic(igList, alpha = alpha)
  pList <- vector("list",length(igList))

  for(i in seq_along(igList)){
    pList[[i]] <- ggraph(igList[[i]],layout="manual",x=xy[[i]][,1],y=xy[[i]][,2])+
        geom_edge_link0(edge_width=0.6,edge_colour="grey66")
    if (tag){
      pList[[i]] <- pList[[i]] + geom_node_point(aes(size=Freq,fill=tag),
        shape=21, show.legend=FALSE)
    } else {
      pList[[i]] <- pList[[i]] + geom_node_point(aes(size=Freq,fill=Freq),
        shape=21, show.legend=FALSE)+
                    scale_fill_gradient(low="blue",high="red",
                           name = "Frequency")
    }
      pList[[i]] <- pList[[i]] +
        geom_node_text(aes(label=name),repel = TRUE,bg.color="white")+
        theme_graph()+
        theme(legend.position="bottom")+
        labs(title=titles[i])
  }
  Reduce("+",pList)
}