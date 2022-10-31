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
#' 
#' @export
#' @examples
#' wc1 <- wcGeneSummary("DDX41", plotType="network")
#' wc2 <- wcGeneSummary("IRF3", plotType="network")
#' compare <- compareWordNet(list(wc1, wc2))
#' @return plot comparing two gene clusters
#' @import ggforce
#' @importFrom stringr str_replace
compareWordNet <- function(listOfNets, titles=NULL,
                           layout="nicely", hull=FALSE, size=4, conc=1,
                           tag=FALSE, tagLevel=1) {
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
  uig <- Reduce(igraph::union, listOfIGs)
  nodeAttr <- names(get.vertex.attribute(uig))
  tagName <- nodeAttr[grepl("tag", nodeAttr)]
  
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


  col <- c()
  for (node in names(V(uig))){
      if (node %in% commonNodes) {
          col <- c(col, "Common")
      } else {
          tmpcol <- c()
          for (e in seq_along(listOfIGs)){
              if (node %in% names(V(listOfIGs[[e]]))) {
                  tmpcol <- c(tmpcol, e)
              }
          }
          tmpcol <- paste(tmpcol, collapse="_")
          col <- c(col, tmpcol)
      }
  }
  V(uig)$col <- col
  
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
  
  
  comNet <- ggraph(uig, layout=layout) +
    geom_edge_link()+
    geom_node_text(
      aes_(label=~name),
      check_overlap=TRUE, repel=TRUE,# size = labelSize,
      bg.color = "white", segment.color="black",
      bg.r = .15, show.legend=FALSE)+
    theme_graph()
  
  if (tag){
    comNet <- comNet + 
      geom_node_point(aes(color=col), size=4)+
        scale_color_discrete(name="Group")
    comNet + geom_mark_hull(
      aes(comNet$data$x,
          comNet$data$y,
          group = contag,
          fill=contag,
          label=contag,
          filter = !is.na(contag)),
      concavity = conc,
      # expand = unit(2, "mm"),
      alpha = 0.25,
      na.rm = TRUE,
      # label.fill="transparent",
      show.legend = FALSE
    )
    
  } else {
  
    if (hull) {  
      comNet + 
        geom_node_point(size=size)+
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
      comNet + geom_node_point(aes(color=col), size=4)+
        scale_color_discrete(name="Group")
    }
  }
}