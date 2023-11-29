#' split_by_ea
#' 
#' used internally for splitting the gene list by EA and returns the list of networks
#' 
#' @noRd
#' 
split_by_ea <- function(args) {
    if (args$keyType!="ENTREZID"){
        geneList <- AnnotationDbi::select(args$orgDb,
            keys = args$geneList, columns = c("ENTREZID"),
            keytype = args$keyType)$ENTREZID
        geneList <- geneList[!is.na(geneList)] |> unique()
    } else {
    	geneList <- args$geneList
    }
    if (args$splitByEA=="kegg") {
    	enr_res <- clusterProfiler::enrichKEGG(geneList)
    } else {
    	enr_res <- ReactomePA::enrichPathway(geneList)
    }
    sig_thresh <- args$genePathPlotSig
    if (dim(subset(enr_res@result, p.adjust<sig_thresh))[1]==0) {
        stop("No enriched term found.")
    }
    enr_genes <- enr_res@result %>% data.frame() %>% 
        filter(p.adjust<sig_thresh)
    gene_list <- lapply(enr_genes$geneID, function(x) strsplit(x, "/"))
    names(gene_list) <- enr_genes$Description
    
    ## Those genes not in the enrichment analysis results
    no_enr <- geneList[!(geneList %in% unique(unlist(gene_list)))]
    if (length(no_enr)!=0) {
	    gene_list[["no_enrichment"]] <- no_enr	
    }
    qqcat("Total of @{length(gene_list)} pathways, including non-enrichment terms\n")
    args2 <- args
    btg_list <- lapply(gene_list, function(tmp_gene_list) {
    	args2$geneList <- tmp_gene_list |> unlist()
    	args2$splitByEA <- NULL
    	args2$keyType <- "ENTREZID"
    	do.call(refseq, args2)
    })
    names(btg_list) <- names(gene_list)
    return(btg_list)
}


#' changeLayout
#' @param g biotext object
#' @param layout_func layout function in igraph
#' @examples refseq(c("IRF3","PNKP","DDX41")) |> changeLayout(igraph::layout_nicely)
#' @export
#' @return biotext class object
changeLayout <- function(g, layout_func) {
  lyt <- do.call(layout_func, list(graph=g@igraph))
  # lyt <- eval(parse(text=layout_func))(g@igraph)
  g@net$data$x <- lyt[,1]
  g@net$data$x <- lyt[,2]
  g
}

#' geom_node_shadowtext
#' 
#' Plot shadowtext at node position
#' 
#' @export
#' @param mapping aes mapping
#' @param data data to plot
#' @param position positional argument
#' @param show.legend whether to show legend
#' @param ... passed to `params` in `layer()` function
#' @return geom
#' @importFrom shadowtext GeomShadowText
#' @examples
#' nodes <- data.frame(name=c("hsa:1029","hsa:4171"),
#'                    x=c(1,1),
#'                    xmin=c(-1,-1),
#'                    xmax=c(2,2),
#'                    y=c(1,1),
#'                    ymin=c(-1,-1),
#'                    ymax=c(2,2))
#' edges <- data.frame(from=1, to=2)
#' graph <- tidygraph::tbl_graph(nodes, edges)
#' plt <- ggraph::ggraph(graph, layout="manual", x=x, y=y) +
#'  geom_node_shadowtext(aes(label=name))
geom_node_shadowtext <- function(mapping = NULL, data = NULL,
                           position = 'identity',
                           show.legend = NA, ...) {
  params <- list(na.rm = FALSE, ...)

  mapping <- c(mapping, aes(x=.data$x, y=.data$y))
  class(mapping) <- "uneval"

  layer(
    data = data, mapping = mapping, stat = StatFilter, geom = GeomShadowText,
    position = position, show.legend = show.legend, inherit.aes = FALSE,
    params = params
  )
}

#' obtainMatrix
#' 
#' obtain matrix of words-to-words relationship
#' 
#' @noRd
#' @return list of graph and osplot object
#' 
obtainMatrix <- function(ret, bn, R, DTM, freqWords,
    corThresh, cooccurrence, onWholeDTM, numWords, autoThresh=FALSE) {

    retList <- list()
    if (bn) {
      qqcat("bn specified, R=@{R}\n")
      # To avoid computaitonal time, subset to numWords
      bnboot <- bnlearn::boot.strength(
          data.frame(
              DTM[, colnames(DTM) %in% freqWords]),
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
      ## TODO: speed up calculation using Rcpp
      if (onWholeDTM) {## To obtain correlation for all the words
          corInput <- DTM
      } else {
          corInput <- DTM[, colnames(DTM) %in% freqWords]
      }
      if (cooccurrence) {
          corData <- t(corInput) %*% corInput
      } else {
          corData <- abs(cor(corInput))
      }
      if (autoThresh) {
        qqcat("Ignoring corThresh, automatically determine the value\n")
        tryVals <- seq(min(corData, na.rm=TRUE), max(corData, na.rm=TRUE),
        	(max(corData, na.rm=TRUE)-min(corData, na.rm=TRUE)) / 10)
        noden_list <- lapply(tryVals, function(tmpCorThresh) {
            tmpCorData <- corData
            tmpCorData[tmpCorData<tmpCorThresh] <- 0
            tmpCoGraph <- graph.adjacency(tmpCorData, weighted=TRUE,
                      mode="undirected", diag = FALSE)
            tmpCoGraph <- induced.subgraph(tmpCoGraph,
                degree(tmpCoGraph) > 0)
            length(V(tmpCoGraph))
        })
        names(noden_list) <- tryVals
        noden_list <- noden_list |> unlist()
        tmp_thresh <- noden_list[noden_list >= numWords] |> names()
        corThresh <- tmp_thresh[length(tmp_thresh)] |> as.numeric()
        qqcat("threshold = @{corThresh}\n")
      }
      ret@corMat <- corData
      ret@corThresh <- corThresh
      ## Set correlation below threshold to zero
      corData[corData<corThresh] <- 0
      coGraph <- graph.adjacency(corData, weighted=TRUE,
                  mode="undirected", diag = FALSE)
  }
  retList[["ret"]] <- ret
  retList[["coGraph"]] <- coGraph
  return(retList)
}

#' appendEdges
#' 
#' function to append edges of ggraph
#' 
#' @noRd
appendEdges <- function(netPlot,
  bn, edgeLink, edgeLabel, showLegend,
  fontFamily) {

    if (bn){
      if (edgeLink){
          if (edgeLabel){
              netPlot <- netPlot +
                          geom_edge_link(
                              aes(width=.data$weight,
                              color=.data$edgeColor,
                              label=.data$weightLabel),
                              angle_calc = 'along',
                              family=fontFamily,
                              label_dodge = unit(2.5, 'mm'),
                              # arrow = arrow(length = unit(4, 'mm')), 
                              # start_cap = circle(3, 'mm'),
                              # end_cap = circle(3, 'mm'),
                              alpha=0.5,
                              show.legend = showLegend)
          } else {
              netPlot <- netPlot +
                          geom_edge_link(aes(width=.data$weight,
                              color=.data$edgeColor),
                              # arrow = arrow(length = unit(4, 'mm')), 
                              # start_cap = circle(3, 'mm'),
                              # end_cap = circle(3, 'mm'),
                              alpha=0.5, show.legend = showLegend)
          }
      } else {
          if (edgeLabel){
              netPlot <- netPlot +
                          geom_edge_diagonal(
                              aes(width=.data$weight,
                              color=.data$edgeColor,
                              label=.data$weightLabel),
                              angle_calc = 'along',
                              family=fontFamily,
                              label_dodge = unit(2.5, 'mm'),
                              # arrow = arrow(length = unit(4, 'mm')), 
                              # start_cap = circle(3, 'mm'),
                              # end_cap = circle(3, 'mm'),
                              alpha=0.5,
                              show.legend = showLegend)
          } else {
              netPlot <- netPlot +
                          geom_edge_diagonal(aes(width=.data$weight,
                              color=.data$edgeColor),
                              # arrow = arrow(length = unit(4, 'mm')), 
                              # start_cap = circle(3, 'mm'),
                              # end_cap = circle(3, 'mm'),                                    
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
                              label=.data$weightLabel),
                              family=fontFamily,
                              angle_calc = 'along',
                              label_dodge = unit(2.5, 'mm'),
                              alpha=0.5,
                              show.legend = showLegend)                        
          } else {
              netPlot <- netPlot +
                          geom_edge_link(aes(width=.data$weight,
                              color=.data$edgeColor),
                              alpha=0.5, show.legend = showLegend)
          }
      } else {
          if (edgeLabel){
              netPlot <- netPlot +
                          geom_edge_diagonal(
                              aes(width=.data$weight,
                              color=.data$edgeColor,
                              label=.data$weightLabel),
                              angle_calc = 'along',
                              label_dodge = unit(2.5, 'mm'),
                              alpha=0.5,
                              family=fontFamily,
                              show.legend = showLegend)                        
          } else {
              netPlot <- netPlot +
                          geom_edge_diagonal(aes(width=.data$weight,
                              color=.data$edgeColor),
                              alpha=0.5, show.legend = showLegend)                
          }
      }
  }
}




#' appendNodesAndTexts
#' 
#' function to append nodes and texts based on parameters
#' to ggraph. Colorize by these combinations based on the parameter
#' 
#' words (tag / community [disc]) - other nodes (category [disc])
#' words (frequency [cont]) - other nodes (category [disc])
#' words (category [disc]) - other nodes (category [disc])
#' 
#' @noRd
appendNodesAndTexts <- function(netPlot,tag,colorize,nodePal,
  showLegend,catColors,pal,fontFamily,colorText,scaleRange,useSeed,ret,tagColors,
  discreteColorWord){

  if (is.null(catColors)) {
      catNum <- length(unique(V(ret@igraph)$nodeCat))
      if (catNum>2) {
        catColors <- RColorBrewer::brewer.pal(catNum, "Dark2")        
      } else {
        catColors <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(catNum)]
      }
      names(catColors) <- unique(V(ret@igraph)$nodeCat)
  }
  if (is.null(tagColors)) {
      tagNum <- length(unique(V(ret@igraph)$tag))
      if (tagNum>2) {
        tagColors <- RColorBrewer::brewer.pal(tagNum, "Dark2")        
      } else {
        tagColors <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(tagNum)]
      }
      names(tagColors) <- unique(V(ret@igraph)$tag)
  }

  if (tag!="none") { ## use pvpick
      useTagColors <- tagColors[ netPlot$data$tag ]
      netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$tag),
                                          show.legend = showLegend)+
                           scale_color_manual(values=tagColors)
  } else {
      if (colorize) {## Colorize by node category (except for words)
          if (discreteColorWord){
            ## All the nodes are discrete (e.g. Words, Genes, ...)
            useCatColors <- catColors[ netPlot$data$nodeCat ]
            netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$nodeCat),
                                          show.legend = showLegend)+
                           scale_color_manual(values=catColors)

          } else {
            ## colorize points of texts
            netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$Freq),
                                                show.legend = showLegend)+
                                scale_color_gradient(low=pal[1],high=pal[2],
                                                  na.value="grey50")
            ## colorize the other points
            useCatColors <- catColors[ netPlot$data[ netPlot$data$nodeCat != "Words", ]$nodeCat ]
            if (sum(is.na(useCatColors))>0) {qqcat("Some color contains NA, proceeding\n")}
            netPlot <- netPlot + geom_node_point(aes(size=.data$Freq,
                                                filter=.data$nodeCat != "Words"),
                                                show.legend=FALSE,
                                                color=useCatColors)
          }
      } else {
          netPlot <- netPlot + geom_node_point(aes(size=.data$Freq, color=.data$Freq),
                                              show.legend = showLegend)+
                               scale_color_gradient(low=pal[1],high=pal[2],na.value="grey50",
                                                    name = "Frequency")
      }
  }
  if (colorText){
      if (tag!="none") {
          netPlot <- netPlot + 
              geom_node_text(aes(label=.data$name, size=.data$Freq, color=.data$tag),
                # color=useTagColors,
                  check_overlap=TRUE, repel=TRUE,
                  bg.color = "white", segment.color="black",family=fontFamily,
                  bg.r = .15, show.legend=showLegend)
      } else {
        if (colorize) {
          if (discreteColorWord) {
            netPlot <- netPlot + 
              geom_node_text(aes(label=.data$name, size=.data$Freq, color=.data$nodeCat),
                  check_overlap=TRUE, repel=TRUE,
                  bg.color = "white", segment.color="black",family=fontFamily,
                  bg.r = .15, show.legend=showLegend)
          } else {
        ## [TODO] discrete and continuous scale in same ggplot is discouraged
        ##        use ggnewscale?
        ## [TODO] repel not work for multiple layers
          netPlot <- netPlot + 
              geom_node_text(aes(
                  label=.data$name,
                  size=.data$Freq,
                  color=.data$Freq),
                  check_overlap=TRUE, repel=TRUE,
                  bg.color = "white", segment.color="black",family=fontFamily,
                  bg.r = .15, show.legend=showLegend)

          layerNum <- length(netPlot$layers)
          geom_param_list <- netPlot$layers[[layerNum]]$geom_params
          build <- ggplot_build(netPlot)$data[[layerNum]]
          build[ netPlot$data$nodeCat !="Words", ]$colour <- useCatColors

          aes_list <- netPlot$layers[[layerNum]]$mapping
          aes_list["filter"] <- NULL
          geom_param_list[["repel"]] <- TRUE

          geom_param_list[["seed"]] <- useSeed
          netPlot$layers[[layerNum]]$geom_params[["seed"]] <- useSeed
          geom_param_list[["show.legend"]] <- FALSE
          geom_param_list["color"] <- NULL;
          aes_list["label"] <- NULL
          aes_list["color"] <- NULL
          aes_list["colour"] <- NULL
          geom_param_list["na.rm"] <- NULL
          aes_params <- netPlot$layers[[layerNum]]$aes_params

          netPlot <- netPlot + do.call(geom_node_text,
                         c(aes_params,
                          list(mapping=aes_list),
                          geom_param_list, list(
                           label=build$label,
                            color=build$colour)))

          # netPlot <- netPlot + geom_node_text(aes(label=.data$name, size=.data$Freq,
          #         filter=.data$tag != "Words"), color=useCatColors,
          #         check_overlap=TRUE, repel=TRUE, seed=useSeed,
          #         bg.color = "white", segment.color="black",family=fontFamily,
          #         bg.r = .15, show.legend=FALSE)
        }

      } else {
          netPlot <- netPlot + 
              geom_node_text(aes(label=.data$name, size=.data$Freq, color=.data$Freq),
                  check_overlap=TRUE, repel=TRUE,# size = labelSize,
                  bg.color = "white", segment.color="black",family=fontFamily,
                  bg.r = .15, show.legend=showLegend)
      }
    }
  } else {
      netPlot <- netPlot +
                  geom_node_text(aes(label=.data$name, size=.data$Freq),
                      check_overlap=TRUE, repel=TRUE,# size = labelSize,
                      color = "black",
                      bg.color = "white", segment.color="black",family=fontFamily,
                      bg.r = .15, show.legend=showLegend) 
    }
  # netPlot <- netPlot + scale_size(range=scaleRange, name="Frequency")
  netPlot
}



#'
#' connectGenes
#'
#' When an interesting word is found in two or more networks,
#' connect the words and gene names and return the graph.
#' Note that genePlot must be set to TRUE in all the networks.
#'
#' @param nets named list of nets (if no names are given, automatically set)
#' @param query_word word to connect
#' @param return_tbl_graph return tbl_graph
#' @param neighbors obtain neighbors words for queried words, default to FALSE
#' @return igraph object
#' @export
connectGenes <- function(nets, query_word, return_tbl_graph=FALSE,
    neighbors=FALSE) {
  words <- NULL
  if (is.null(names(nets))) {
    names(nets) <- paste0("Net",seq_len(length(nets)))
  }
  k <- 1
  for (net in nets) {
    if (dim(net@geneMap)[1]==0) {
      stop("No gene mapping present in object")
    }
    if (neighbors) {
        el <- induced.subgraph(net@igraphRaw, neighbors(net@igraphRaw, query_word)) |> as_edgelist()
    }
    
    mp <- net@geneMap
    mp <- data.frame(mp[mp[,1]==query_word,])
    mp$type <- names(nets)[k]
    
    words <- rbind(words, mp[,c(1,2)])
    words <- rbind(words, mp[,c(2,3)]|>`colnames<-`(colnames(mp[,c(1,2)])))
    if (neighbors) {
        nbs <- NULL
        for (i in c(unique(el[,2],el[,1]))) {
            nbs <- rbind(nbs, c(i, query_word))
        }
        nbs <- nbs |> data.frame() |> `colnames<-`(colnames(mp[,c(1,2)]))
        words <- rbind(words, nbs)
    }
    k <- k+1
  }
  g <- graph_from_edgelist(as.matrix(words), directed = FALSE)
  if (return_tbl_graph) {
    return(g |> as_tbl_graph())
  } else {
    return(g)
  }
}

#'
#' returnQuanteda
#' 
#' simple function to generate dfm and export to wc* functions.
#' It did not utilize full functionality of quanteda, and 
#' will add better support.
#' @param ret osplot object
#' @param quantedaArgs args to passed to tokens function
#' @param numWords number of words
#' @param ngram tokens_ngrams
#' @param tfidf use dfm_tfidf or not
#' @param filterWords filtered words
#' @param additionalRemove user-specified filter words
#' @param collapse collapse sentences
#' @return osplot object after filtering using quanteda
returnQuanteda <- function(ret, quantedaArgs,numWords,ngram,
                           filterWords,additionalRemove, tfidf, collapse) {
    
    if (length(quantedaArgs)==0) {
        quantedaArgs <- list(
            remove_punct=TRUE,
            remove_symbols=TRUE,
            remove_numbers=TRUE,
            remove_url=TRUE,
            remove_separators=TRUE,
            split_hyphens=FALSE
        )
    }
    if (collapse) {
      docs <- quanteda::corpus(paste(ret@rawText$text, collapse=" "))
    } else {
      docs <- quanteda::corpus(ret@rawText$text)      
    }
    ret@corpusQuanteda <- docs
    
    ## case_insensitive=TRUE
    
    quantedaArgs[["x"]] <- docs
    tokens <- do.call(quanteda::tokens, quantedaArgs) |>
        quanteda::tokens_remove(quanteda::stopwords("english"))

    if (length(additionalRemove[!is.na(additionalRemove)])!=0) {
        tokens <- quanteda::tokens_remove(tokens, additionalRemove)
    }    
    if (length(filterWords[!is.na(filterWords)])!=0) {
        tokens <- quanteda::tokens_remove(tokens, filterWords)
    }
    if (ngram!=1) {
        tokens <- quanteda::tokens_ngrams(tokens, n=ngram)
    }

    freqWordsDFM <- quanteda::dfm(tokens)
    if (tfidf) {
      freqWordsDFM <- quanteda::dfm_tfidf(freqWordsDFM)
    }
    ret@dfm <- freqWordsDFM
    # freqWordsAll <- sort(quanteda::featfreq(freqWordsDFM),
    #                   decreasing=TRUE)
    # freqWords <- names(freqWordsAll[1:numWords])
    return(ret)
}


#' convertMetaCyc
#'
#' this function needs taxonomizr package to be installed,
#' and needs download of databases by prepareDatabase(getAccessions=FALSE)
#' 
#' @param ids tax ids from metacyc
#' @param onlySpecies parse only species
#' @return coverted species name
#' @examples
#' \dontrun{convertMetaCyc("TAX-9606")}
#' @export
convertMetaCyc <- function (ids, onlySpecies=FALSE) {
  convIds <- sapply(ids, function(x) if (grepl("TAX-",x)) unlist(strsplit(x,"-")[[1]])[2] else x)
  if (onlySpecies) {
    parsed <- apply(taxonomizr::getTaxonomy(convIds), 1, function(x) x[7])
  } else {
    parsed <- apply(taxonomizr::getTaxonomy(convIds), 1, function(x) paste(x, collapse=";"))
  }
  as.character(parsed)
}


#' clearPath
#' 
#' clear MetaCyc pathway strings
#' 
#' @param ex vector of text
#' @noRd
clearPath <- function (ex) {
    ex <- gsub("|FRAME: ", "", ex,useBytes = TRUE)
    ex <- gsub("FRAME: ", "", ex,useBytes = TRUE)
    ex <- gsub("|CITS: ", "", ex,useBytes = TRUE)
    ex <- gsub("CITS: ", "", ex,useBytes = TRUE)
    ex <- gsub("\\[[^][]*]", "", ex,useBytes = TRUE)
    ex <- gsub("\"", "", ex,useBytes = TRUE)
    ## Clean HTML tags
    ex <- gsub("<.*?>", "", ex,useBytes = TRUE)
    ## lastly remove \\|
    ex <- gsub("\\|","",ex,useBytes = TRUE)
    ex
}


#'
#' parseMetaCycPathway
#' 
#' parse MetaCyc "pathways.dat"
#' 
#' 
#' @param file path to pathways.dat
#' @param candSp species to grepl
#' @param withTax parse taxonomy information
#' @param noComma no comma separated when taxonomy parsing
#' @param clear delete HTML tags and some symbols
#' @return data.frame of MetaCyc pathway information
#' @examples
#' file <- "pathways.dat"
#' \dontrun{parseMetaCycPathway(file, candSp="all")}
#' @export
#' 
parseMetaCycPathway <- function(file, candSp, withTax=FALSE, noComma=FALSE, clear=FALSE) {
  flg <- FALSE
  allFlg <- FALSE
  allmeta <- NULL
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (startsWith(line,"UNIQUE-ID - ")) {
      com <- NA
      commn <- NA
      spec <- NULL
      taxr <- NULL
      pwy <- gsub("UNIQUE-ID - ","",line,useBytes = TRUE)
      flg <- TRUE
    }
    if (flg) {
      if (startsWith(line, "/")) {
        com <- c(com, line)
      }
      if (startsWith(line, "COMMON-NAME")) {
        commn <- gsub("COMMON-NAME - ","",line,useBytes = TRUE)
      }
      if (startsWith(line, "SPECIES - ")) {
        spec <- c(spec, gsub("SPECIES - ","",line,useBytes = TRUE))
      }
      if (startsWith(line, "TAXONOMIC-RANGE - ")) {
        taxr <- c(taxr, gsub("TAXONOMIC-RANGE - ","",line,useBytes = TRUE))
      }
      if (startsWith(line,"//")) {
        coms <- paste(com[!is.na(com)], collapse=" ")
        coms <- gsub("/","",coms,useBytes = TRUE)

        if (!noComma) {
          if (length(spec)!=0) {spec <- paste0(spec, collapse=",")} else {spec <- ""}
          if (length(taxr)!=0) {taxr <- paste0(taxr, collapse=",")} else {taxr <- ""}          
        } else {
          if (length(spec)!=0) {} else {spec <- ""}
          if (length(taxr)!=0) {} else {taxr <- ""}
        }
        if (length(candSp)==1) {
          if (candSp=="all") {
             allFlg <- TRUE
          }
        }

        if (!allFlg) {
          if (grepl(paste(candSp,collapse="|"),coms)) {
            if (withTax) {
              if (noComma) {
                for (sp in spec) {
                  for (tax in taxr) {
                    allmeta <- rbind(allmeta, c(pwy, coms, commn, sp, tax))
                  }
                }
              } else {
                allmeta <- rbind(allmeta, c(pwy, coms, commn, spec, taxr))
              }
            } else {
              allmeta <- rbind(allmeta, c(pwy, coms, commn))
            }
          }
        } else {
            if (withTax) {
              if (noComma) {
                for (sp in spec) {
                  for (tax in taxr) {
                    allmeta <- rbind(allmeta, c(pwy, coms, commn, sp, tax))
                  }
                }
              } else {
                allmeta <- rbind(allmeta, c(pwy, coms, commn, spec, taxr))
              }
            } else {
              allmeta <- rbind(allmeta, c(pwy, coms, commn))
            }
        }
        flg <- FALSE
      }
    }
  }
  close(con)
  if (withTax) {
    allmeta <- data.frame(allmeta) |> `colnames<-`(c("pathwayID","text","commonName","species","taxonomicRange"))
  } else {
    allmeta <- data.frame(allmeta) |> `colnames<-`(c("pathwayID","text","commonName"))
  }
  # allmeta |> dim()
  if (!allFlg) {
    queries <- NULL
    for (q in candSp) {
      queries <- cbind(queries, grepl(q, allmeta$text))
    }
    allmeta$query <- apply(queries, 1, function(x) paste0(candSp[x], collapse=","))
  }

  if (clear) {
    allmeta$text <- clearPath(allmeta$text)
  }
  allmeta$text <- iconv(allmeta$text, "ASCII", "UTF-8")
  return(allmeta)
}

#' retFiltWords
#' return filtered words
#' @param useFil use filter
#' @param filType "above" or "below"
#' @param filNum number
#' @noRd

retFiltWords <- function(useFil, filType, filNum) {

    data_env <- new.env(parent = emptyenv())
    load(system.file("extdata", "sysdata.rda", package = "biotextgraph"),
        envir=data_env)

    if (filType=="above" | filType==">") {
        if (useFil=="GS_TfIdf") {

          allTfIdfGeneSummary <- data_env[["allTfIdfGeneSummary"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummary[
            allTfIdfGeneSummary$tfidf > filNum,]$word
            
        } else if (useFil=="BSDB_TfIdf"){

          allTfIdfBSDB <- data_env[["allTfIdfBSDB"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDB[
            allTfIdfBSDB$tfidf > filNum,]$word

        } else if (useFil=="GS_TfIdf_Max"){

          allTfIdfGeneSummaryMax <- data_env[["allTfIdfGeneSummaryMax"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummaryMax[
            allTfIdfGeneSummaryMax$tfidf > filNum,]$word

        } else if (useFil=="BSDB_TfIdf_Max"){

          allTfIdfBSDBMax <- data_env[["allTfIdfBSDBMax"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDBMax[
            allTfIdfBSDBMax$tfidf > filNum,]$word

        } else if (useFil=="GS_Freq"){

          allFreqGeneSummary <- data_env[["allFreqGeneSummary"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allFreqGeneSummary[
            allFreqGeneSummary$freq > filNum,]$word

        } else if (useFil=="BSDB_Freq"){

          allFreqBSDB <- data_env[["allFreqBSDB"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allFreqBSDB[
            allFreqBSDB$freq > filNum,]$word

        } else {
          stop("Please specify useFil")
        }
    } else {
        if (useFil=="GS_TfIdf") {

          allTfIdfGeneSummary <- data_env[["allTfIdfGeneSummary"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummary[
            allTfIdfGeneSummary$tfidf < filNum,]$word

        } else if (useFil=="BSDB_TfIdf"){

          allTfIdfBSDB <- data_env[["allTfIdfBSDB"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDB[
            allTfIdfBSDB$tfidf < filNum,]$word

        } else if (useFil=="GS_TfIdf_Max"){

          allTfIdfGeneSummaryMax <- data_env[["allTfIdfGeneSummaryMax"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummaryMax[
            allTfIdfGeneSummaryMax$tfidf < filNum,]$word

        } else if (useFil=="BSDB_TfIdf_Max"){

          allTfIdfBSDBMax <- data_env[["allTfIdfBSDBMax"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDBMax[
            allTfIdfBSDBMax$tfidf < filNum,]$word

        } else if (useFil=="GS_Freq"){
          allFreqGeneSummary <- data_env[["allFreqGeneSummary"]]
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allFreqGeneSummary[
            allFreqGeneSummary$freq < filNum,]$word

        } else if (useFil=="BSDB_Freq"){

          allFreqBSDB <- data_env[["allFreqBSDB"]]
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allFreqBSDB[
            allFreqBSDB$freq < filNum,]$word

        } else {
          stop("Please specify useFil")
        }
    }
    return(filterWords)
}


#' removeAloneNumbers
#' 
#' @noRd
#' 
removeAloneNumbers <- 
    function (x) PlainTextDocument(
        gsub('\\s*(?<!\\B|-)\\d+(?!\\B|-)\\s*', " ",
                      x, perl=TRUE))

#' preserveDict
#' 
#' 
#' @param docs corpus
#' @param ngram ngram params
#' @param numOnly remove number only
#' @param stem use stemming
#' 
#' @noRd
#' 
preserveDict <- function(docs, ngram, numOnly, stem) {
    NgramTokenizer <- function(x)
    unlist(lapply(ngrams(words(x), ngram),
        paste, collapse = " "),
        use.names = FALSE)
  ## Changes results if perform tolower here
  # ldocs <- docs
  ldocs <- docs %>%
      tm_map(FUN=content_transformer(tolower))
  if (numOnly) {
    ldocs <- ldocs %>% tm_map(FUN=removeAloneNumbers)
    docs <- docs %>% tm_map(FUN=removeAloneNumbers)
  } else {
    ldocs <- ldocs %>% tm_map(FUN=removeNumbers)
    docs <- docs %>% tm_map(FUN=removeNumbers)
  }

  docs <- docs %>%
    tm_map(removeWords, stopwords::stopwords("english",
                                           "stopwords-iso")) %>%
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
  ldocs <- ldocs %>%
    tm_map(removeWords, stopwords::stopwords("english",
                                           "stopwords-iso")) %>%
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
  if (stem) {
    docs <- docs %>% tm_map(stemDocument) 
    ldocs <- ldocs %>% tm_map(stemDocument)
  }

  if (ngram!=1){
      docs <- TermDocumentMatrix(docs,
                                 control = list(tokenize = NgramTokenizer,
                                                tolower=FALSE))
      ldocs <- TermDocumentMatrix(ldocs,
                                    control = list(tokenize = NgramTokenizer,
                                                   tolower=TRUE))
  } else {
      docs <- TermDocumentMatrix(docs, control=list(tolower=FALSE))
      ldocs <- TermDocumentMatrix(ldocs, control=list(tolower=TRUE))
  }

  rawrn <- row.names(docs)
  lrn <- row.names(ldocs)
  dic <- list()

  for (i in rawrn) {
    if (tolower(i) %in% lrn) {
      dic[[tolower(i)]] <- i      
    }
  }
  rawWords <- names(dic)
  newWords <- unlist(dic)
  names(newWords) <- rawWords
  return(newWords)
}

#' getPubMed
#' 
#' obtain pubmed information
#' 
#' @param ret osplot object
#' @param searchQuery search query
#' @param rawQuery raw query
#' @param type abstract or title
#' @param apiKey api key
#' @param retMax retmax
#' @param perQuery query per searchQuery.
#' apiKey should be enabled as it recursively query API.
#' @import rentrez
#' 
#' @noRd

getPubMed <- function(ret, searchQuery, rawQuery,
    type="abstract", apiKey=NULL, retMax=10, sortOrder="relevance",
    perQuery=FALSE) {
    if (is.null(apiKey)){
        qqcat("Proceeding without API key\n")
        if (perQuery) {
      	    qqcat("It is strongly recommended to use API key when using",
      		    " multiple queries without a delimiter\n")
        }
    } else {
        set_entrez_key(apiKey)
    }
    if (type=="abstract") {
        obt <- "Abstract"
    } else {
        obt <- "ArticleTitle"
    }
    if (perQuery) {
    	obtainedText <- NULL
    	pmids <- NULL
    	for (tmp_query in searchQuery) {
    		qqcat("  Querying @{tmp_query}\n")
		    pubmedSearch <- entrez_search("pubmed",
		                                  term = tmp_query, 
		                                  retmax = retMax,
		                                  sort = sortOrder)
		    ret@pmids <- c(ret@pmids, pubmedSearch$ids)
		    searchResults <- entrez_fetch(db="pubmed",
		                                  pubmedSearch$ids, rettype="xml", 
		                                  parsed=FALSE)
		    parsedXML <- xmlTreeParse(as.character(searchResults))
		    charset <- xmlElementsByTagName(parsedXML$doc$children$PubmedArticleSet,
		                                    obt,
		                                    recursive = TRUE)
		    obtainedText <- c(obtainedText, as.character(xmlValue(charset)))
		    Sys.sleep(1) ## Make sure not querying in one second
    	}
    } else {
	    pubmedSearch <- entrez_search("pubmed", term = searchQuery, 
	                                  retmax = retMax, sort = sortOrder)
	    ret@pmids <- pubmedSearch$ids
	    searchResults <- entrez_fetch(db="pubmed",
	                                  pubmedSearch$ids, rettype="xml", 
	                                  parsed=FALSE)
	    parsedXML <- xmlTreeParse(as.character(searchResults))



	    charset <- xmlElementsByTagName(parsedXML$doc$children$PubmedArticleSet,
	                                    obt,
	                                    recursive = TRUE)

	    obtainedText <- as.character(xmlValue(charset))
    }

    ## Query flag is obtained by grep the query in the whole text,
    ## thus the multiple queries can be found in one article.
    
    incs <- c()
    for (i in rawQuery){
      li <- tolower(i)
      inc <- grepl(li, tolower(obtainedText), fixed=TRUE)
      inc[inc] <- i
      # inct <- grepl(li, tolower(titletext), fixed=TRUE)
      # inct[inct] <- i
      incs <- cbind(incs, inc)
    }
    obtainedDf <- data.frame(cbind(
      seq_len(length(obtainedText)),
      obtainedText,
      apply(incs, 1, function(x) paste(unique(x[x!="FALSE"]), collapse=","))
    )) |> `colnames<-`(c("ID","text","query"))
    ret@rawText <- obtainedDf
    return(ret)
}


#' findTerm
#' 
#' find queried terms in list of gene clusters and return frequency
#' 
#' @param query query words
#' @param listOfGenes named list of genes
#' @param split split the query by space, default to FALSE
#' @param ngram use ngram or not
#' @param tfidf use tfidf when creating TDM or not
#' @param keyType key type of listOfGenes
#' @param calc "sum", "mean" or "highest"
#' @param argList passed to refseq
#' @return named list of frequency
#' @export
#' @examples
#' query <- "DNA repair"
#' lg <- list()
#' lg[["sample"]] <- c("ERCC1","ERCC2")
#' findTerm(query, lg)
#' 
findTerm <- function (query, listOfGenes, split=FALSE, ngram=1,
                      tfidf=TRUE, calc="sum", keyType="SYMBOL", argList=list()) {
    qqcat("Finding query in @{length(listOfGenes)} clusters ...\n")
    if (split) {
        querySplit <- tolower(unlist(strsplit(query, " ")))
    } else {
        querySplit <- tolower(query)
    }
    frq <- list()
    for (clus in names(listOfGenes)) {
        argList[["geneList"]] <- listOfGenes[[clus]]
        argList[["keyType"]] <- keyType
        argList[["ngram"]] <- ngram
        argList[["tfidf"]] <- tfidf
        argList[["onlyTDM"]] <- TRUE
        tmptdm <- do.call("refseq", argList)
        # tmptdm <- refseq(listOfGenes[[clus]],
        #                         keyType = keyType, ngram=ngram,
        #                         tfidf=tfidf, onlyTDM=TRUE)
        querytdm <- t(as.matrix(tmptdm[Terms(tmptdm) %in% querySplit, ]))
        tmp <- rep(0, length(querySplit))
        names(tmp) <- querySplit
        if (calc=="sum"){
            tmpfrq <- apply(querytdm, 2, sum)# / dim(querytdm)[1]
        } else if (calc=="mean") {
            tmpfrq <- apply(querytdm, 2, mean)# / dim(querytdm)[1]
        } else if (calc=="highest") {
            tmpfrq <- apply(querytdm, 2, max)
        } else {
            stop("please specify sum, mean or highest")
        }
        for (i in names(tmpfrq)){
            tmp[names(tmp)==i] <- tmpfrq[i]
        }
        frq[[clus]] <- tmp
    }
    frq
}

#' returnSim
#' 
#' return similarity matrix of cluster
#' based on intersection over union of words
#' 
#' @param cllist cluster list (named vector or list)
#' @param keyType keytype
#' @param numLimit threshold for gene number limit
#'                 default to 5000
#' @param target target to query
#' @param argList parameters to pass to refseq()
#' @return similarity matrix
#' @export
#' @examples
#' ex <- returnExample()
#' returnSim(ex$color, keyType="ENSEMBL")
#' @importFrom GetoptLong qqcat
returnSim <- function (cllist, keyType="ENTREZID", numLimit=5000,
  target="refseq", argList=list()) {
    store <- list()
    if (!is.list(cllist)){
        converted <- list()
        clname <- unique(cllist)
        for (c in clname) {
            converted[[as.character(c)]] <-
                names(cllist)[cllist==c]
        }
    } else {
        converted <- cllist
    }
    qqcat("Number of clusters: @{length(converted)}\n")
    for (i in names(converted)) {
        ## Time consuming cluster
        if (length(converted[[i]])<numLimit){
            qqcat("@{i}\n")
            argList[["geneList"]] <- converted[[i]]
            argList[["keyType"]] <- keyType
            store[[i]] <-
                do.call("refseq", argList)@freqDf
        }
    }
    sim <- sapply(store, function(x) sapply(store,
                    function(y)
                        length(intersect(x$word, y$word)) /
                        length(union(x$word, y$word))))
    return(as.matrix(sim))
}

#' makeBar
#' 
#' Makeing a barplot of word frequency from queried genes
#' 
#' @examples
#' geneList <- c("DDX41")
#' makeBar(geneList)
#' @param queries gene IDs
#' @param keyType default to SYMBOL
#' @param top how many numbers of words to be shown
#' @param grad use gradient of frequency
#' @param pal palette used in barplot
#' @param textSize text size in barplot
#' @param reord order by frequency or not
#' @param flip flip the barplot (gene name in y-axis)
#' @param orgDb orgDb
#' @param retList return result of refseq
#' @param argList passed to refseq()
#' @import org.Hs.eg.db
#' @return barplot of word frequency
#' @importFrom stats reorder
#' @export
#' 
makeBar <- function(queries, top=10, keyType="SYMBOL",
                    pal=NULL, textSize=20, reord=TRUE, orgDb=org.Hs.eg.db,
                    flip=FALSE, grad=FALSE, retList=FALSE, argList=list()) {
  if (is.null(pal)) {
    # palNum <- sample(1:151,1)
    # pal <- pokepal(palNum)
    pal <- palette()
    if (length(pal)<top){
      pal <- rep(pal, ceiling(top/length(pal)))
    }
  }

  argList[["geneList"]] <- queries
  argList[["keyType"]] <- keyType
  argList[["orgDb"]] <- orgDb
  argList[["plotType"]] <- "wc"
  argList[["madeUpper"]] <- c("dna","rna",
                                  tolower(AnnotationDbi::keys(orgDb,
                                                              keytype="SYMBOL")))
  wc <- do.call("refseq",argList)
  barp <- utils::head(wc@freqDf, n=top)
  ## Need rewrite
  if (reord){
    if (grad) {
      plt <- ggplot(barp, aes(x=reorder(barp$word, barp$freq),
                              y=barp$freq, fill=barp$freq))+ scale_fill_gradient(low="blue",high="red", guide="none")
    } else {
      plt <- ggplot(barp, aes(x=reorder(barp$word, barp$freq),
                              y=barp$freq, fill=barp$word))+ scale_fill_manual(values=pal, guide="none")
    }
  } else {
    if (grad) {
      plt <- ggplot(barp, aes(x=barp$word,
                              y=barp$freq, fill=barp$freq))+ scale_fill_gradient(low="blue",high="red", guide="none")
    } else {
      plt <- ggplot(barp, aes(x=barp$word,
                              y=barp$freq, fill=barp$word))+ scale_fill_manual(values=pal, guide="none")
    }
  }     
  ## Need rewrite
  if (flip) {
    plt <- plt +    
      geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
      theme_minimal() + 
      theme(axis.text = element_text(size = textSize))+coord_flip()
  } else {
    plt <- plt +    
      geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
      theme_minimal() + 
      theme(axis.text = element_text(size = textSize, angle=90))
  }
  if (retList){
    return(list(wc=wc, plot=plt))
  } else {
    return(plt)
  }
}


#' exportCyjs
#' 
#' Export Cytoscape.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @param bubble default to FALSE, if node attribute "group" is available,
#' use the grouping to draw bubblesets using
#' `https://github.com/upsetjs/cytoscape.js-bubblesets`
#' @param withEdge connect edges or not in bubblesets
#' @import jsonlite
#' @importFrom cyjShiny dataFramesToJSON
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \dontrun{exportCyjs(g, "./", "net")}
#' @export
#' 
exportCyjs <- function(g, rootDir, netDir,
  bubble=FALSE, withEdge=TRUE) {
    
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    
    nodes <- data.frame(
        id=names(V(g)),
        label=names(V(g)),
        image=V(g)$image,
        size=V(g)$size,
        shape=V(g)$shape
    )
    
    if (!is.null(V(g)$group)) {
      nodes$group <- V(g)$group
    }

    edgeList <- as_edgelist(g)
    edges <- data.frame(source=edgeList[,1],
        target=edgeList[,2], interaction=NA)
    edges$strength <- E(g)$strength
    
    pret <- prettify(dataFramesToJSON(edges, nodes))
    pret <- substr(pret, 18, nchar(pret)-3)
    pret
    

    if (bubble) {
      js <- 'import cytoscapeBubblesets from "https://cdn.skypack.dev/cytoscape-bubblesets@3.1.0";'
    } else {
      js <- ""
    }

    js <- paste0(js, "
        var cy = window.cy = cytoscape({
            container: document.getElementById('cy'),
              style: cytoscape.stylesheet()
                .selector('node')
                .css({
                          'content': 'data(label)',
                          'shape' : 'data(shape)',
                          'background-image': 'data(image)',
                          'text-valign': 'bottom',
                          'background-color': '#FFF',
                          'background-fit': 'cover',
                          'width': 'data(size)',
                          'height': 'data(size)',
                          'font-size' : 'mapData(size, 0, 100, 1, 20)',
                          'text-outline-width': 1,
                          'text-outline-color': '#FFF',
                          'border-color' : '#555',
                          'border-width': 1
                      })
                .selector('edge')
                .css({
                        'width' : '4',
                        'target-arrow-shape': 'triangle',
                        'curve-style': 'bezier',
                        'width' : 'mapData(strength, 0.5, 1, 0, 5)'
                      }),
            'elements':
        ", pret, ",
            layout:{
                  name: 'cola',
                  padding: 0.5,
                  avoidOverlap: true, 
                  nodeSpacing: function( node ){ return 0.1; },
                  nodeDimensionsIncludeLabels: true
              }
            });
        ")

    if (bubble) {
      if (withEdge) {
        js <- paste0(js, 
               "cy.ready(() => {
                 let groups = cy.nodes().map(nodes => nodes.data().group).flat().filter(Boolean);
                 let uniqGroups = Array.from(new Set(groups));
                 const bb = cy.bubbleSets();
                 
                 for (let i = 0; i < uniqGroups.length; i++) {
                   var tmpgr = uniqGroups[i];
                   const nodes = cy.nodes().filter(function(nodes){
                     let gr = nodes.data().group;
                     if (gr!=null) {
                       let int = gr.includes(tmpgr);
                       return(int);
                     }
                   });
                   let gr = nodes.data().group;
                   console.log(nodes.length)
                   const edges = nodes.connectedEdges();
                   bb.addPath(nodes, edges);
                 }
               });"
        )
      } else {
        js <- paste0(js, 
                     "cy.ready(() => {
                 let groups = cy.nodes().map(nodes => nodes.data().group).flat().filter(Boolean);
                 let uniqGroups = Array.from(new Set(groups));
                 const bb = cy.bubbleSets();
                 
                 for (let i = 0; i < uniqGroups.length; i++) {
                   var tmpgr = uniqGroups[i];
                   const nodes = cy.nodes().filter(function(nodes){
                     let gr = nodes.data().group;
                     if (gr!=null) {
                       let int = gr.includes(tmpgr);
                       return(int);
                     }
                   });
                   let gr = nodes.data().group;
                   console.log(nodes.length)
                   const edges = nodes.connectedEdges();
                   bb.addPath(nodes);
                 }
               });"
        )
      }
    }

    ## Using cola layout by default.
    if (bubble) {
      html <- '
        <!DOCTYPE html>
        <html lang="en">
        
        <head>
            <meta charset="UTF-8">
            <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
            <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
            <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
            <link rel="stylesheet" type="text/css" href="style.css">
        </head>
        
        <body>
            <div id="cy"></div>
            <script type="module" src="script.js"></script>
        </body>
        
        </html>
        '
      
    } else {
      html <- '
        <!DOCTYPE html>
        <html lang="en">
        
        <head>
            <meta charset="UTF-8">
            <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
            <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
            <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
            <link rel="stylesheet" type="text/css" href="style.css">
        </head>
        
        <body>
            <div id="cy"></div>
            <script src="script.js"></script>
        </body>
        
        </html>
        '
    }
    style <- "
      body {
          font-family: helvetica, sans-serif;
          font-size: 14px;
      }
      
      #cy {
          position: absolute;
          left: 0;
          top: 0;
          right: 0;
          bottom: 0;
          z-index: 999;
      }
      
      h1 {
          opacity: 0.5;
          font-size: 1em;
      }"


    write(js, file = paste0(rootDir, "/", netDir, "/script.js"))
    write(style, file = paste0(rootDir, "/", netDir, "/style.css"))
    write(html, file = paste0(rootDir, "/", netDir, "/index.html"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' exportVisjs
#' 
#' Export vis.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @import jsonlite
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \dontrun{exportVisjs(g, "./", "net")}
#' @export
#' 
exportVisjs <- function(g, rootDir, netDir){
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    if (unique(V(g)$shape)=="rectangle"){
        visjsShape <- "image"
    } else {
        visjsShape <- "circularImage"
    }
    
    nodejson <- toJSON(data.frame(
            id=names(V(g)),
            label=names(V(g)),
            image=V(g)$image,
            shape=visjsShape,
            size=V(g)$size
        ))
    
    edgeList <- as_edgelist(g)
    edgejson <- toJSON(data.frame(from=edgeList[,1], to=edgeList[,2],
        width=E(g)$strength))
    
    # Make JS
    js <- paste0("
    var nodes = null;
    var edges = null;
    var network = null;
    
    function draw() {
    nodes = ", nodejson, ";
    edges = ", edgejson, ";
      var container = document.getElementById('mynetwork');
      var data = {
        nodes: nodes,
        edges: edges,
      };
      var options = {
        nodes: {
          borderWidth: 2,
          size: 30,
          color: {
            border: '#222222',
    background: 'white',
    },
    font: { color: 'black' },
    },
    edges: {
        length: 200,
        color: 'lightgray',
        arrows: { to: {enabled: true} }
    },
    layout: {
        improvedLayout: true
    }
    };
    network = new vis.Network(container, data, options);
    }
    
    window.addEventListener('load', () => {
        draw();
    });
    ")
    
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="mynetwork"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- '
    body {
      font: 10pt arial;
    }
    #mynetwork {
      width: 1000px;
      height: 1000px;
      border: 1px solid lightgray;
      background-color: white;
    }
    '
    
    write(js, file = paste0(rootDir, "/",netDir, "/script.js"))
    write(html, file = paste0(rootDir,"/",netDir, "/index.html"))
    write(style, file = paste0(rootDir, "/",netDir, "/style.css"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' returnExample
#' 
#' return an example dataset used in the analysis
#' 
#' @import org.Hs.eg.db
#' @return return example MEs and colors
#' @examples returnExample()
#' @export
#' 
returnExample <- function() {
    ## Simulate WGCNA results (three modules)
    ccls <- c()
    for (i in c(1,2,3,4,5,6,7,8,9)){
        ccls <- c(ccls, paste0("CCL",i))
    }
    cxcls <- c()
    for (i in c(1,2,3,5,6,8,9,10,11,12,13,14,16)){
        cxcls <- c(cxcls, paste0("CXCL",i))
    }
    erccs <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")

    CCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=ccls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    CXCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=cxcls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    ERCCensg <- AnnotationDbi::select(org.Hs.eg.db, keys=erccs,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL

    CCLensg <- CCLensg[!is.na(CCLensg)]
    CXCLensg <- CXCLensg[!is.na(CXCLensg)]
    ERCCensg <- ERCCensg[!is.na(ERCCensg)]

    CCLcol <- rep(1, length(CCLensg))
    names(CCLcol) <- CCLensg
    CXCLcol <- rep(2, length(CXCLensg))
    names(CXCLcol) <- CXCLensg
    ERCCcol <- rep(3, length(ERCCensg))
    names(ERCCcol) <- ERCCensg
    modColors <- c(CCLcol, CXCLcol, ERCCcol)
    ensg <- names(modColors)

    MEs <- data.frame(
        ME1 = 1:10,
        ME2 = 2:11,
        ME3 = 10:1
    )

    mod <- list()
    mod[["MEs"]] <- MEs
    mod[["colors"]] <- modColors
    mod
}

#' makeCorpus
#' 
#' Clean-up the corpus
#' 
#' @param docs corpus to clean
#' @param filterWords words to filter based on frequency
#' @param additionalRemove words to filter
#' @param numOnly delete number only
#' @param stem use stem or not
#' @param lower transform lower for corpus
#' 
#' @return cleaned corpus
#' @import tm
#' 
#' 
makeCorpus <- function (docs, filterWords, additionalRemove,
    numOnly, stem, lower=TRUE) {
    if (lower) {
        docs <- docs %>%
            tm_map(FUN=content_transformer(tolower))
    }
    if (numOnly) {
        docs <- docs %>% tm_map(FUN=removeAloneNumbers)
    } else {
        docs <- docs %>% tm_map(FUN=removeNumbers)
    }
    docs <- docs %>%
        tm_map(removeWords, stopwords::stopwords("english",
            "stopwords-iso")) %>%
        tm_map(removeWords, filterWords) %>% 
        tm_map(FUN=removePunctuation) %>%
        tm_map(FUN=stripWhitespace)
    if (prod(is.na(additionalRemove))!=1){
        docs <- docs %>% tm_map(removeWords, additionalRemove)
    }
    if (stem) {
        docs <- docs %>% tm_map(stemDocument)   
    }
    return(docs)
}


#' getUPtax
#' 
#' Obtain the list of UniProt organism identification codes, by querying taxonomy or UniProt codes.
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt
#' The file above must be downloaded and specified to file argument.
#' 
#' @param file downloaded file
#' @param candUP candidate UniProt organism identification codes
#' @param candTax candidate taxonomy name
#' @return data.frame consisting of taxonomy name and UniProt IDs
#' @examples
#' file <- "speclist.txt"
#' \dontrun{getUPtax(file, candUP="all")}
#' @export
#' 
getUPtax <- function(file, candUP, candTax=NULL) {
  tmp <- NULL
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (grepl("N=",line)) {
      code <- strsplit(line, " ")[[1]][1]
      if (code=="Code") {next}
      sn <- unlist(strsplit(line, "N="))[2]
      if (length(candUP)==1) {
        if (candUP=="all") {
          if (!is.null(candTax)) {
            for (ct in candTax) {
              if (grepl(ct, sn)){
                tmp <- rbind(tmp, c(code, sn))
              }
            }
          } else {
            tmp <- rbind(tmp, c(code, sn))
          }
        }
      }
      if (code %in% candUP) {
        tmp <- rbind(tmp, c(code, sn))
      }
    }
  }
  close(con)
  if (!is.null(tmp)) {
    tmp <- data.frame(tmp) |> `colnames<-`(c("UPID","Taxonomy"))
    return(tmp)
  } else {
    return(NULL)
  }
}




#' exportCyjsWithoutImage
#' 
#' Export Cytoscape.js script, HTML and stylesheet for the graph without image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @param edgeWidth attribute name for edgeWidth
#' @param nodeColor attribute name for node color
#' @param nodeSize attribute name for node size
#' @param nodeColorDiscretePal color mapping palette for discrete variables
#' @param sizeMin minimum size for scaling node size
#' @param sizeMax maximum size for scaling node size
#' @import jsonlite
#' @importFrom cyjShiny dataFramesToJSON
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$size <- c(1,1)
#' \dontrun{exportCyjsWithoutImage(g, "./", "net")}
#' @export
#' 
exportCyjsWithoutImage <- function(g, rootDir, netDir,
  edgeWidth="weight",nodeColor="tag",nodeSize="Freq",nodeColorDiscretePal="RdBu",
  sizeMin=10, sizeMax=50) {
    
    if (is.null(V(g)$shape)){qqcat("No node shape specified, set to 'circle'\n")
      V(g)$shape <- rep("circle", length(V(g)))}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}


    ## Make node color
    nc <- get.vertex.attribute(g, nodeColor)
    nodeColors <- RColorBrewer::brewer.pal(length(unique(nc)), nodeColorDiscretePal)
    names(nodeColors) <- unique(nc)

    ## Make node size
    V(g)$size <- get.vertex.attribute(g, nodeSize)
    rawMin <- min(V(g)$size)
    rawMax <- max(V(g)$size)
    scf <- (sizeMax-sizeMin)/(rawMax-rawMin)
    V(g)$size <- scf * V(g)$size + sizeMin - scf * rawMin

    nodes <- data.frame(
        id=names(V(g)),
        label=names(V(g)),
        size=V(g)$size,
        color=nodeColors[nc],
        shape=V(g)$shape
    )
    
    edgeList <- as_data_frame(g)
    edgeListRename <- colnames(edgeList)
    edgeListRename[1] <- "source";edgeListRename[2] <- "target"
    colnames(edgeList) <- edgeListRename
    edgeList$interaction <- rep(NA, nrow(edgeList))
    edgeList$width <- edgeList[[edgeWidth]]
    # edges <- data.frame(source=edgeList[,1],
    #     target=edgeList[,2], interaction=NA)
    # edges$strength <- E(g)$strength
    
    pret <- prettify(dataFramesToJSON(edgeList, nodes))
    pret <- substr(pret, 18, nchar(pret)-3)
    pret
    
    js <- paste0("
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
          style: cytoscape.stylesheet()
            .selector('node')
            .css({
                      'content': 'data(label)',
                      'shape' : 'data(shape)',
                      'text-valign': 'bottom',
                      'background-color': 'data(color)',
                      'background-fit': 'cover',
                      'width': 'data(size)',
                      'height': 'data(size)',
                      'font-size' : 'mapData(size, 0, 100, 10, 20)',
                      'text-outline-width': 1,
                      'text-outline-color': '#FFF',
                      'border-color' : '#555',
                      'border-width': 1
                  })
            .selector('edge')
            .css({
                    'curve-style': 'bezier',
                    'width' : 'mapData(width, 0.5, 1, 1, 5)'
                  }),
        'elements':
    ", pret, ",
        layout:{
              name: 'cola',
              padding: 0.5,
              avoidOverlap: true, 
              nodeSpacing: function( node ){ return 0.1; },
              nodeDimensionsIncludeLabels: true
          }
        });
    ")
    
    ## Using cola layout by default.
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="cy"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- "
    body {
        font-family: helvetica, sans-serif;
        font-size: 14px;
    }
    
    #cy {
        position: absolute;
        left: 0;
        top: 0;
        right: 0;
        bottom: 0;
        z-index: 999;
    }
    
    h1 {
        opacity: 0.5;
        font-size: 1em;
    }"
    
    
    write(js, file = paste0(rootDir, "/",netDir, "/script.js"))
    write(style, file = paste0(rootDir, "/",netDir, "/style.css"))
    write(html, file = paste0(rootDir, "/",netDir, "/index.html"))
    # message(paste0("Exported to ",rootDir,netDir))
}


#' obtainTextPosition
#' 
#' obtain text position in biofabric layout
#' 
#' @return tbl_graph
#' @param ig igraph
#' @param sort.by the argument to be passed to fabric layout function
#' @param verbose show logs
#' @examples refseq(c("PNKP","DDX41"))@igraph |> obtainTextPosition()
#' @export
#' 
obtainTextPosition <- function(ig, sort.by=node_rank_fabric(),
  verbose=FALSE) {
  fab <- ggraph(ig,
       "fabric",
        sort.by = sort.by)

  textPos <- fab$data
  sortedTextPos <- textPos[order(textPos$x),]

  # Obtain raw position index
  rawpos <- NULL
  for (i in seq_len(nrow(sortedTextPos))) {
      cur <- sortedTextPos$name[i]
      nex <- sortedTextPos$name[i+1]
      tmpCur <- textPos[textPos$name==cur,]
      tmpNex <- textPos[textPos$name==nex,]
      max1 <- tmpCur$xmax |> as.numeric()
      max2 <- tmpNex$xmax |> as.numeric()
      if (sum(is.na(max2))==0) {
        if (max2 < max1) {
          rawpos <- c(rawpos, tmpNex$.ggraph.orig_index)
        }
      }
  }

  rawpos <- NULL
  cur_max <- 0
  flag <- NULL
  adjFlag <- NULL
  for (i in seq_len(nrow(sortedTextPos))) {
      cur <- sortedTextPos$name[i]
      tmpCur <- textPos[textPos$name==cur,]
      tmp_max <- tmpCur$xmax |> as.numeric()
      if (tmp_max < cur_max) {
        flag <- c(flag, i)
      } else {
        cur_max <- tmp_max
        adjFlag <- c(adjFlag, i)
      }
  }

  raws <- sortedTextPos[flag,]
  raws$center <- raws$xmax

  adjs <- sortedTextPos[adjFlag,]

  if (verbose) {
    qqcat("Raw: @{dim(raws)[1]}, Position to be adjusted: @{dim(adjs)[1]}")
  }

  nextpos <- NULL
  for (i in seq_len(nrow(adjs))) {
      cur <- adjs$name[i]
      nex <- adjs$name[i+1]
      tmpCur <- textPos[textPos$name==cur,]
      tmpNex <- textPos[textPos$name==nex,]
      max1 <- tmpCur$xmax |> as.numeric()
      max2 <- tmpNex$xmax |> as.numeric()
      nextpos <- c(nextpos, max1 + (max2 - max1)/2)
  }
  nextpos <- nextpos[!is.na(nextpos)]
  nextpos <- c(adjs$x[1] |> as.numeric(), nextpos)
  adjs$center <- nextpos
  adjs <- rbind(adjs, raws)
  adjs <- adjs[as.character(seq_len(nrow(adjs))),]
  ret <- ig |> as_tbl_graph() |> mutate(center=adjs$center)
  ret
}

#' plot_biofabric
#' 
#' plot the network in biofabric layout
#' Internally, obtainTextPosition() is used and 
#' the layers are stacked.
#' For BioFabric layout, please refer to: https://biofabric.systemsbiology.net/
#' Citation is: https://doi.org/10.1186/1471-2105-13-275
#' 
#' @param res biotext object
#' @param sort.by passed to fabric layout
#' @param highlight which category to highlight, default to NULL
#' @param highlight_color highlight color of the category
#' @param size mapping to use for size of text on x-axis
#' default to "rank", which uses default ordering,
#' other accepted options are those in the node table of tbl_graph object
#' @param textScaleRange decide text scale range on x-axis
#' @param end_shape end shape on `geom_edge_span`
#' @param color_map passed to aes mapping in geom_node_range.
#' Override `highlight` option.
#' @return ggplot2
#' @export
#' @examples refseq(c("DDX41","PNKP")) |> plot_biofabric()
#'
#'
plot_biofabric <- function(res,
                           sort.by = node_rank_fabric(),
                           size="rank",
                           end_shape="circle",
                           color_map=NULL,
                           highlight=NULL,
                           highlight_color="tomato",
                           textScaleRange=c(1.5,2)) {
    if (size=="rank") {size <- "xmin"}
    tbl <- res@igraph |> obtainTextPosition(sort.by=sort.by)
    if (!is.null(color_map)) {
        g <- ggraph(tbl, "fabric", sort.by = sort.by)+
            geom_node_range(aes(color=!!sym(color_map))) +
            geom_edge_span(end_shape = end_shape) +
            geom_node_shadowtext(aes(x=.data$xmin-4, label=.data$name), color="grey20",size=2, bg.colour="white")+
            geom_node_shadowtext(aes(x=.data$center, size=!!sym(size), y=.data$y+1, label=.data$name),
                            bg.colour="white", color="grey20")+
            theme_graph()+
            scale_size(trans = 'reverse', range=textScaleRange)+
            guides(size="none")
        return(g)
    }
    if (!is.null(highlight)) {
      g <- ggraph(tbl, "fabric", sort.by = sort.by)+
        geom_node_range(aes(color=.data$nodeCat)) +
        geom_edge_span(end_shape = end_shape) +
        geom_node_shadowtext(aes(x=.data$xmin-4, label=.data$name, filter=.data$type==highlight),
                             size=2, color=highlight_color, bg.colour="white")+
        geom_node_text(aes(x=.data$xmin-4, label=.data$name, filter=.data$type!=highlight), size=2)+
        geom_node_shadowtext(aes(x=.data$center, size=!!sym(size),
          filter=.data$type!=highlight, color=.data$nodeCat, y=.data$y+1, label=.data$name), bg.colour="white")+
        geom_node_shadowtext(aes(x=.data$center, size=!!sym(size),
          filter=.data$type==highlight, color=.data$nodeCat,
                                 y=.data$y+1, label=.data$name), bg.colour="white")+
        scale_color_manual(values=c(highlight_color, "grey20"),name="Category")+
        theme_graph()+
        scale_size(trans = 'reverse', range=textScaleRange)+
        guides(size="none")
    } else {
      g <- ggraph(tbl, "fabric", sort.by = sort.by)+
        geom_node_range() +
        geom_edge_span(end_shape = end_shape) +
        geom_node_shadowtext(aes(x=.data$xmin-4, label=.data$name), color="grey20",size=2, bg.colour="white")+
        geom_node_shadowtext(aes(x=.data$center, size=!!sym(size), y=.data$y+1, label=.data$name),
                             bg.colour="white", color="grey20")+
        theme_graph()+
        scale_size(trans = 'reverse', range=textScaleRange)+
        guides(size="none")
    }
    g
}