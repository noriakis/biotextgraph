#' biotextgraph
#' 
#' @description wrapper for functions refseq, pubmed, enzyme, and bugsigdb
#' @details The function calls the various types of databases (refseq, pubmed, ...)
#' for summarizing the textual data.
#' 
#' @param target "pubmed", "bugsigdb", "refseq", "ec"
#' @param argList passed to each function
#' @return list of data frames and ggplot2 object
#' 
#' @examples
#' geneList <- c("DDX41","PNKP","IRF3","XRCC1","ERCC2")
#' biotextgraph("refseq", argList=list(geneList=geneList))
#' @export
#' 
biotextgraph <- function(target, argList) {
	if (target=="pubmed"){
		return(do.call("pubmed", argList))
	} else if (target=="bugsigdb"){
		return(do.call("bugsigdb", argList))
	} else if (target=="refseq"){
		return(do.call("refseq", argList))
	} else if (target=="ec") {
		return(do.call("enzyme", argList))
	} else {
		stop("Please specify pubmed, bugsigdbr, ec, or refseq")
	}
}

setOldClass("pvclust")
setOldClass("igraph")
setOldClass("VCorpus")
setOldClass("corpus")
setOldClass("TermDocumentMatrix")
setOldClass("gg")
setOldClass("ggraph")
setOldClass("dfm")
setOldClass("communities")
setClass("biotext", slots=list(
	date="POSIXct",
    query="character",
    delim="character",
    type="character",
    model="character",
    tag="character",
    filtered="character",
    pmids="character",
    retMax="numeric",
    excludeFreq="numeric",
    excludeTfIdf="numeric",
    numWords="numeric",
    corThresh="numeric",
    ora="vector",
    enrichResults="data.frame",
    rawText="data.frame",
    rawAnnot="data.frame",
    rawTextBSDB="data.frame",
    TDM="TermDocumentMatrix",
    dfm="dfm",
    corpus="VCorpus",
    corpusQuanteda="corpus",
    freqDf="data.frame",
    pvclust="pvclust",
    pvpick="list",
    strength="data.frame",
    corMat="matrix",
    igraphRaw="igraph",
    igraph="igraph",
    geneCount="table",
    geneMap="matrix",
    net="ggraph",
    wc="gg",
    ec="data.frame",
    wholeFreq="numeric",
    dic="vector",
    sortOrder="character",
    numOnly="logical",
    stem="logical",
    ngram="numeric",
    curate="logical",
    communities="communities"
))

#' @importFrom utils object.size
setMethod("show",
  signature(object="biotext"),
  function(object) {
    qqcat("Type: @{object@type}\n")
    qqcat("Number of words: @{object@numWords}\n")
    if (length(object@query)<10) {
      cat(paste0("Query: ",paste(object@query, collapse="/")));cat("\n")
    } else {
      cat(paste0("Query: ",paste0(paste(object@query[seq_len(10)],
        collapse="/"), "/truncated")));cat("\n")
    }
    deg <- NULL; vnum <- NULL; enum <- NULL;
    if (is.igraph(object@igraphRaw)) {
      deg <- degree(object@igraphRaw)
      vnum <- length(V(object@igraphRaw)); enum <- length(E(object@igraphRaw))
      ord <- V(object@igraphRaw)$name[order(deg, decreasing=TRUE)]
    }
    if (is.igraph(object@igraph)) {
      deg <- degree(object@igraph)
      vnum <- length(V(object@igraph)); enum <- length(E(object@igraph))
      ord <- V(object@igraph)$name[order(deg, decreasing=TRUE)]
    }
    if (!is.null(deg)) {
      showdeg <- paste0(paste0(ord[1:5],
        "(",deg[order(deg, decreasing=TRUE)][1:5],")"),
        collapse="/")
      qqcat("Graph: V(@{vnum}), E(@{enum})\n")
      qqcat("Degree: @{showdeg}\n")
    }
    print(object.size(object), units="auto")
  })


#' @importFrom grDevices adjustcolor
setMethod("plot",
          signature = "biotext",
          definition = function(x) {
            retSc <- function(va, min=4,max=9){
              (max-min) * ((va-min(va)) / 
                             (max(va)-min(va))) + min
              
            }
            g <- x@igraph
            
            fillna <- V(g)$Freq
            fillna[is.na(fillna)] <- min(fillna[!is.na(fillna)])
            V(g)$Freq <- fillna
            
            if (length(x@pvpick)!=0) {
              pal <- colorRampPalette(brewer.pal(8,"Set2"))
              gradn <- adjustcolor(
                pal(unique(length(V(g)$tag)))[as.numeric(factor(V(g)$tag))], 0.8
                )
            } else {
              pal <- colorRampPalette(c("blue","red"))
              gradn <- adjustcolor(pal(length(V(g)))[V(g)$Freq],0.5)
            }
            vs <- retSc(V(g)$Freq, 4,9)
            tsz <- retSc(V(g)$Freq, 1,2)
            
            plot(g,
                 vertex.color=gradn,
                 # vertex.label.color=gradn,
                 vertex.size=vs,
                 vertex.label.cex=tsz,
                 vertex.label.dist=1,
                 vertex.label.family="arial",
                 edge.curved=0)
          })

#' plotNet
#' 
#' @description Plot the network stored in the biotext object with changing the visualization parameters.
#' @details The function accepts the already calculated biotext object and outputs the visualization based on
#' the specified parameters.
#' 
#' @rdname plotnet
#' @param x biotextgraph object
#' @param layout the layout for the network, defaul to "nicely".
#' It can be one of the layouts implemented in `igraph` and `ggraph`, such as
#' `kk` (Kamada-Kawai), `nicely` (automatic selection of algorithm), `drl` (the force-directed DrL layout).
#' The options are available at: https://igraph.org/r/doc/layout_.html
#' @param edgeLink if FALSE, use geom_edge_diagonal. if TRUE, geom_edge_link. Default to TRUE.
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param pal palette for color gradient in correlation network.
#' should be a vector of length two like c("red","blue").
#' @param showLegend whether to show legend in the network
#' @param colorText color text label based on frequency in the network
#' @param tagPalette tag palette when `tag` is TRUE. It is also used for dependency network
#' using udpipe, and tagging colorization for word cloud.
#' Default to NULL, which indicates automatically set.
#' @param naEdgeColor edge colors for NA values (linking query with the category other than text)
#' @param fontFamily font family to use, default to "sans".
#' @param colorize color the word nodes by their frequency, and the other nodes by their category.
#' if colorize=FALSE and addFreqToGene=TRUE, gene nodes are colorized according to the minimum frequency 
#' of the words in the network
#' @param discreteColorWord colorize words by "Words" category, not frequency.
#' @param catColors colors for words and texts when colorize is TRUE and discreteColorWord is TRUE
#' @param scaleEdgeWidth scale for edge width
#' @param asis plot the original network (default to FALSE)
#' @param queryColor color for associated queries with words
#' @param useSeed random seed
#' @param scaleRange scale for label and node size in the network.
#' @param autoScale scale the label and node size automatically for the large network.
#' @export
#' @examples
#' library(ggraph)
#' geneList <- c("DDX41","PNKP","ERCC1","IRF3","XRCC1")
#' test <- refseq(geneList)
#' plotNet(test, asis=TRUE)
#' @return biotext object with network visualization changed
setGeneric("plotNet",
    function(x, layout="nicely", edgeLink=TRUE,
    	edgeLabel=FALSE, showLegend=FALSE, fontFamily="sans",
    	tagPalette=NULL, catColors=NULL, queryColor="grey",
    	pal=c("blue","red"), colorize=FALSE,
    	discreteColorWord=FALSE, useSeed=42, autoScale=FALSE,
    	scaleRange=c(5,10), scaleEdgeWidth=c(1,3),
    	naEdgeColor="grey", colorText=FALSE, asis=FALSE)
    standardGeneric("plotNet"))
#' @rdname plotnet
setMethod("plotNet", "biotext",
    function(x, layout="nicely", edgeLink=FALSE,
    	edgeLabel=FALSE, showLegend=FALSE, fontFamily="sans",
    	tagPalette=NULL, catColors=NULL, queryColor="grey",
    	pal=c("blue","red"), colorize=FALSE,
    	discreteColorWord=FALSE, useSeed=42, autoScale=FALSE,
    	scaleRange=c(5,10), scaleEdgeWidth=c(1,3),
    	naEdgeColor="grey", colorText=FALSE, asis=FALSE) {
    	
    	if (x@type=="combine") {
    		asis <- TRUE
    	}
        if (x@type=="udpipe") {
            asis <- TRUE
        }
    	if (asis) {
    		return(x@net)
    	}
    	    	
		allnodecat <- V(x@igraph)$nodeCat
		allnodecat <- allnodecat[allnodecat!="Words"] |> unique()
    	coGraph <- x@igraph

        netPlot <- ggraph(coGraph, layout=layout)

        netPlot <- appendEdges(netPlot, FALSE, edgeLink,
            edgeLabel, showLegend, fontFamily)


        if (!is.null(names(x@pvpick))) { ## Obtain tag coloring
            if (is.null(tagPalette)) {
                cols <- V(coGraph)$tag |> unique()
                if (length(cols)>2) {
                    tagPalette <- RColorBrewer::brewer.pal(8, "Dark2")
                    tagPalette <- colorRampPalette(tagPalette)(length(cols))
                } else {
                    tagPalette <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(length(cols))]
                }
                names(tagPalette) <- cols
                tagPalette[allnodecat] <- queryColor
            }
        }

        if (is.null(catColors)) { ## Obtain category coloring
            catLen <- length(unique(V(coGraph)$nodeCat))
            if (catLen>2) {
                catColors <- RColorBrewer::brewer.pal(catLen, "Dark2")
            } else {
                catColors <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(catLen)]
            }
            names(catColors) <- unique(V(coGraph)$nodeCat)
            catColors[allnodecat] <- queryColor
        }
        if ("tag" %in% slotNames(x)) {
            tag <- x@tag        
        } else {
            tag <- "none"
        }
        if (identical(tag, character(0))) {tag <- "none"}

        netPlot <- appendNodesAndTexts(netPlot, tag, colorize, tagPalette,
                          showLegend, catColors, pal, fontFamily, colorText,
                          scaleRange, useSeed, ret=x, tagColors=tagPalette,
                          discreteColorWord=discreteColorWord)
        if (autoScale) {
        	scaleRange <- c((500 * (1 / x@numWords))/2.5,
        		500 * (1 / x@numWords))
        }
        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=scaleEdgeWidth, name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation", na.value=naEdgeColor)+
            theme_graph()

        if (dim(x@enrichResults)[1]!=0) {
            netPlot <- netPlot + ggforce::geom_mark_hull(
                aes(netPlot$data$x,
                    netPlot$data$y,
                    group = netPlot$data$grp,
                    label=netPlot$data$grp,
                    fill=netPlot$data$grp,
                    filter = !is.na(netPlot$data$grp)),
                concavity = 4,
                expand = unit(2, "mm"),
                alpha = 0.25,
                na.rm = TRUE,
                show.legend = FALSE
            )
        }
        netPlot
    }
)


#' plotWC
#' 
#' @description Plot the wordcloud changing the visualization parameters.
#' @details The function accepts the already calculated biotext object and outputs the visualization based on
#' the specified parameters.
#' @rdname plotwc
#' @param x biotext object
#' @param tagPalette tag palette when `tag` is TRUE. It is also used for dependency network
#' using udpipe, and tagging colorization for word cloud.
#' Default to NULL, which indicates automatically set.
#' @param madeUpper make these words uppercase in resulting plot,
#' default to c("rna" and "dna")
#' @param preserve Try to preserve original characters.
#' @param scaleFreq default to NULL, scale the value if specified
#' @param numWords the number of words to be shown in the plot.
#' When `autoThresh` is TRUE, the number of this value will be shown.
#' @param useggwordcloud default to TRUE, otherwise use `wordcloud` function.
#' @param wcScale scaling size for ggwordcloud
#' @param argList parameters to pass to wordcloud() or ggwordcloud()
#' @param asis plot the original network (default to FALSE)
#' @param fontFamily font family to use, default to "sans".
#' @export
#' @examples
#' geneList <- c("DDX41","PNKP","ERCC1","IRF3","XRCC1")
#' test <- refseq(geneList, plotType="wc")
#' plotWC(test, asis=TRUE)
#' @return wordcloud visualization
setGeneric("plotWC",
    function(x, tagPalette=NULL, madeUpper=c("dna","rna"),
    	preserve=FALSE, scaleFreq=NULL, fontFamily="sans", numWords=NULL,
    	wcScale=10, argList=list(), useggwordcloud=TRUE, asis=FALSE)
    standardGeneric("plotWC"))
#' @rdname plotwc
setMethod("plotWC", "biotext",
    function(x, tagPalette=NULL, madeUpper=c("dna","rna"),
    	preserve=FALSE, scaleFreq=NULL, fontFamily="sans", numWords=NULL,
    	wcScale=10, argList=list(), useggwordcloud=TRUE, asis=FALSE) {
    	
    	if (asis) {
    		return(x@wc)
    	}
    	if (is.null(numWords)) {
    		numWords <- x@numWords	
    	}
	    matSorted <- x@wholeFreq
	    if (length(matSorted) < numWords) {
	    	numWords <- length(matSorted)
	    }
	    matSorted <- matSorted[1:numWords]
	    
        if ("tag" %in% slotNames(x)) {
            tag <- x@tag        
        } else {
            tag <- "none"
        }
        if (identical(tag, character(0))) {tag <- "none"}

	    docs <- x@TDM
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
            
        wcCol <- NULL
        returnDf$wcColor <- "black"
        
        genePlot <- FALSE
        if (dim(x@geneMap)[1]!=0) {
        	genePlot <- TRUE
        	genemap <- x@geneMap
            genemap <- data.frame(genemap) |> `colnames<-`(c("word","gene"))
            collapsed_genemap <- genemap %>%
                group_by(.data$word) %>%
                summarise(gene_name=paste0(.data$gene, collapse=","))
            returnDf <- merge(returnDf, collapsed_genemap, by="word")
        }

        
        if (!is.null(names(x@pvpick))) {

        	pvc <- x@pvclust
        	pvcl <- x@pvpick
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
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        if (preserve) {
        	pdic <- x@dic
            for (nm in unique(returnDf$word)) {
                if (nm %in% names(pdic)) {
                    returnDf[returnDf$word == nm, "word"] <- pdic[nm]
                }
            }
        }
        
        if (!is.null(scaleFreq)) {
            showFreq <- returnDf$freq*scaleFreq
            returnDf$freq <- returnDf$freq*scaleFreq
        } else {
            showFreq <- returnDf$freq
        }
        
        if (!("min.freq" %in% names(argList))) {
            argList[["min.freq"]] <- 3
        }
        returnDf$wcColor <- wcCol
        returnDf <- returnDf[returnDf$freq > argList[["min.freq"]], ]

        if (tag!="none"){
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- returnDf$freq
            argList[["family"]] <- fontFamily
            argList[["colors"]] <- returnDf$wcColor
            argList[["random.order"]] <- FALSE
            argList[["ordered.colors"]] <- TRUE
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                if (genePlot) {
                    argList[["label_content"]] <- 
                    sprintf("%s<span style='font-size:7.5pt'><br>(%s)</span>",
                        returnDf$word, returnDf$gene_name)
                }
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill="transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        } else {
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- returnDf$freq
            argList[["family"]] <- fontFamily
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                if (genePlot) {
                    argList[["label_content"]] <- 
                    sprintf("%s<span style='font-size:7.5pt'><br>(%s)</span>",
                        returnDf$word, returnDf$gene_name)
                }
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill = "transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        }
        wc 	
    }
    
)


#' getSlot
#' 
#' @description get the slot value from biotext object
#' @details get the slot value from biotext object
#' 
#' @param x biotext object
#' @param slot slot name
#' @export
#' @examples
#' n1 <- refseq(c("IRF3","PNKP","DDX41","ERCC1","ERCC2","XRCC1")) 
#' getSlot(n1, "igraphRaw")
#' @return attribute value
setGeneric("getSlot",
    function(x, slot) standardGeneric("getSlot"))
#' getSlot
#' 
#' @description get the slot value from biotext object
#' @details get the slot value from biotext object
#' 
#' @param x biotext object
#' @param slot slot name
#' @export
#' @examples
#' n1 <- refseq(c("IRF3","PNKP","DDX41","ERCC1","ERCC2","XRCC1")) 
#' getSlot(n1, "igraphRaw")
#' @return attribute value
setMethod("getSlot", "biotext",
    function(x, slot) attr(x, slot))


#' plotORA
#' 
#' @description plot volcano-plot like plot for ORA results
#' @details Plot the volcano-plot like plot for the ORA results
#' using the biotext object. The ORA should be performed beforehand
#' by specifying ora option to TRUE.
#' 
#' @param x biotext object
#' @param thresh hline to draw in plot
#' @examples
#' testgenes <- c("DDX41","IRF3","XRCC1","ERCC1","ERCC2","ERCC3")
#' net1 <- refseq(testgenes, ora=TRUE)
#' plotORA(net1)
#' @return volcano plot for ORA results
#' @importFrom ggrepel geom_text_repel
#' @export
#' 
plotORA <- function(x, thresh=0.001) {
  subr <- intersect(tolower(x@freqDf$word), names(x@ora))
  vp <- x@freqDf[tolower(x@freqDf$word) %in% subr, ]
  vp$p <- -log10(x@ora[subr])
  
  ggplot(vp, aes(x=.data$freq,y=.data$p, fill=.data$p))+
    geom_point(shape=21,size=3,show.legend=FALSE) +
    geom_text_repel(aes(color=.data$p, label=.data$word),bg.color = "white",
                    segment.color="black",size=3,max.overlaps = Inf,
                    bg.r = .15, show.legend=FALSE)+
    scale_color_gradient(low="blue",high="red")+
    scale_fill_gradient(low="blue",high="red")+
    geom_hline(yintercept=-log10(thresh), linetype="dashed")+
    xlab("Frequency")+ylab("-log10(p.adjust)")+
    theme_minimal()
}

