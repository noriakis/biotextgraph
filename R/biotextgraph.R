#' biotextgraph
#' 
#' wrapper for functions wcGeneSummary, wcAbst, and wcBSDB
#' 
#' @param target "pubmed", "bugsigdb", "refseq", "ec"
#' @param argList passed to each function
#' @return list of data frames and ggplot2 object
#' 
#' @examples
#' geneList <- c("DDX41","PNKP")
#' biotextgraph("refseq", argList=list(geneList=geneList))
#' @export
#' 
biotextgraph <- function(target, argList) {
	if (target=="pubmed"){
		return(do.call("wcAbst", argList))
	} else if (target=="bugsigdb"){
		return(do.call("wcBSDB", argList))
	} else if (target=="refseq"){
		return(do.call("wcGeneSummary", argList))
	} else if (target=="ec") {
		return(do.call("wcEC", argList))
	} else {
		stop("Please specify pubmed, bugsigdbr, or refseq")
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
setClass("biotext", slots=list(
        query="character",
        delim="character",
        type="character",
        model="character",
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
        igraph="igraph",
        geneCount="table",
        geneMap="matrix",
        net="ggraph",
        wc="gg",
        ec="data.frame",
        wholeFreq="numeric",
        dic="vector"
        ))

#' @importFrom utils object.size
setMethod("show",
  signature(object="biotext"),
  function(object) {
    qqcat("Type: @{object@type}\n")
    qqcat("Number of words: @{object@numWords}\n")
    if (length(object@query)<10) {
      cat(paste(object@query, collapse="/"));cat("\n")
    } else {
      cat(paste0(paste(object@query[1:10],
        collapse="/"), "/truncated"));cat("\n")
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
              gradn <- adjustcolor(pal(unique(length(V(g)$tag)))[as.numeric(factor(V(g)$tag))], 0.8)
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


#' plotORA
#' 
#' plot volcano-plot like plot for ORA results
#' 
#' @param x biotext object
#' @param thresh hline to draw in plot
#' @examples
#' net1 <- wcGeneSummary(c("DDX41","IRF3"), ora=TRUE)
#' plotORA(net1)
#' @return volcano plot for ORA results
#' @importFrom ggrepel geom_text_repel
#' @export
#' 
plotORA <- function(x, thresh=0.001) {
    subr <- intersect(row.names(x@freqDf),names(x@ora))
  vp <- x@freqDf[subr, ]
  vp$p <- -log10(x@ora[subr])
  
  ggplot(vp, aes(x=vp$freq,y=vp$p, fill=vp$p))+
    geom_point(shape=21,size=3,show.legend=FALSE) +
    geom_text_repel(aes(color=vp$p, label=vp$word),bg.color = "white",
                    segment.color="black",size=3,max.overlaps = Inf,
                    bg.r = .15, show.legend=FALSE)+
    scale_color_gradient(low="blue",high="red")+
    scale_fill_gradient(low="blue",high="red")+
    geom_hline(yintercept=-log10(thresh), linetype="dashed")+
    xlab("Frequency")+ylab("-log10(p.adjust)")+
    theme_minimal()
}

