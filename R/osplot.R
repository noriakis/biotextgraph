#' osplot
#' 
#' wrapper for functions wcGeneSummary, wcAbst, and wcBSDB
#' 
#' @param target "pubmed", "bugsigdb", "refseq", "ec"
#' @param argList passed to each function
#' @return list of data frames and ggplot2 object
#' 
#' @examples
#' geneList <- c("DDX41")
#' osplot("refseq", geneList)
#' @export
#' 
osplot <- function(target, argList) {
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
setOldClass("TermDocumentMatrix")
setOldClass("gg")
setOldClass("ggraph")
setClass("osplot", slots=list(
        query="character",
        delim="character",
        type="character",
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
        rawTextBSDB="data.frame",
        TDM="TermDocumentMatrix",
        corpus="VCorpus",
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
        oraPlot="gg",
        dic="vector"
        ))
setMethod("show",
  signature(object="osplot"),
  function(object) {
    qqcat("type: @{object@type}\n")
    qqcat("Number of words: @{object@numWords}\n")
    if (length(object@query)<10) {
      cat(paste(object@query, collapse="/"));cat("\n")
    } else {
      cat(paste0(paste(object@query[1:10],
        collapse="/"), "/truncated"));cat("\n")
    }
    print(object.size(object), units="auto")
  })


setMethod("plot",
          signature = "osplot",
          definition = function(x) {
            retSc <- function(x, min=4,max=9){
              (max-min) * ((x-min(x)) / 
                             (max(x)-min(x))) + min
              
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
                 vertex.label.color=gradn,
                 vertex.size=vs,
                 vertex.label.cex=tsz,
                 vertex.label.dist=1,
                 vertex.label.family="arial",
                 edge.curved=0)
          })
