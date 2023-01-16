#' osplot
#' 
#' wrapper for functions wcGeneSummary, wcAbst, and wcBSDB
#' 
#' @param target "pubmed", "bugsigdb", "refseq", "ec"
#' @param ... passed to each function
#' @return list of data frames and ggplot2 object
#' 
#' @examples
#' geneList <- c("DDX41")
#' osplot("refseq", geneList)
#' @export
#' 
osplot <- function(target, ...) {
	if (target=="pubmed"){
		return(wcAbst(...))
	} else if (target=="bugsigdb"){
		return(wcBSDB(...))
	} else if (target=="refseq"){
		return(wcGeneSummary(...))
	} else if (target=="ec") {
		return(wcEC(...))
	} else {
		stop("Please specify pubmed, bugsigdbr, or refseq")
	}
}

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
        rawTextBSDB="character",
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
        ec="data.frame"
        ))
setMethod("show",
  signature(object="osplot"),
  function(object) {
    qqcat("type: @{object@type}\n")
    qqcat("Number of words: @{object@numWords}\n")
    cat(paste(object@query, collapse="/"));cat("\n")
    print(object.size(object), units="auto")
  })