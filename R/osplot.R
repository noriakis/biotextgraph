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
        filtered="character",
        excludeFreq="numeric",
        excludeTfIdf="numeric",
        enrichResults="data.frame",
        rawText="data.frame",
        TDM="TermDocumentMatrix",
        corpus="VCorpus",
        freqDf="data.frame",
        pvclust="list",
        pvpick="list",
        corMat="matrix",
        graph="igraph",
        geneCount="vector",
        geneMap="matrix",
        net="list"
        ))
setMethod("show",
  signature(object="osplot"),
  function(object) {
    qqcat("Query: @{object@query}\n")
    print(object.size(object), units="auto")
  })