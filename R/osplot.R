#' osplot
#' 
#' wrapper for functions wcGeneSummary, wcAbst, and wcBSDB
#' 
#' @param target "abstract", "bugsigdb", "refseq"
#' @param ... passed to each function
#' @return list of data frames and ggplot2 object
#' 
#' @examples
#' geneList <- c("DDX41")
#' osplot("refseq", geneList)
#' @export
#' 
osplot <- function(target, ...) {
	if (target=="abstract"){
		return(wcAbst(...))
	} else if (target=="bugsigdb"){
		return(wcBSDB(...))
	} else if (target=="refseq"){
		return(wcGeneSummary(...))
	} else {
		stop("Please specify abstract, bugsigdbr, or refseq")
	}
}