#' textORA
#' 
#' Return p-values from hypergeometric distribution
#' 
#' @param queries gene list (Entrez ID)
#' @return p-values for the words
#' @import GeneSummary tm
#' @return p-values for the words
#' @examples textORA(c("2067","2068","2071","2072"))
#' @export
#' 
textORA <- function(queries) {
    tb <- loadGeneSummary()
    
    fil <- tb |> dplyr::filter(Gene_ID %in% queries)
    
    ## Make corpus for queried genes
    docs <- VCorpus(VectorSource(fil$Gene_summary))
    docs <- docs |>
        tm_map(FUN=content_transformer(tolower)) |> 
        tm_map(FUN=removeNumbers) |>
        tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) |>
        tm_map(FUN=removePunctuation) |>
        tm_map(FUN=stripWhitespace)
    
    docs <- TermDocumentMatrix(docs)
    mat <- as.matrix(docs)
    matSorted <- sort(rowSums(mat), decreasing=TRUE)
    
    pvs <- vapply(names(matSorted),
        function (x) returnP(x, matSorted, allFreqGeneSummary),
        FUN.VALUE=1)
    pvs
}

#' returnP
#' 
#' Calculate p-values
#' @param name name of a gene
#' @param matSorted words frequency in the cluster
#' @param allFreqGeneSummary words frequency in the whole description
#' @return p-values for the words
#' 
returnP <- function(name, matSorted, allFreqGeneSummary){
    query <- as.numeric(matSorted[name])
    noquery <- sum(matSorted) - query
    queryAll <- as.numeric(allFreqGeneSummary[name, "freq"])
    allwords <- sum(allFreqGeneSummary$freq) - queryAll
    
    ## p-value
    return(sum(dhyper(query:sum(matSorted), queryAll,
        allwords, sum(matSorted))))
}