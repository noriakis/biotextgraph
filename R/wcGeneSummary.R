#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' 
#' @param geneList ENTREZID list
#' @param excludeFreq exclude words with overall frequency above excludeFreq, default to 5000
#' @param additionalRemove NA
#' @param madeUpper make the words uppercase in resulting plot
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import dplyr
#' @import GeneSummary
#' @import wordcloud
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
#' 
#' @examples wcGeneSummary(geneList)
#' @export
#' 
wcGeneSummary <- function (geneList,
                           excludeFreq=5000,
                           additionalRemove=NA,
                           madeUpper=c("dna","rna"),
                           organism=9606,
                           ...) {
    returnList <- list()
    tb <- loadGeneSummary(organism = organism)

    load("allFreqGeneSummary.rda") ## Already performed
    filterWords <- allFreqGeneSummary[allFreqGeneSummary$freq>excludeFreq,]$word
    fil <- tb %>% filter(Gene_ID %in% geneList)
    
    docs <- VCorpus(VectorSource(fil$Gene_summary))
    docs <- docs %>%
        tm_map(FUN=content_transformer(tolower)) %>% 
        tm_map(FUN=removeNumbers) %>%
        tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
        tm_map(removeWords, filterWords) %>% 
        tm_map(FUN=removePunctuation) %>%
        tm_map(FUN=stripWhitespace)
    if (!is.na(additionalRemove)){
        docs <- docs %>% tm_map(removeWords, additionalRemove)
    }
    mat <- as.matrix(TermDocumentMatrix(docs))
    matSorted <- sort(rowSums(mat),decreasing=TRUE)
    returnDf <- data.frame(word = names(matSorted),freq=matSorted)
    for (i in madeUpper) {
        # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
        returnDf[returnDf$word == i,"word"] = toupper(i)
    }
    

    wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, freq = returnDf$freq, ...)))
    returnList[["df"]] <- returnDf
    returnList[["wc"]] <- wc
    return(returnList)
}
