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
#' @importFrom igraph graph.adjacency
#' 
#' @examples wcGeneSummary(geneList)
#' @export
#' 
wcGeneSummary <- function (geneList, excludeFreq=5000, additionalRemove=NA, madeUpper=c("dna","rna"), organism=9606,
                           palette=c("blue","red"), numWords=15, labelSize=5,
                           plotType="wc", corThresh=0.6, layout="nicely", ...) {
    returnList <- list()
    tb <- loadGeneSummary(organism = organism)

    # load("allFreqGeneSummary.rda") ## Already performed
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
    
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    docs <- TermDocumentMatrix(docs)
    
    if (plotType=="network"){
        mft <- findMostFreqTerms(docs, n = numWords)
        freq <- mft[[1]]
        freqWords <- names(mft[[1]])
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        corData <- cor(freqWordsDTM)
        corData[corData<corThresh] <- 0
        coGraph <- graph.adjacency(corData, weighted=TRUE, diag = FALSE)
        V(coGraph)$Freq <- freq[V(coGraph)$name]
        nodeName <- V(coGraph)$name
        for (i in madeUpper) {
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        netPlot <- ggraph(coGraph, layout=layout) +
            geom_edge_link(aes(width=weight, color=weight), alpha=0.5, show.legend = F)+
            geom_node_point(aes(size=Freq, color=Freq), show.legend = F)+
            geom_node_text(aes(label=name), check_overlap=TRUE, repel=TRUE, size = labelSize,
                           color = "black",
                           bg.color = "white", segment.color="black",
                           bg.r = .15)+
            scale_size(range=c(2,10))+
            scale_edge_width(range=c(1,3))+
            scale_color_gradient(low=palette[1],high=palette[2])+
            scale_edge_color_gradient(low=palette[1],high=palette[2])+
            theme_graph()
        return(netPlot)
    } else {
        mat <- as.matrix(docs)
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
}
