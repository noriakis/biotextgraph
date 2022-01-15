#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' 
#' @param geneList ENTREZID list
#' @param excludeFreq exclude words with overall frequency above excludeFreq, default to 5000
#' @param additionalRemove specific words to be excluded
#' @param madeUpper make the words uppercase in resulting plot
#' @param palette palette for color gradient in correlation network
#' @param numWords the number of words to be shown in correlation network
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param corThresh the correlation threshold
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param ORA perform over-representation analysis
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import GeneSummary
#' @import wordcloud
#' @import igraph
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency
#' @import ggraph ggplot2
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
#' 
#' @examples wcGeneSummary(geneList)
#' @export
#' 
wcGeneSummary <- function (geneList, excludeFreq=5000, additionalRemove=NA, madeUpper=c("dna","rna"), organism=9606,
                           palette=c("blue","red"), numWords=15, scaleRange=c(5,10), ORA=FALSE, showLegend=FALSE,
                           plotType="wc", corThresh=0.6, layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, ...) {
    returnList <- list()
    ## Load from GeneSummary
    tb <- loadGeneSummary(organism = organism)

    # load("allFreqGeneSummary.rda") # Already performed, and automatically loaded

    ## Filter high frequency words
    filterWords <- allFreqGeneSummary[allFreqGeneSummary$freq>excludeFreq,]$word
    filterWords <- c(filterWords, "pmids", "geneid") # 'PMIDs' is excluded by default
    fil <- tb %>% filter(Gene_ID %in% geneList)
    
    ## Make corpus
    docs <- VCorpus(VectorSource(fil$Gene_summary))
    docs <- docs %>%
        tm_map(FUN=content_transformer(tolower)) %>% 
        tm_map(FUN=removeNumbers) %>%
        tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
        tm_map(removeWords, filterWords) %>% 
        tm_map(FUN=removePunctuation) %>%
        tm_map(FUN=stripWhitespace)
    if (prod(is.na(additionalRemove))!=1){
        docs <- docs %>% tm_map(removeWords, additionalRemove)
    }
    
    if (ORA){
        message("This would take time ...")
        allDocs <- VCorpus(VectorSource(tb$Gene_summary))
        allDocs <- allDocs %>%
            tm_map(FUN=content_transformer(tolower)) %>% 
            tm_map(FUN=removeNumbers) %>%
            tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
            tm_map(removeWords, filterWords) %>% 
            tm_map(FUN=removePunctuation) %>%
            tm_map(FUN=stripWhitespace)
        if (prod(is.na(additionalRemove))!=1){
            allDocs <- allDocs %>% tm_map(removeWords, additionalRemove)
        }
        matAll <- as.matrix(TermDocumentMatrix(allDocs))
        matAllSorted <- sort(rowSums(matAll), decreasing=TRUE)
        
        returnp <- function(name){
            query <- as.numeric(matSorted[name])
            noquery <- sum(matSorted) - query
            queryAll <- as.numeric(matAllSorted[name])
            allwords <- sum(matAllSorted) - queryAll
            return(sum(dhyper(query:sum(matSorted), queryAll, allwords, sum(matSorted))))
        }
    }


    ## Set parameters for correlation network
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    docs <- TermDocumentMatrix(docs)
    mat <- as.matrix(docs)
    matSorted <- sort(rowSums(mat), decreasing=TRUE)
    if (ORA){
        ps <- sapply(names(matSorted), returnp)
        returnList[["ps"]] <- ps
    }

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        freqWords <- names(matSorted)
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        
        ## Check correlation
        corData <- cor(freqWordsDTM)

        ## Set correlation below threshold to zero
        corData[corData<corThresh] <- 0
        coGraph <- graph.adjacency(corData, weighted=TRUE, diag = FALSE)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]
        nodeName <- V(coGraph)$name
        for (i in madeUpper) {
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }

        ## Main plot
        netPlot <- ggraph(coGraph, layout=layout)
        if (edgeLink){
            netPlot <- netPlot + geom_edge_link(aes(width=weight, color=weight), alpha=0.5, show.legend = showLegend)
        } else {
            netPlot <- netPlot + geom_edge_diagonal(aes(width=weight, color=weight), alpha=0.5, show.legend = showLegend)
        }
        netPlot <- netPlot + scale_edge_width(range=c(1,3), name = "Correlation")

        ## Add others
        netPlot <- netPlot + geom_node_point(aes(size=Freq, color=Freq), show.legend = showLegend)+
            geom_node_text(aes(label=name, size=Freq), check_overlap=TRUE, repel=TRUE,# size = labelSize,
                           color = "black",
                           bg.color = "white", segment.color="black",
                           bg.r = .15)+
            scale_size(range=scaleRange, name="Frequency")+
            scale_color_gradient(low=palette[1],high=palette[2], name = "Frequency")+
            scale_edge_color_gradient(low=palette[1],high=palette[2], name = "Correlation")+
            theme_graph()
        returnList[["net"]] <- netPlot
        return(returnList)
    } else {
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
