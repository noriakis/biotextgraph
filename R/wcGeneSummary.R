#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' 
#' @param geneList ENTREZID list
#' @param excludeFreq exclude words with overall frequency above excludeFreq
#'                    default to 5000
#' @param keyType default to SYMBOL
#' @param additionalRemove specific words to be excluded
#' @param madeUpper make the words uppercase in resulting plot
#' @param palette palette for color gradient in correlation network
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param corThresh the correlation threshold
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @param colorText color text label based on frequency in correlation network
#' @param orgDb default to org.Hs.eg.db
#' @param organism organism ID to use
#' @param enrich currently, only 'reactome' and 'kegg' is supported
#' @param topPath how many pathway descriptions are included in text analysis
#' @param ora perform ora or not (experimental)
#' @param ngram default to NA (1)
#' @param ... parameters to pass to wordcloud()
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import GeneSummary
#' @import org.Hs.eg.db
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @importFrom dplyr filter
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichKEGG
#' 
#' @examples
#' geneList <- c("DDX41")
#' wcGeneSummary(geneList)
#' @export
#' 
wcGeneSummary <- function (geneList, keyType="SYMBOL",
                            excludeFreq=5000, additionalRemove=NA,
                            madeUpper=c("dna","rna"), organism=9606,
                            palette=c("blue","red"), numWords=15,
                            scaleRange=c(5,10), showLegend=FALSE,
                            orgDb=org.Hs.eg.db, edgeLabel=FALSE,
                            ngram=NA, plotType="wc",
                            colorText=FALSE, corThresh=0.6,
                            layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, 
                            enrich=NULL, topPath=10, ora=FALSE, ...) {

    qqcat("Input genes: @{length(geneList)}\n")
    if (keyType!="ENSEMBL"){
        geneList <- AnnotationDbi::select(orgDb,
            keys = geneList, columns = c("ENTREZID"),
            keytype = keyType)$ENTREZID
        geneList <- geneList[!is.na(geneList)]
        qqcat("Converted input genes: @{length(geneList)}\n")
    }

    returnList <- list()

    ## If specified pathway option
    if (!is.null(enrich)) {
        message("performing enrichment analysis ...")
        if (enrich=="reactome"){
            pathRes <- enrichPathway(geneList)
        } else if (enrich=="kegg"){
            pathRes <- enrichKEGG(geneList)
        } else {
            ## Currently only supports some pathways
            stop("please specify 'reactome' or 'kegg'")
        }
        ## Make corpus
        docs <- VCorpus(VectorSource(pathRes@result$Description[1:topPath]))
        docs <- makeCorpus(docs, additionalRemove, additionalRemove)
    } else {
        ## Load from GeneSummary
        tb <- loadGeneSummary(organism = organism)

        ## Already performed, and automatically loaded
        # load("allFreqGeneSummary.rda") 

        ## Filter high frequency words
        filterWords <- allFreqGeneSummary[
                                allFreqGeneSummary$freq > excludeFreq,]$word
        qqcat("filtered @{length(filterWords)} words ...\n")
        filterWords <- c(filterWords, "pmids", "geneid") 
        ## Excluded by default

        if (ora){
            message("performing ORA\n")
            sig <- textORA(geneList)
            sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
            qqcat("filtered @{length(sigFilter)} words ...\n")
            filterWords <- c(filterWords, sigFilter)
        }

        fil <- tb %>% filter(Gene_ID %in% geneList)
        
        ## Make corpus
        docs <- VCorpus(VectorSource(fil$Gene_summary))
        docs <- makeCorpus(docs, filterWords, additionalRemove)
    }

    ## Set parameters for correlation network
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    if (!is.na(ngram)){
        NgramTokenizer <- function(x)
            unlist(lapply(ngrams(words(x), ngram),
                paste, collapse = " "),
                use.names = FALSE)
        docs <- TermDocumentMatrix(docs,
            control = list(tokenize = NgramTokenizer))
    } else {
        docs <- TermDocumentMatrix(docs)
    }
    mat <- as.matrix(docs)
    matSorted <- sort(rowSums(mat), decreasing=TRUE)
    returnList[["rawfrequency"]] <- matSorted
    returnList[["TDM"]] <- docs

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        freqWords <- names(matSorted)
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        
        ## Check correlation
        corData <- cor(freqWordsDTM)
        returnList[["corMat"]] <- corData

        ## Set correlation below threshold to zero
        corData[corData<corThresh] <- 0
        coGraph <- graph.adjacency(corData, weighted=TRUE,
                    mode="undirected", diag = FALSE)
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
            if (edgeLabel){
                netPlot <- netPlot +
                            geom_edge_link(
                                aes(width=weight,
                                color=weight,
                                label=round(weight,3)),
                                angle_calc = 'along',
                                label_dodge = unit(2.5, 'mm'),
                                alpha=0.5,
                                show.legend = showLegend)
            } else {
                netPlot <- netPlot +
                            geom_edge_link(aes(width=weight, color=weight),
                                alpha=0.5, show.legend = showLegend)
            }
        } else {
            if (edgeLabel){
                netPlot <- netPlot +
                            geom_edge_diagonal(
                                aes(width=weight,
                                color=weight,
                                label=round(weight,3)),
                                angle_calc = 'along',
                                label_dodge = unit(2.5, 'mm'),
                                alpha=0.5,
                                show.legend = showLegend)
            } else {
                netPlot <- netPlot +
                            geom_edge_diagonal(aes(width=weight, color=weight),
                                alpha=0.5, show.legend = showLegend)                
            }
        }
        netPlot <- netPlot + geom_node_point(aes(size=Freq, color=Freq),
                                            show.legend = showLegend)
        if (colorText){
            netPlot <- netPlot +
                        geom_node_text(
                            aes_(label=~name, size=~Freq, color=~Freq),
                            check_overlap=TRUE, repel=TRUE,# size = labelSize,
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend)
        } else {
            netPlot <- netPlot +
                        geom_node_text(aes_(label=~name, size=~Freq),
                            check_overlap=TRUE, repel=TRUE,# size = labelSize,
                            color = "black",
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend) 
        }
        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=c(1,3), name = "Correlation")+
            scale_color_gradient(low=palette[1],high=palette[2],
                name = "Frequency")+
            scale_edge_color_gradient(low=palette[1],high=palette[2],
                name = "Correlation")+
            theme_graph()
        returnList[["net"]] <- netPlot
    } else {
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
        for (i in madeUpper) {
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        wc <- as.ggplot(as_grob(
                ~wordcloud(words = returnDf$word, freq = returnDf$freq, ...)))
        returnList[["df"]] <- returnDf
        returnList[["wc"]] <- wc
    }

    if (!is.null(enrich)){
        returnList[["enrich"]] <- pathRes
    }

    return(returnList)
}
