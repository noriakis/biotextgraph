#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' 
#' @param mbList microbe list
#' @param excludeFreq exclude words with overall frequency above excludeFreq
#'                    default to 0
#' @param additionalRemove specific words to be excluded
#' @param mbPlot plot microbe names
#' @param disPlot plot diseases
#' @param madeUpper make the words uppercase in resulting plot
#' @param pal palette for color gradient in correlation network and tag coloring
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
#' @param ngram default to NA (1)
#' @param tag perform pvclust on words and colorlize them in wordcloud
#' @param excludeTfIdf exclude based on tfidf (default: NA)
#' @param ... parameters to pass to wordcloud()
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import bugsigdbr
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @import pvclust
#' @import ggforce
#' @importFrom dplyr filter
#' @importFrom stats dist
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency graph_from_edgelist
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' 
#' @examples
#' mbList <- c("Veillonella dispar")
#' wcBSDB(mbList)
#' @export
#' 
wcBSDB <- function (mbList,
                    excludeFreq=1000, excludeTfIdf=NA,
                    additionalRemove=NA,
                    madeUpper=c("dna","rna"),
                    pal=c("blue","red"), numWords=15,
                    scaleRange=c(5,10), showLegend=FALSE,
                    edgeLabel=FALSE, mbPlot=FALSE,
                    ngram=NA, plotType="wc", disPlot=NA,
                    colorText=FALSE, corThresh=0.6, tag=FALSE,
                    layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, ...) {

    qqcat("Input microbes: @{length(mbList)}\n")
    returnList <- list()

    ## Load from GeneSummary
    tb <- importBugSigDB()

    ## Filter high frequency words
    filterWords <- allFreqBSDB[
        allFreqBSDB$freq > excludeFreq,]$word
    if (!is.na(excludeTfIdf)){
        filterWords <- c(filterWords,
            allTfIdfBSDB[
                allTfIdfBSDB$tfidf > excludeTfIdf,]$word)
    }
    qqcat("filtered @{length(filterWords)} words (frequency) ...\n")
    # filterWords <- c(filterWords, "pmids", "geneid") 
    ## Excluded by default

    # if (ora){
    #     qqcat("performing ORA\n")
    #     sig <- textORA(mbList)
    #     sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
    #     qqcat("filtered @{length(sigFilter)} words (ORA) ...\n")
    #     filterWords <- c(filterWords, sigFilter)
    #     returnList[["ora"]] <- sig
    # }
    
    subTb <- c()        
    for (m in mbList) {
        tmp <- tb[grepl(m, tb$`MetaPhlAn taxon names`,
                 fixed = TRUE),]
        tmp$query <- m
        subTb <- rbind(subTb, tmp)
    }
    returnList[["subsetBSDB"]] <- subTb
    fil <- subTb[!duplicated(subTb$Title),]


    ## Make corpus
    docs <- VCorpus(VectorSource(fil$Title))
    docs <- makeCorpus(docs, filterWords, additionalRemove)

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

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }

    returnList[["rawfrequency"]] <- matSorted
    returnList[["TDM"]] <- docs

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        freqWords <- names(matSorted)
        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        
        if (tag) {
            pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
            pvcl <- pvpick(pvc)
            returnList[["pvcl"]] <- pvcl
        }

        ## genePlot: plot associated genes
        if (!is.na(disPlot)) {mbPlot <- TRUE}
        if (mbPlot) {
            row.names(freqWordsDTM) <- fil$query
        }
        
        if (disPlot) {
            dis <- c()
            for (i in seq_len(nrow(subTb))){
                dis <- rbind(dis,
                c(subTb[i, "Condition"], subTb[i, "query"]))
            }
            dis <- dis[!is.na(dis[,1]),]
            dis <- dis[!is.na(dis[,2]),]
            dmg <- simplify(graph_from_data_frame(dis, directed=FALSE))
        }

        ## Check correlation
        corData <- cor(freqWordsDTM)
        returnList[["corMat"]] <- corData

        ## Set correlation below threshold to zero
        corData[corData<corThresh] <- 0
        coGraph <- graph.adjacency(corData, weighted=TRUE,
                    mode="undirected", diag = FALSE)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }

        nodeName <- V(coGraph)$name
        dtmCol <- colnames(freqWordsDTM)
        for (i in madeUpper) {
            dtmCol[dtmCol == i] <- toupper(i)
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        colnames(freqWordsDTM) <- dtmCol

        if (mbPlot) {
            mbmap <- c()
            for (rn in nodeName){
                tmp <- freqWordsDTM[ ,rn]
                for (nm in names(tmp[tmp!=0])){
                    mbmap <- rbind(mbmap, c(rn, nm))
                }
            }
            mbmap <- simplify(graph_from_edgelist(mbmap, directed = FALSE))
            coGraph <- igraph::union(coGraph, mbmap)
            
            if (disPlot) {
                coGraph <- igraph::union(coGraph, dmg)
            }
            
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
                corThreshMbPlot <- corThresh - 0.1}
            tmpW[is.na(tmpW)] <- corThreshMbPlot
            E(coGraph)$weight <- tmpW
        }
        


        if (tag) {
            netCol <- tolower(names(V(coGraph)))
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    netCol[netCol==j] <- paste0("cluster",i)
            }
            netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
            V(coGraph)$tag <- netCol
        }

        returnList[["ig"]] <- coGraph

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

        if (tag) {
            netPlot <- netPlot + geom_node_point(aes(size=Freq, color=tag),
                                                show.legend = showLegend)
        } else { 
            netPlot <- netPlot + geom_node_point(aes(size=Freq, color=Freq),
                                                show.legend = showLegend)+
                                 scale_color_gradient(low=pal[1],high=pal[2],
                                                      name = "Frequency")
        }

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
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation")+
            theme_graph()

        returnList[["net"]] <- netPlot
    } else {
        ## WC
        matSorted <- matSorted[1:numWords]
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)

        if (tag) {

            freqWords <- names(matSorted)
            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
            pvcl <- pvpick(pvc)
            returnList[["pvcl"]] <- pvcl

            wcCol <- returnDf$word
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- pal[i]
            }
            wcCol[!wcCol %in% pal] <- "grey"

        }
        for (i in madeUpper) {
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        if (tag){
            wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, 
                                               freq = returnDf$freq,
                                               colors = wcCol,
                                               random.order=FALSE,
                                               ordered.colors = TRUE)))
        } else {
            wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, 
                                               freq = returnDf$freq, ...)))
        }
        returnList[["df"]] <- returnDf
        returnList[["wc"]] <- wc
    }

    return(returnList)
}
