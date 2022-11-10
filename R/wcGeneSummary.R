#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' 
#' @param geneList gene ID list
#' @param excludeFreq exclude words with overall frequency above excludeFreq
#'                    default to 5000
#' @param keyType default to SYMBOL
#' @param additionalRemove specific words to be excluded
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
#' @param orgDb default to org.Hs.eg.db
#' @param organism organism ID to use
#' @param enrich currently, only 'reactome' and 'kegg' is supported
#' @param topPath how many pathway descriptions are included in text analysis
#' @param ora perform ora or not (experimental)
#' @param ngram default to NA (1)
#' @param genePlot plot associated genes (default: FALSE)
#' @param genePlotNum number of genes to be plotted
#' @param genePathPlot plot associated genes and pathways (default: FALSE)
#'                     "kegg" or "reactome"
#' @param genePathPlotSig threshold for adjusted p-values (default: 0.05)
#' @param tag perform pvclust on words and colorlize them in wordcloud
#' @param pvclAlpha alpha for pvpick()
#' @param excludeTfIdf exclude based on tfidf (default: NA)
#' @param onlyTDM return only TDF
#' @param onlyCorpus return only corpus
#' @param tfidf use TfIdf when making TDM
#' @param mergeCorpus specify multiple corpus if intend to combine them
#' @param numOnly delete number only
#' @param bn perform bootstrap-based Bayesian network inference instead of correlation
#' @param R how many bootstrap when bn is stated
#' @param ... parameters to pass to wordcloud()
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import GeneSummary
#' @import org.Hs.eg.db
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @import pvclust
#' @import ggforce
#' @importFrom dplyr filter
#' @importFrom bnlearn boot.strength averaged.network as.igraph
#' @importFrom stats dist
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' @importFrom ReactomePA enrichPathway
#' @importFrom clusterProfiler enrichKEGG setReadable
#' 
#' @examples
#' geneList <- c("DDX41")
#' wcGeneSummary(geneList)
#' @export
#' 
wcGeneSummary <- function (geneList, keyType="SYMBOL",
                            excludeFreq=2000, excludeTfIdf=NA,
                            tfidf=FALSE, genePlotNum=10,
                            additionalRemove=NA, onlyCorpus=FALSE,
                            madeUpper=c("dna","rna"), organism=9606,
                            pal=c("blue","red"), numWords=15,
                            scaleRange=c(5,10), showLegend=FALSE,
                            orgDb=org.Hs.eg.db, edgeLabel=FALSE,
                            pvclAlpha=0.95, bn=FALSE, R=20,
                            ngram=NA, plotType="wc", onlyTDM=FALSE,
                            colorText=FALSE, corThresh=0.2, genePlot=FALSE,
                            genePathPlot=NA, genePathPlotSig=0.05, tag=FALSE,
                            layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, 
                            enrich=NULL, topPath=10, ora=FALSE,
                            mergeCorpus=NULL, numOnly=TRUE, ...) {
    if (is.null(mergeCorpus)) {
        qqcat("Input genes: @{length(geneList)}\n")
        if (keyType!="ENTREZID"){
            geneList <- AnnotationDbi::select(orgDb,
                keys = geneList, columns = c("ENTREZID"),
                keytype = keyType)$ENTREZID
            geneList <- geneList[!is.na(geneList)]
            qqcat("Converted input genes: @{length(geneList)}\n")
        }

        returnList <- list()

        ## Filter high frequency words if needed
        filterWords <- allFreqGeneSummary[
                                allFreqGeneSummary$freq > excludeFreq,]$word
        if (!is.na(excludeTfIdf)){
            filterWords <- c(filterWords,
                allTfIdfGeneSummary[
                                allTfIdfGeneSummary$tfidf > excludeTfIdf,]$word)
        }
        filterWords <- c(filterWords, "pmids", "geneid") ## Excluded by default

        qqcat("filtered @{length(filterWords)} words (frequency | tfidf) ...\n")

        ## If specified pathway option
        if (!is.null(enrich)) {
            if (genePlot) {stop("genePlot cant be performed in enrichment analysis mode")}
            qqcat("performing enrichment analysis ...")
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
            docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly)
        } else {
            ## Load from GeneSummary
            tb <- loadGeneSummary(organism = organism)

            ## Already performed, and automatically loaded
            # load("allFreqGeneSummary.rda") 

            if (ora){
                qqcat("performing ORA\n")
                sig <- textORA(geneList)
                sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
                qqcat("filtered @{length(sigFilter)} words (ORA) ...\n")
                filterWords <- c(filterWords, sigFilter)
                returnList[["ora"]] <- sig
            }

            fil <- tb %>% filter(Gene_ID %in% geneList)
            fil <- fil[!duplicated(fil$Gene_ID),]

            returnList[["rawtext"]] <- fil
            ## Make corpus
            docs <- VCorpus(VectorSource(fil$Gene_summary))
            docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly)
        }

        if (onlyCorpus){
            return(docs)
        }
    } else {
        if (length(mergeCorpus)<2){
            stop("Please provide multile corpus")
        }
        returnList <- list()
        docs <- mergeCorpus
    }

    ## Set parameters for correlation network
    if (is.na(corThresh)){corThresh<-0.6}
    if (is.na(numWords)){numWords<-10}
    if (!is.na(ngram)){
        NgramTokenizer <- function(x)
            unlist(lapply(ngrams(words(x), ngram),
                paste, collapse = " "),
                use.names = FALSE)
        if (tfidf) {
            docs <- TermDocumentMatrix(docs,
                control = list(tokenize = NgramTokenizer,
                    weighting = weightTfIdf))
        } else {
            docs <- TermDocumentMatrix(docs,
                control = list(tokenize = NgramTokenizer))
        }
    } else {
        if (tfidf) {
            docs <- TermDocumentMatrix(docs,
                control = list(weighting = weightTfIdf))
        } else {
            docs <- TermDocumentMatrix(docs)
        }
    }
    mat <- as.matrix(docs)
    matSorted <- sort(rowSums(mat), decreasing=TRUE)

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }

    returnList[["rawfrequency"]] <- matSorted
    returnList[["TDM"]] <- docs

    if (onlyTDM) {
        return(docs)
    }

    if (plotType=="network"){
        matSorted <- matSorted[1:numWords]
        freqWords <- names(matSorted)
        # freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        ## TODO: before or after?
        freqWordsDTM <- t(as.matrix(docs))
        
        if (tag) {
            ## TODO: tagging based on cluster_walktrap
            pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
            pvcl <- pvpick(pvc, alpha=pvclAlpha)
            returnList[["pvcl"]] <- pvcl
        }

        ## genePlot: plot associated genes
        if (!is.na(genePathPlot)) {genePlot <- TRUE}
        if (genePlot) {
            if (!is.null(mergeCorpus)) {
                stop("cannot perform genePlot when merging corpus ...")
            }
            revID <- AnnotationDbi::select(orgDb,
                keys = as.character(fil$Gene_ID), 
                columns = c("SYMBOL"),
                keytype = "ENTREZID")$SYMBOL
            row.names(freqWordsDTM) <- revID
        }

        ## genePathPlot: plot associated genes and pathways
        if (!is.na(genePathPlot)) {
            
            if (genePathPlot == "reactome") {
                pathRes <- enrichPathway(geneList)
            }
            else if (genePathPlot == "kegg") {
                pathRes <- enrichKEGG(geneList)
            }
            else {
                stop("please specify 'reactome' or 'kegg'")
            }
            
            if (dim(subset(pathRes@result, p.adjust<0.05))[1]==0) {
                stop("no enriched term found.")
            } else {
                qqcat("found @{dim(subset(pathRes@result, p.adjust<0.05))[1]} enriched term ...\n")
            }
            if (genePathPlot=="kegg"){pathRes@keytype <- "ENTREZID"}
            returnList[["pathRes"]] <- pathRes
            sigPath <- subset(setReadable(pathRes, orgDb)@result, p.adjust<genePathPlotSig)
            pathGraph <- c()
            for (i in 1:nrow(sigPath)){
                pa <- sigPath[i, "Description"]
                for (j in unlist(strsplit(sigPath[i, "geneID"], "/"))){
                    pathGraph <- rbind(pathGraph, c(pa, j))
                }
            }
        }

        if (bn) {
            qqcat("bn specified, R=@{R} ...\n")
            # To avoid computaitonal time, subset to numWords
            bnboot <- boot.strength(
                data.frame(freqWordsDTM[,colnames(freqWordsDTM) %in% freqWords]),
                algorithm = "hc", R=R)
            returnList[["strength"]] <- bnboot
            av <- averaged.network(bnboot)
            avig <- as.igraph(av)
            el <- data.frame(as_edgelist(avig))
            colnames(el) <- c("from","to")
            mgd <- merge(el, bnboot, by=c("from","to"))
            colnames(mgd) <- c("from","to","weight","direction")
            coGraph <- graph_from_data_frame(mgd, directed=TRUE)
        } else {
            ## Check correlation
            corData <- cor(freqWordsDTM)
            returnList[["corMat"]] <- corData

            ## Set correlation below threshold to zero
            corData[corData<corThresh] <- 0
            coGraph <- graph.adjacency(corData, weighted=TRUE,
                        mode="undirected", diag = FALSE)
        }
        ## before or after?
        coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
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

        if (genePlot) {
            genemap <- c()
            for (rn in nodeName){
                tmp <- freqWordsDTM[ ,rn]
                for (nm in names(tmp[tmp!=0])){
                    genemap <- rbind(genemap, c(rn, nm))
                }
            }

            gcnt <- table(genemap[,2])
            gcnt <- gcnt[order(gcnt, decreasing=TRUE)]
            returnList[["geneCount"]] <- gcnt
            
            incGene <- names(gcnt)[1:genePlotNum]
            genemap <- genemap[genemap[,2] %in% incGene,]

            genemap <- simplify(igraph::graph_from_edgelist(genemap, directed = FALSE))
            coGraph <- igraph::union(coGraph, genemap)
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
                corThreshGenePlot <- corThresh - 0.1}
            tmpW[is.na(tmpW)] <- corThreshGenePlot
            E(coGraph)$weight <- tmpW
        }


        if (!is.na(genePathPlot)) {

            withinCoGraph <- intersect(pathGraph[,2], V(coGraph)$name)
            withinCoGraphPathGraph <- pathGraph[ pathGraph[,2] %in% withinCoGraph,]

            grp <- c()
            for (i in V(coGraph)$name) {
                if (i %in% withinCoGraphPathGraph[,2]){
                    tmpMap <- withinCoGraphPathGraph[withinCoGraphPathGraph[,2] %in% i,]
                    if (is.vector(tmpMap)) {
                        grp <- c(grp, tmpMap[1])
                    } else {
                        tmpMap <- tmpMap[order(tmpMap[,1]),]
                        grp <- c(grp, paste(tmpMap[,1], collapse="\n"))
                    }
                } else {
                    grp <- c(grp, NA)
                }
            }
            V(coGraph)$grp <- grp
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

        if (bn){
            if (edgeLink){
                if (edgeLabel){
                    netPlot <- netPlot +
                                geom_edge_link(
                                    aes(width=weight,
                                    color=weight,
                                    label=round(weight,3)),
                                    angle_calc = 'along',
                                    label_dodge = unit(2.5, 'mm'),
                                    arrow = arrow(length = unit(4, 'mm')), 
                                    start_cap = circle(3, 'mm'),
                                    end_cap = circle(3, 'mm'),
                                    alpha=0.5,
                                    show.legend = showLegend)
                } else {
                    netPlot <- netPlot +
                                geom_edge_link(aes(width=weight, color=weight),
                                    arrow = arrow(length = unit(4, 'mm')), 
                                    start_cap = circle(3, 'mm'),
                                    end_cap = circle(3, 'mm'),
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
                                    arrow = arrow(length = unit(4, 'mm')), 
                                    start_cap = circle(3, 'mm'),
                                    end_cap = circle(3, 'mm'),
                                    alpha=0.5,
                                    show.legend = showLegend)
                } else {
                    netPlot <- netPlot +
                                geom_edge_diagonal(aes(width=weight, color=weight),
                                    arrow = arrow(length = unit(4, 'mm')), 
                                    start_cap = circle(3, 'mm'),
                                    end_cap = circle(3, 'mm'),                                    
                                    alpha=0.5, show.legend = showLegend)                
                }
            }
        } else {
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
            if (tag) {
                netPlot <- netPlot + 
                    geom_node_text(aes(label=name, size=Freq, color=tag),
                        check_overlap=TRUE, repel=TRUE,# size = labelSize,
                        bg.color = "white", segment.color="black",
                        bg.r = .15, show.legend=showLegend)
            } else {
                netPlot <- netPlot + 
                    geom_node_text(aes(label=name, size=Freq, color=Freq),
                        check_overlap=TRUE, repel=TRUE,# size = labelSize,
                        bg.color = "white", segment.color="black",
                        bg.r = .15, show.legend=showLegend)
            }
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

        if (!is.na(genePathPlot)) {
            netPlot <- netPlot + geom_mark_hull(
                aes(netPlot$data$x,
                    netPlot$data$y,
                    group = grp,
                    label=grp, fill=grp,
                    filter = !is.na(grp)),
                concavity = 4,
                expand = unit(2, "mm"),
                alpha = 0.25,
                na.rm = TRUE,
                # label.fill="transparent",
                show.legend = FALSE
            )
        }
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

        if (tfidf) {
            showFreq <- returnDf$freq*10
        } else {
            showFreq <- returnDf$freq
        }

        if (tag){
            wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, 
                                               freq = showFreq,
                                               colors = wcCol,
                                               random.order=FALSE,
                                               ordered.colors = TRUE)))
        } else {
            wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, 
                                               freq = showFreq, ...)))
        }
        returnList[["df"]] <- returnDf
        returnList[["wc"]] <- wc
    }

    if (!is.null(enrich)){
        returnList[["enrich"]] <- pathRes
    }

    return(returnList)
}
