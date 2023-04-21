#' wcGeneSummary
#' 
#' Text mining RefSeq description obtained by GeneSummary
#' 
#' @param geneList gene ID list
#' @param exclude "frequency" or "tfidf",
#' @param excludeFreq default to 5000
#' @param excludeType ">" or "<"
#' @param keyType default to SYMBOL
#' @param additionalRemove specific words to be excluded
#' @param madeUpper make the words uppercase in resulting plot
#' @param madeUpperGenes make genes upper case automatically (default to TRUE)
#' @param pre remove words "pmids", "geneids"
#' @param pal palette for color gradient in correlation network
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param cooccurrence default to FALSE, if TRUE, use cooccurrence instead of correlation
#' @param corThresh the correlation (cooccurrence) threshold
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
#' @param filterMax use pre-calculated filter based on max-values when excludeTfIdf is not null
#' @param genePathPlotSig threshold for adjusted p-values (default: 0.05)
#' @param tag perform pvclust on words and colorlize them in wordcloud
#' @param tagWhole whether to perform pvclust on whole matrix or subset
#' @param pvclAlpha alpha for pvpick()
#' @param onlyTDM return only TDF
#' @param onlyCorpus return only corpus
#' @param tfidf use TfIdf when making TDM
#' @param mergeCorpus specify multiple corpus if intend to combine them.
#'                    like PubMed information and RefSeq summary
#' @param numOnly delete number only (not deleting XXX123)
#' @param bn perform bootstrap-based Bayesian network inference instead of correlation using bnlearn
#' @param R how many bootstrap when bn is stated
#' @param onWholeDTM calculate correlation network
#'                   on whole dataset or top-words specified by numWords
#' @param useUdpipe use udpipe to make a network
#' @param udpipeModel udpipe model file name
#' @param cl for parPvclust, parallel clustering can be performed
#' @param stem whether to use stemming
#' @param tagPalette tag palette when tag is TRUE
#' @param preserve preserve original characters
#' @param takeMax take max values for each term in term-document matrix
#' @param filterMax use pre-calculated filter based on max-values when excluding TfIdf
#' Otherwise take sum.
#' @param collapse default to FALSE, collapse all the sentences
#' @param argList parameters to pass to wordcloud()
#' @param normalize sum normalize the term frequency document-wise
#' @param takeMean take mean values for each term in term-document matrix
#' @param naEdgeColor edge colors for NA values (linking query with the category other than text)
#' @param fontFamily font family to use, default to "sans"
#' @param useggwordcloud default to TRUE
#' @param wcScale scaling size for ggwordcloud
#' @param addFreqToGene add pseudo frequency to gene in genePlot
#' @param colorize color the word nodes by their frequency, and the other nodes by their category
#' if colorize=FALSE and addFreqToGene=TRUE, gene nodes are colorized according to the minimum frequency 
#' of the words in the network
#' @param geneColor color for associated genes with words
#' @param scaleFreq default to NULL, scale the value if specified
#' @param useSeed seed
#' @param scaleEdgeWidth scale for edge width
#' @return list of data frame and ggplot2 object
#' @import tm
#' @import GeneSummary
#' @import org.Hs.eg.db
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @importFrom pvclust pvclust pvpick
#' @import methods
#' @importFrom dplyr filter
#' @importFrom stats dist
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' 
#' @export
#' @examples
#' geneList <- c("DDX41","PNKP")
#' wcGeneSummary(geneList)
#' 
wcGeneSummary <- function (geneList, keyType="SYMBOL",
                            excludeFreq=2000, exclude="frequency",
                            filterMax=FALSE, excludeType=">",
                            tfidf=FALSE, genePlotNum=10,
                            preserve=TRUE, takeMax=FALSE,
                            additionalRemove=NA, onlyCorpus=FALSE,
                            madeUpper=c("dna","rna"), organism=9606,
                            pal=c("blue","red"), numWords=15,
                            scaleRange=c(5,10), showLegend=FALSE,
                            orgDb=org.Hs.eg.db, edgeLabel=FALSE,
                            naEdgeColor="grey50", cooccurrence=FALSE,
                            pvclAlpha=0.95, bn=FALSE, R=20, cl=FALSE,
                            ngram=NA, plotType="network", onlyTDM=FALSE, stem=FALSE,
                            colorText=FALSE, corThresh=0.2, genePlot=FALSE,
                            genePathPlot=NA, genePathPlotSig=0.05, tag=FALSE,
                            layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, 
                            enrich=NULL, topPath=10, ora=FALSE, tagWhole=FALSE,
                            mergeCorpus=NULL, numOnly=TRUE, madeUpperGenes=TRUE,
                            onWholeDTM=FALSE, pre=TRUE, takeMean=FALSE,
                            tagPalette=palette(), collapse=FALSE, addFreqToGene=FALSE,
                            useUdpipe=FALSE, normalize=FALSE, fontFamily="sans",
                            udpipeModel="english-ewt-ud-2.5-191206.udpipe",
                            scaleFreq=NULL, colorize=FALSE, geneColor="grey",
                            argList=list(), useggwordcloud=TRUE, wcScale=10,
                            useSeed=42, scaleEdgeWidth=c(1,3)) {
    ## Make class object to store results
    ret <- new("biotext")
    ret@query <- geneList
    ret@type <- "refseq"

    if (useUdpipe) {
        qqcat("Using udpipe mode\n")
        plotType="network"
        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)}

    if (madeUpperGenes){
        madeUpper <- c(madeUpper, tolower(keys(orgDb, "SYMBOL")))
    }

    if (is.null(mergeCorpus)) {
        qqcat("Input genes: @{length(geneList)}\n")
        if (keyType!="ENTREZID"){
            geneList <- AnnotationDbi::select(orgDb,
                keys = geneList, columns = c("ENTREZID"),
                keytype = keyType)$ENTREZID
            geneList <- geneList[!is.na(geneList)] |> unique()
            qqcat("  Converted input genes: @{length(geneList)}\n")
        }

        ## Filter high frequency words if needed
        if (exclude=="frequency") {
            pref = "GS_Freq"
        } else {
            pref = "GS_TfIdf"
        }
        if (filterMax) {
            pref <- paste0(pref, "_Max")
        }
        filterWords <- retFiltWords(useFil=pref, filType=excludeType, filNum=excludeFreq)

        if (pre){
            filterWords <- c(filterWords, "pmids", "geneid", "pmid", "geneids") ## Excluded by default
        }
        if (length(filterWords)!=0 | length(additionalRemove)!=0){
            allfils <- c(filterWords, additionalRemove)
            allfils <- allfils[!is.na(allfils)]
            if (length(allfils)!=0) {
                ret@filtered <- allfils
            }
        }
        qqcat("Filtered @{length(filterWords)} words (frequency and/or tfidf)\n")

        if (!is.null(enrich)) { ## mining pathway enrichment analysis text
            if (genePlot) {stop("genePlot can't be performed in enrichment analysis mode")}
            qqcat("Performing enrichment analysis\n")
            if (enrich=="reactome"){
                pathRes <- ReactomePA::enrichPathway(geneList)
                pathRes@result$Description <- gsub("Homo sapiens\r: ",
                                "",
                                pathRes@result$Description)                
            } else if (enrich=="kegg"){
                pathRes <- clusterProfiler::enrichKEGG(geneList)
            } else {
                ## Currently only supports some pathways
                stop("Please specify 'reactome' or 'kegg'")
            }
            ret@enrichResults <- pathRes@result
            ## Make corpus
            if (collapse) {
                docs <- VCorpus(VectorSource(paste(pathRes@result$Description[1:topPath], collapse=" ")))
            } else {
                docs <- VCorpus(VectorSource(pathRes@result$Description[1:topPath]))
            }
            if (preserve) {
                pdic <- preserveDict(docs, ngram, numOnly, stem)
                ret@dic <- pdic
            }
            docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly, stem)
        } else { ## Mining RefSeq text 
            tb <- loadGeneSummary(organism = organism)
            fil <- tb %>% filter(tb$Gene_ID %in% geneList)
            fil <- fil[!duplicated(fil$Gene_ID),]
            ret@rawText <- fil

            if (ora){
                qqcat("Performing ORA\n")
                sig <- textORA(geneList)
                sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
                qqcat("Filtered @{length(sigFilter)} words (ORA)\n")
                filterWords <- c(filterWords, sigFilter)
                ret@ora <- sig
            }

            ## Make corpus
            if (collapse) {
                docs <- VCorpus(VectorSource(paste(fil$Gene_summary, collapse=" ")))
            } else {
                docs <- VCorpus(VectorSource(fil$Gene_summary))
            }
            if (preserve) {pdic <- preserveDict(docs, ngram, numOnly, stem)
                        ret@dic <- pdic}
            docs <- makeCorpus(docs, filterWords, additionalRemove, numOnly, stem)
        }

        if (onlyCorpus){
            return(docs)
        }
    } else {
        if (length(mergeCorpus)<2){
            stop("Please provide multile corpus")
        }
        docs <- mergeCorpus
    }

    if (useUdpipe) {
        fil$text <- fil$Gene_summary
        fil$ID <- fil$Gene_ID
        ret <- retUdpipeNet(ret=ret, texts=fil,udmodel_english=udmodel_english,
            orgDb=orgDb, filterWords=filterWords, additionalRemove=additionalRemove,
            colorText=colorText,edgeLink=edgeLink,queryPlot=genePlot, layout=layout,
            pal=pal, showNeighbors=NULL, showFreq=NULL, nodePal=tagPalette)
        return(ret)
    }
    ret@corpus <- docs

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
    if (normalize) {
        mat <- sweep(mat, 2, colSums(mat), `/`)
    }

    if (takeMax & takeMean) {stop("Should specify either of takeMax or takeMean")}
    if (takeMax) {
        perterm <- apply(mat, 1, max, na.rm=TRUE)
    } else {
        if (takeMean) {
            perterm <- apply(mat,1,mean)
        } else {
            perterm <- rowSums(mat)
        }
    }

    matSorted <- sort(perterm, decreasing=TRUE)
    ret@wholeFreq <- matSorted

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }

    ret@TDM <- docs

    if (onlyTDM) {
        return(docs)
    }

    ## Subset to numWords
    ret@numWords <- numWords
    matSorted <- matSorted[1:numWords]
    freqWords <- names(matSorted)

    if (plotType=="network"){
        
        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
        for (i in madeUpper) {
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        ## TODO: Force all functions to return lower freqDf.
        ret@freqDf <- returnDf

        DTM <- t(as.matrix(docs))
        
        if (tag) {
            ret <- tag_words(ret, cl, pvclAlpha, whole=tagWhole, num_words=ret@numWords)
            pvc <- ret@pvclust
            pvcl <- ret@pvpick
        }

        ## genePlot: plot associated genes
        if (!is.na(genePathPlot)) {genePlot <- TRUE}
        if (genePlot) {
            if (!is.null(mergeCorpus)) {
                stop("Cannot perform genePlot when merging corpus")
            }
            revID <- AnnotationDbi::select(orgDb,
                keys = as.character(fil$Gene_ID), 
                columns = c("SYMBOL"),
                keytype = "ENTREZID")$SYMBOL
            row.names(DTM) <- revID
        }

        ## genePathPlot: plot associated genes and pathways
        if (!is.na(genePathPlot)) {
            
            if (genePathPlot == "reactome") {
                pathRes <- ReactomePA::enrichPathway(geneList)
                pathRes@result$Description <- gsub("Homo sapiens\r: ",
                                "",
                                pathRes@result$Description)
            }
            else if (genePathPlot == "kegg") {
                pathRes <- clusterProfiler::enrichKEGG(geneList)
            }
            else {
                stop("Please specify 'reactome' or 'kegg'")
            }
            
            if (dim(subset(pathRes@result, p.adjust<0.05))[1]==0) {
                stop("No enriched term found.")
            } else {
                qqcat("Found @{dim(subset(pathRes@result, p.adjust<0.05))[1]} enriched term\n")
            }
            if (genePathPlot=="kegg"){pathRes@keytype <- "ENTREZID"}
            ret@enrichResults <- pathRes@result
            sigPath <- subset(clusterProfiler::setReadable(pathRes,
                orgDb)@result, p.adjust<genePathPlotSig)
            pathGraph <- c()
            for (i in 1:nrow(sigPath)){
                pa <- sigPath[i, "Description"]
                for (j in unlist(strsplit(sigPath[i, "geneID"], "/"))){
                    pathGraph <- rbind(pathGraph, c(pa, j))
                }
            }
        }

        matrixs <- obtainMatrix(ret, bn, R, DTM, freqWords,
            corThresh, cooccurrence, onWholeDTM)
        coGraph <- matrixs$coGraph
        ret <- matrixs$ret

        ret@igraphRaw <- coGraph

        coGraph <- induced.subgraph(coGraph,
            names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph,
                degree(coGraph) > 0)
        }

        ## Make uppercase
        nodeName <- V(coGraph)$name
        dtmCol <- colnames(DTM)
        for (i in madeUpper) {
            dtmCol[dtmCol == i] <- toupper(i)
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        colnames(DTM) <- dtmCol

        if (genePlot) {
            genemap <- c()
            for (rn in nodeName){
                tmp <- DTM[, rn]
                for (nm in names(tmp[tmp!=0])){
                    genemap <- rbind(genemap, c(rn, nm))
                }
            }

            gcnt <- table(genemap[,2])
            gcnt <- gcnt[order(gcnt, decreasing=TRUE)]

            if (length(gcnt)!=1) {ret@geneCount <- gcnt}

            
            incGene <- names(gcnt)[1:genePlotNum]
            genemap <- genemap[genemap[,2] %in% incGene,]
            ret@geneMap <- genemap
            genemap <- simplify(igraph::graph_from_edgelist(genemap,
                directed = FALSE))
            coGraph <- igraph::union(coGraph, genemap)

            E(coGraph)$edgeColor <- E(coGraph)$weight
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
                corThreshGenePlot <- corThresh - 0.1}
            tmpW[is.na(tmpW)] <- corThreshGenePlot
            E(coGraph)$weight <- tmpW
        } else {
            E(coGraph)$edgeColor <- E(coGraph)$weight
        }


        if (!is.na(genePathPlot)) {

            withinCoGraph <- intersect(pathGraph[,2], V(coGraph)$name)
            withinCoGraphPathGraph <- pathGraph[
            pathGraph[,2] %in% withinCoGraph,]

            grp <- c()
            for (i in V(coGraph)$name) {
                if (i %in% withinCoGraphPathGraph[,2]){
                    tmpMap <- withinCoGraphPathGraph[
                    withinCoGraphPathGraph[,2] %in% i,]
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
        if (colorize) {addFreqToGene<-TRUE}

        if (addFreqToGene) {
            ## Set pseudo freq as min value of freq
            fre <- V(coGraph)$Freq
            fre[is.na(fre)] <- min(fre, na.rm=TRUE)
            V(coGraph)$Freq <- fre
        }

        ## Assign node category
        if (genePlot) {
            nodeN <- NULL
            genes <- ret@geneMap[,2] |> unique()
            for (nn in V(coGraph)$name) {
                if (nn %in% genes) {
                    nodeN <- c(nodeN, "Genes")
                } else {
                    nodeN <- c(nodeN, "Words")
                }
            }
            V(coGraph)$nodeCat <- nodeN
        } else {
            nodeN <- rep("Words", length(V(coGraph)))
            V(coGraph)$nodeCat <- nodeN
        }
        names(nodeN) <- V(coGraph)$name

        if (tag) {
            ## If tag=TRUE, significant words are assigned `cluster`
            ## The other words and gene nodes are assigned their category.
            netCol <- tolower(names(V(coGraph)))
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    netCol[netCol==j] <- paste0("cluster",i)
            }
            netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
            V(coGraph)$tag <- netCol

            if (!is.null(nodeN)) {
              addC <- V(coGraph)$tag
              for (nn in seq_along(names(V(coGraph)))) {
                if (V(coGraph)$tag[nn]!="not_assigned"){next}
                if (names(V(coGraph))[nn] %in% names(nodeN)) {
                  addC[nn] <- nodeN[names(V(coGraph))[nn]]
                } else {
                  next
                }
              }
              V(coGraph)$tag <- addC
            }
        }

        if (preserve) {
            newGname <- NULL
            for (nm in names(V(coGraph))) {
              if (nm %in% names(pdic)) {
                newGname <- c(newGname, pdic[nm])
              } else {
                newGname <- c(newGname, nm)
              }
            }
            coGraph <- set.vertex.attribute(coGraph, "name", value=newGname)
        }

        ret@igraph <- coGraph

        ## Main plot
        E(coGraph)$weightLabel <- round(E(coGraph)$weight, 3)
        netPlot <- ggraph(coGraph, layout=layout)

        netPlot <- appendEdges(netPlot, bn, edgeLink,
            edgeLabel, showLegend, fontFamily)

        if (tag) { ## Obtain tag coloring
            cols <- V(coGraph)$tag |> unique()
            tagPalette <- tagPalette[1:length(cols)]
            names(tagPalette) <- cols
            tagPalette["Genes"] <- geneColor
        }

        if (colorize) {
            netPlot <- netPlot + 
                geom_node_point(aes(
                    size=.data$Freq,
                    color=.data$Freq),
                    show.legend = showLegend)+
                scale_color_gradient(low=pal[1],high=pal[2],
                                    name = "Frequency")+
                geom_node_point(aes(size=.data$Freq,
                    filter=.data$name %in% incGene),
                show.legend = FALSE, color=geneColor)

        } else {
            if (tag) {
                netPlot <- netPlot + geom_node_point(aes(size=.data$Freq,
                    color=.data$tag), show.legend = showLegend)+
                scale_color_manual(values=tagPalette)
            } else { 
                netPlot <- netPlot + geom_node_point(aes(size=.data$Freq,
                    color=.data$Freq), show.legend = showLegend)+
                scale_color_gradient(low=pal[1],high=pal[2],
                                    name = "Frequency")
            }
        }

        if (colorText){
            if (colorize) {
                    ## This obtain ggrepel text position and 
                    ## show text based on these potisions
                    netPlot <- netPlot + 
                        geom_node_text(aes(label=.data$name, size=.data$Freq,
                            color=.data$Freq),
                            check_overlap=TRUE, repel=TRUE,
                            family=fontFamily,
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=FALSE)

                    layerNum <- length(netPlot$layers)
                    geom_param_list <- netPlot$layers[[layerNum]]$geom_params
                    build <- ggplot_build(netPlot)$data[[layerNum]]
                    build[ netPlot$data$name %in% incGene, ]$colour <- geneColor

                    aes_list <- netPlot$layers[[layerNum]]$mapping
                    aes_list["filter"] <- NULL
                    geom_param_list[["repel"]] <- TRUE

                    geom_param_list[["seed"]] <- useSeed
                    netPlot$layers[[layerNum]]$geom_params[["seed"]] <- useSeed
                    geom_param_list[["show.legend"]] <- FALSE
                    geom_param_list["color"] <- NULL;
                    aes_list["label"] <- NULL
                    aes_list["color"] <- NULL
                    aes_list["colour"] <- NULL
                    geom_param_list["na.rm"] <- NULL
                    aes_params <- netPlot$layers[[layerNum]]$aes_params

                    netPlot <- netPlot + do.call(geom_node_text,
                                    c(aes_params,
                                      list(mapping=aes_list),
                                      geom_param_list, list(
                                       label=build$label,
                                        color=build$colour)))

            } else {
                if (tag) {
                    netPlot <- netPlot + 
                        geom_node_text(aes(label=.data$name,
                            size=.data$Freq, color=.data$tag),
                            family=fontFamily,
                            check_overlap=TRUE, repel=TRUE,
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend)
                } else {
                    netPlot <- netPlot + 
                        geom_node_text(aes(label=.data$name, size=.data$Freq,
                            color=.data$Freq),
                            check_overlap=TRUE, repel=TRUE,
                            family=fontFamily,
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend)
                }
            }
        } else {
            netPlot <- netPlot +
                        geom_node_text(aes(label=.data$name, size=.data$Freq),
                            check_overlap=TRUE, repel=TRUE,# size = labelSize,
                            color = "black",
                            family=fontFamily,
                            bg.color = "white", segment.color="black",
                            bg.r = .15, show.legend=showLegend) 
        }

        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=scaleEdgeWidth, name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation", na.value=naEdgeColor)+
            theme_graph()

        if (!is.na(genePathPlot)) {
            netPlot <- netPlot + ggforce::geom_mark_hull(
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
        ret@net <- netPlot
    } else {
        ## WC

        returnDf <- data.frame(word = names(matSorted),freq=matSorted)

        if (tag) {
            freqWords <- names(matSorted)
            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            if (tagWhole){
                pvc <- pvclust(as.matrix(dist(as.matrix(docs))), parallel=cl)
            } else {
                pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))), parallel=cl)
            }
            pvcl <- pvpick(pvc)
            ret@pvclust <- pvc
            ret@pvpick <- pvcl

            wcCol <- returnDf$word
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- tagPalette[i]
            }
            wcCol[!wcCol %in% tagPalette] <- "grey"

        }
        for (i in madeUpper) {
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        if (preserve) {
            for (nm in unique(returnDf$word)) {
                if (nm %in% names(pdic)) {
                    returnDf[returnDf$word == nm, "word"] <- pdic[nm]
                }
            }
        }
        
        if (!is.null(scaleFreq)) {
            showFreq <- returnDf$freq*scaleFreq
        } else {
            showFreq <- returnDf$freq
        }

        if (tag){
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- showFreq
            argList[["family"]] <- fontFamily
            argList[["colors"]] <- wcCol
            argList[["random.order"]] <- FALSE
            argList[["ordered.colors"]] <- TRUE
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill="transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        } else {
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- showFreq
            argList[["family"]] <- fontFamily
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill = "transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        }
        ret@freqDf <- returnDf
        ret@wc <- wc
    }

    return(ret)
}
