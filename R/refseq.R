#' refseq, alliance, pubmed, manual, bugsigdb
#' 
#' Text mining RefSeq description, PubMed, BugSigDB and the other manually curated data
#' 
#' @param geneList gene ID list
#' @param queries query ID list
#' @param mbList microbe ID list
#' @param df manual document data.frame (must have column `text`) or vector of text.
#' If `query` column and other columns are present, regards them as category related to the text on the same row.
#' 
#' @param alliancePath path to The Alliance of Genome Resources gene description file.
#' default to "GENE-DESCRIPTION-TSV_HUMAN.tsv"
#' @param plotType "wc" or "network", default to "network"
#' @param exclude "frequency" or "tfidf",
#' @param excludeFreq default to 5000
#' @param excludeType ">" or "<", combined with `exclude` and `excludeFreq`,
#' e.g. filter the words with the pre-calculated frequency > 5000
#' @param keyType default to SYMBOL
#' @param additionalRemove specific words to be excluded
#' @param madeUpper make these words uppercase in resulting plot,
#' default to c("rna" and "dna")
#' @param madeUpperGenes make genes upper case automatically (default to TRUE)
#' This uses the `SYMBOL` key in `orgDb`. 
#' @param pre remove preset filtering words.
#' @param numWords the number of words to be shown in the plot.
#' When `autoThresh` is TRUE, the number of this value will be shown.
#' @param scaleRange scale for label and node size in the network.
#' @param autoScale scale the label and node size automatically for the large network.
#' @param cooccurrence default to FALSE, if TRUE, use cooccurrence instead of correlation.
#' @param corThresh the correlation (cooccurrence) threshold.
#' @param deleteZeroDeg delete zero degree node from plot in the network
#' @param orgDb the database used to convert identifiers, default to org.Hs.eg.db.
#' @param organism organism ID to use in `GeneSummary`
#' @param enrich currently, only 'reactome' and 'kegg' is supported.
#' @param topPath how many pathway descriptions are included in text analysis,
#' sorted by p-values in the results.
#' @param ora perform over-representation analysis or not (experimental)
#' @param ngram N-gram specification, default to 1.
#' @param genePlot plot associated genes (default: FALSE)
#' This option first calculates the high-frequent words, and subsequently calculates
#' the occurrences of these words in each gene, and prioritize the genes by them.
#' @param genePlotNum number of genes to be plotted (default: 10)
#' @param genePathPlot plot associated genes and pathways (default: NULL)
#'                     Must be "kegg" or "reactome", automatically set genePlot to TRUE.
#' @param genePathPlotSig threshold for adjusted p-values (default: 0.05)
#' @param filterMax Use pre-calculated filter based on max-values when excluding TfIdf
#' Otherwise take sum.
#' @param tag perform pvclust on words and colorlize them in wordcloud or network
#' argument of those accepted in pvclust `method.dist` option, like "correlation".
#' Default to "none", which performs no tagging.
#' @param tagWhole whether to perform pvclust on whole matrix or subset of the matrix.
#' @param pvclAlpha alpha value for pvpick()
#' @param onlyTDM return only TDM (tm).
#' @param onlyCorpus return only corpus (tm).
#' @param tfidf use TfIdf when making TDM, default to FALSE.
#' @param mergeCorpus specify multiple corpus if intend to combine them.
#'                    like PubMed information and RefSeq summary
#' @param numOnly delete number only (not deleting XXX123, but delete only the number)
#' @param onWholeDTM calculate correlation network
#'                   on whole dataset or top-words specified by numWords.
#' @param autoThresh automatically choose thresholding value to show the `numWords`,
#' when deleteZeroDeg (deleting no-connected words) is TRUE, which is default.
#' @param autoNumWords determine the number of words to be shown by ORA, default to FALSE.
#' @param useUdpipe use udpipe to make a dependency network.
#' @param udpipeModel udpipe model file name.
#' @param cl for parPvclust, parallel clustering can be performed
#' @param stem whether to use stemming when making corpus.
#' @param preserve Try to preserve original characters.
#' @param takeMax Take max values for each term in term-document matrix
#' @param collapse default to FALSE, collapse all the sentences.
#' @param normalize sum normalize the term frequency document-wise.
#' @param takeMean take mean values for each term in term-document matrix.
#' @param useSeed random seed
#' 
#' @param useggwordcloud default to TRUE, otherwise use `wordcloud` function.
#' @param wcScale scaling size for ggwordcloud
#' @param argList parameters to pass to wordcloud() or ggwordcloud()
#' 
#' @param layout the layout for the network, defaul to "nicely".
#' It can be one of the layouts implemented in `igraph` and `ggraph`, such as
#' `kk` (Kamada-Kawai), `nicely` (automatic selection of algorithm), `drl` (the force-directed DrL layout).
#' The options are available at: https://igraph.org/r/doc/layout_.html
#' @param edgeLink if FALSE, use geom_edge_diagonal. if TRUE, geom_edge_link. Default to TRUE.
#' @param edgeLabel if TRUE, plot the edge label (default: FALSE)
#' @param pal palette for color gradient in correlation network.
#' should be a vector of length two like c("red","blue").
#' @param showLegend whether to show legend in the network
#' @param colorText color text label based on frequency in the network
#' @param tagPalette tag palette when `tag` is TRUE. It is also used for dependency network
#' using udpipe, and tagging colorization for word cloud.
#' Default to NULL, which indicates automatically set.
#' @param naEdgeColor edge colors for NA values (linking query with the category other than text)
#' @param fontFamily font family to use, default to "sans".
#' @param addFreqToGene add pseudo frequency to gene in genePlot, default to FALSE.
#' @param colorize color the word nodes by their frequency, and the other nodes by their category.
#' if colorize=FALSE and addFreqToGene=TRUE, gene nodes are colorized according to the minimum frequency 
#' of the words in the network
#' @param discreteColorWord colorize words by "Words" category, not frequency.
#' @param catColors colors for words and texts when colorize is TRUE and discreteColorWord is TRUE
#' @param geneColor color for associated genes with words (used when tag or colorize option is TRUE)
#' @param scaleFreq default to NULL, scale the value if specified
#' @param scaleEdgeWidth scale for edge width
#' @param splitByEA automatically split the genes based on significant enrichment analysis results,
#' and returns the list of object for each term. Default to NULL. Must be 'kegg' or 'reactome',
#' in which the function performs over-representation analysis by enrichKEGG or enrichPathway in
#' clusterProfiler and ReactomePA.
#' @param filterByGO filter the results to the words obtained from GO terms,
#' while preserving the number of words to be shown
#' @param docsum if TRUE, convert the term-document matrix to binary.
#' @param absolute calculate absolute correlation value
#' @param corOption passed to `cor` function, like list("method"="kendall")
#' @param useRawQuery if you would like to send the query as is, please set this option to TRUE.
#' @param redo if plot in other parameters, input the previous list
#' @param apiKey api key for eutilities
#' @param perQuery search for the queries one by one recursively, not using `delim`.
#' @param retMax how many items are to be retlieved?
#' @param quote whether to quote the queries
#' @param sortOrder sort order, passed to rentrez function
#' @param onlyGene plot only the gene symbol
#' (orgDb with SYMBOL key can be used)
#' @param distinguish_query if TRUE, distinguish query and returned texts
#' by appending (Q) on query
#' @param dateRange if specified, restrict the range of publication date.
#' Must be the two-length vector, like `c("2013/1/1", "2023/1/1")`
#' 
#' @param mbPlot plot microbe names
#' @param disPlot plot diseases
#' @param target "title" or "abstract"
#' @param metab tibble of metabolite - taxon association
#' @param metThresh threshold of association
#' @param metCol metabolite data frame column name in the order of
#' "candidate taxon", "metabolite", "quantitative values for thresholding"
#' @param curate include articles in bugsigdb
#' @param abstArg passed to PubMed function when using curate=FALSE
#' @param mbColor color for Microbes when tagPalette or catColors 
#' is not specified
#' @param ecPlot plot link between enzyme and microbes
#' this option requires two files to be passed to enzyme() and getUPTax().
#' @param ecFile enzyme database file
#' @param upTaxFile UniProt taxonomy file
#' @param addFreqToMB add pseudo frequency to microbes in mbPlot
#' @param udpipeOnlyFreq when using udpipe, include only high-frequent words
#' @param udpipeOnlyFreqNB when using udpipe, include only the neighbors of
#' high-frequent words
#' @param useFil filter based on "GS_TfIdf" (whole gene summary tf-idf)
#'  or "BSDB_TfIdf" (whole bugsigdb tf-idf)
#' @param filNum specify filter tfidf
#' @param filType "above" or "below"
#' @param useQuanteda use quanteda functions to generate
#' @param quantedaArgs list of arguments to be passed to tokens()
#' @param addFreqToNonWords add pseudo-frequency corresponding to minimum
#' frequency of the words to nodes other than words
#' @param delim delimiter for queries
#' @param onlyDf return only the raw data.frame of searching PubMed
#' @param addFreqToQuery add pseudo-frequency to query node
#' @param asis plot the original network (default to FALSE)
#' @param queryColor color for associated queries with words
#' @param queryPlot plot the query in the graph in relation with the words
#' @param x biotext class object
#' @return `biotext` class object
#' 
#' @name generalf
NULL

#' refseq
#' @rdname generalf
#' @import tm
#' @import GeneSummary
#' @import org.Hs.eg.db
#' @import wordcloud
#' @import igraph
#' @importFrom tidygraph graph_join as_tbl_graph activate is.tbl_graph
#' @import ggraph ggplot2
#' @importFrom pvclust pvclust pvpick
#' @import methods
#' @importFrom dplyr filter
#' @importFrom stats dist na.omit
#' @importFrom grDevices palette
#' @importFrom stats as.dendrogram cor dhyper p.adjust
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom NLP ngrams words
#' @importFrom ggplotify as.ggplot
#' @export
#' @examples
#' geneList <- c("DDX41","PNKP","ERCC1","IRF3","XRCC1")
#' refseq(geneList)
#' 
refseq <- function (geneList, keyType="SYMBOL",
    excludeFreq=2000, exclude="frequency",
    filterMax=FALSE, excludeType=">",
    tfidf=FALSE, genePlotNum=10,
    preserve=FALSE, takeMax=FALSE,
    additionalRemove=NA, onlyCorpus=FALSE,
    madeUpper=c("dna","rna"), organism=9606,
    pal=c("blue","red"), numWords=30,
    scaleRange=c(5,10), autoScale=FALSE,
    showLegend=FALSE,
    orgDb=org.Hs.eg.db, edgeLabel=FALSE,
    naEdgeColor="grey50", cooccurrence=FALSE,
    pvclAlpha=0.95, cl=FALSE,
    ngram=1, plotType="network", onlyTDM=FALSE, stem=FALSE,
    colorText=FALSE, corThresh=0.2, genePlot=FALSE,
    autoThresh=TRUE, autoNumWords=FALSE,
    genePathPlot=NULL, genePathPlotSig=0.05, tag="none",
    layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, 
    enrich=NULL, topPath=10, ora=FALSE, tagWhole=FALSE,
    mergeCorpus=NULL, numOnly=TRUE, madeUpperGenes=TRUE,
    onWholeDTM=FALSE, pre=TRUE, takeMean=FALSE,
    tagPalette=NULL, collapse=FALSE, addFreqToGene=FALSE,
    useUdpipe=FALSE, normalize=FALSE, fontFamily="sans",
    udpipeModel="english-ewt-ud-2.5-191206.udpipe",
    scaleFreq=NULL, colorize=FALSE, geneColor="grey",
    argList=list(), useggwordcloud=TRUE, wcScale=10,
    catColors=NULL, discreteColorWord=FALSE,
    useSeed=42, scaleEdgeWidth=c(1,3), splitByEA=NULL,
    filterByGO=FALSE, docsum=FALSE, absolute=TRUE,
    corOption=list()) {
    
    if (!is.null(splitByEA)) {
        if (length(splitByEA)!=1) {
            stop("Please specify kegg or reactome to splitByEA")}
        if ((splitByEA=="kegg")|(splitByEA=="reactome")) {
            return(split_by_ea(as.list(environment())))     
        } else {
            stop("Please specify kegg or reactome to splitByEA")
        }
    }
    
    if (useUdpipe) {
        qqcat("Using udpipe mode\n")
        plotType="network"
        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)
    }

    if (madeUpperGenes){
        madeUpper <- c(madeUpper, tolower(keys(orgDb, "SYMBOL")))
    }

    if (is.null(mergeCorpus)) {
        if (!is.null(enrich)) { ## mining pathway enrichment analysis text
            if (genePlot) {stop("genePlot can't be performed in enrichment analysis mode")}

            ret <- obtain_enrich(geneList, keyType=keyType, enrich=enrich,
                                 org_db=orgDb, top_path=topPath)
            ret <- ret |> set_filter_words(exclude_by=exclude,
                exclude_type=excludeType, exclude="GS",
                exclude_number=excludeFreq, filterMax=filterMax,
                additional_remove=additionalRemove,
                pre=pre, pre_microbe=FALSE)
            ret <- ret |> make_corpus(collapse=collapse,
                                      num_only=numOnly,
                                      stem=stem, preserve=preserve,
                                      ngram=ngram)
            docs <- ret@corpus

        } else { ## Mining RefSeq text 
            ret <- obtain_refseq(geneList, keyType=keyType, organism=organism, org_db=orgDb)
            ret <- ret |> set_filter_words(exclude_by=exclude,
                exclude_type=excludeType, exclude="GS",
                exclude_number=excludeFreq, filterMax=filterMax,
                additional_remove=additionalRemove,
                pre=pre, pre_microbe=FALSE)
            fil <- ret@rawText
            if (autoNumWords) {
                if (keyType!="ENTREZID"){
                    geneList <- AnnotationDbi::select(orgDb,
                        keys = geneList, columns = c("ENTREZID"),
                        keytype = keyType)$ENTREZID
                    geneList <- geneList[!is.na(geneList)] |> unique()
                }
                sig <- textORA(geneList)
                numWords <- names(sig)[p.adjust(sig, "bonferroni")<0.05] |> length()
                if (numWords>50) {numWords <- 50}
                qqcat("numWords is set to @{numWords}\n")
            }
            if (ora){
                qqcat("Performing ORA\n")
                if (!autoNumWords) {
                    if (keyType!="ENTREZID"){
                        geneList <- AnnotationDbi::select(orgDb,
                            keys = geneList, columns = c("ENTREZID"),
                            keytype = keyType)$ENTREZID
                        geneList <- geneList[!is.na(geneList)] |> unique()
                    }
                    sig <- textORA(geneList)                    
                }
                sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
                qqcat("Filtered @{length(sigFilter)} words (ORA)\n")
                ret@filtered <- c(ret@filtered, sigFilter)
                ret@ora <- sig
            }
            ret <- ret |> make_corpus(collapse=collapse,
                                      num_only=numOnly,
                                      stem=stem, preserve=preserve,
                                      ngram=ngram)
            docs <- ret@corpus
            filterWords <- ret@filtered
            ## Udpipe mode can only be used with default RefSeq
            if (useUdpipe) {
                fil$text <- fil$Gene_summary
                fil$ID <- fil$Gene_ID
                ret <- retUdpipeNet(ret=ret, texts=fil,udmodel_english=udmodel_english,
                    orgDb=orgDb, filterWords=filterWords, additionalRemove=additionalRemove,
                    colorText=colorText, edgeLink=edgeLink, queryPlot=genePlot, layout=layout,
                    pal=pal, showNeighbors=NULL, showFreq=NULL, nodePal=tagPalette)
                ret@type <- "udpipe"
                return(ret)
            }
        }
        if (onlyCorpus){
            return(docs)
        }
    } else {
        if (length(mergeCorpus)<2){
            stop("Please provide multile corpus")
        }
        ret <- new("biotext")
        ret@query <- geneList
        ret@type <- "merged"
        docs <- mergeCorpus
    }
    ret@tag <- tag
    ret@corpus <- docs
    pdic <- ret@dic
    ret <- ret |> make_TDM(tfidf=tfidf,
                           normalize=normalize,
                           takeMean=takeMean,
                           takeMax=takeMax,
                           docsum=docsum)

    matSorted <- ret@wholeFreq

    if (numWords > length(matSorted)){
        numWords <- length(matSorted)
    }

    docs <- ret@TDM

    if (onlyTDM) {
        return(docs)
    }

    ## Set parameters for the network
    if (!is.numeric(corThresh)){corThresh<-0.6}
    if (!is.numeric(numWords)){numWords<-30}

    ## If filter by GO terms
    if (filterByGO) {
        qqcat("`filterByGO` option enabled. Filtering by GO terms ...\n")
        if (ngram==1) {
            filtered_by_GO <- names(matSorted)[tolower(names(matSorted)) %in% goWords]
            matSorted <- matSorted[filtered_by_GO]
        } else if (ngram==2) {
            filtered_by_GO <- names(matSorted)[tolower(names(matSorted)) %in% goWords2gram]
            matSorted <- matSorted[filtered_by_GO]          
        } else {# Do nothing
        }
    }

    ## Subset to numWords
    ret@numWords <- numWords
    matSorted <- matSorted[1:numWords]
    freqWords <- names(matSorted)

    if (plotType=="network"){
        incGene <- NULL
        
        returnDf <- data.frame(word = names(matSorted),freq=matSorted) |>
            na.omit()
        for (i in madeUpper) {
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        ## TODO: Force all functions to return lower freqDf.
        ret@freqDf <- returnDf

        DTM <- t(as.matrix(docs))
        
        if (tag!="none") {
            ret <- tag_words(ret, cl, pvclAlpha, whole=tagWhole, num_words=ret@numWords, method=tag)
            pvc <- ret@pvclust
            pvcl <- ret@pvpick
        }

        ## genePlot: plot associated genes
        if (!is.null(genePathPlot)) {genePlot <- TRUE}
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
        if (!is.null(genePathPlot)) {

            if (keyType!="ENTREZID"){
                geneList <- AnnotationDbi::select(orgDb,
                    keys = geneList, columns = c("ENTREZID"),
                    keytype = keyType)$ENTREZID
                geneList <- geneList[!is.na(geneList)] |> unique()
                # qqcat("  Converted input genes: @{length(geneList)}\n")
            }
            
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
            
            if (dim(subset(pathRes@result, p.adjust<genePathPlotSig))[1]==0) {
                stop("No enriched term found.")
            } else {
                qqcat("Found @{dim(subset(pathRes@result, p.adjust<genePathPlotSig))[1]} enriched term\n")
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

        matrixs <- obtainMatrix(ret, FALSE, NULL, DTM, freqWords,
            corThresh, cooccurrence, onWholeDTM, numWords, autoThresh, absolute)

        
        coGraph <- matrixs$coGraph
        ret <- matrixs$ret
        
        # if (tag=="cor") {
        #     ret <- tag_words(ret, cl,
        #         pvclAlpha, whole=tagWhole,
        #         num_words=ret@numWords,
        #         corMat=TRUE, mat=matrixs$ret@corMat)
        #     pvc <- ret@pvclust
        #     pvcl <- ret@pvpick
        # }


        ret@igraphRaw <- coGraph

        coGraph <- induced.subgraph(coGraph,
            names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph,
                degree(coGraph) > 0)
        }

        nodeName <- V(coGraph)$name

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

            coGraph <- tidygraph::graph_join(as_tbl_graph(coGraph),
                as_tbl_graph(genemap), by="name")
            coGraph <- coGraph |> activate("nodes") |>
                mutate(type=ifelse(is.na(.data$Freq),"Genes","Words"))
            # coGraph <- igraph::union(coGraph, genemap)
            E(coGraph)$edgeColor <- E(coGraph)$weight
            tmpW <- E(coGraph)$weight
            if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
                corThreshGenePlot <- corThresh - 0.1}
            tmpW[is.na(tmpW)] <- corThreshGenePlot
            E(coGraph)$weight <- tmpW
        } else {
            coGraph <- as_tbl_graph(coGraph) |> activate("nodes") |>
                mutate(type=ifelse(is.na(.data$Freq),"Genes","Words"))
            E(coGraph)$edgeColor <- E(coGraph)$weight
        }


        if (!is.null(genePathPlot)) {

            withinCoGraph <- intersect(pathGraph[,2], V(coGraph)$name)

            withinCoGraphPathGraph <- pathGraph[
            pathGraph[,2] %in% withinCoGraph,]

            if (is.vector(withinCoGraphPathGraph)) {
                withinCoGraphPathGraph <- data.frame(withinCoGraphPathGraph[1],
                    withinCoGraphPathGraph[2])
            }
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
        nodeN <- (coGraph |> activate("nodes") |> data.frame())$type
        names(nodeN) <- names(V(coGraph))
        V(coGraph)$nodeCat <- nodeN

        if (tag!="none") {
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
            nodeDf <- coGraph |> activate("nodes") |> data.frame()
            V(coGraph)$name <- apply(nodeDf,
                  1,
                  function(x) {ifelse(x["type"]=="Words",
                    ifelse(is.na(pdic[x["name"]]), x["name"], pdic[x["name"]]),
                    x["name"])})
        }

        ## Make uppercase (after preserve)
        nodeName <- V(coGraph)$name
        # dtmCol <- colnames(DTM)
        for (i in madeUpper) {
            # dtmCol[dtmCol == i] <- toupper(i)
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        # colnames(DTM) <- dtmCol

        if (!is.tbl_graph(coGraph)) {
            ret@igraph <- coGraph
        } else {
            ret@igraph <- as.igraph(coGraph)
        }
        E(coGraph)$weightLabel <- round(E(coGraph)$weight, 3)
        
        ## Main plot (probably these functions should be moved to plotNet())
        netPlot <- ggraph(coGraph, layout=layout)

        netPlot <- appendEdges(netPlot, FALSE, edgeLink,
            edgeLabel, showLegend, fontFamily)


        if (tag!="none") { ## Obtain tag coloring
            if (is.null(tagPalette)) {
                cols <- V(coGraph)$tag |> unique()
                if (length(cols)>2) {
                    tagPalette <- RColorBrewer::brewer.pal(8, "Dark2")
                    tagPalette <- colorRampPalette(tagPalette)(length(cols))
                } else {
                    tagPalette <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(length(cols))]
                }
                names(tagPalette) <- cols
                tagPalette["Genes"] <- geneColor
            }
        }

        if (is.null(catColors)) { ## Obtain category coloring
            catLen <- length(unique(V(coGraph)$nodeCat))
            if (catLen>2) {
                catColors <- RColorBrewer::brewer.pal(catLen, "Dark2")
            } else {
                catColors <- RColorBrewer::brewer.pal(3,"Dark2")[seq_len(catLen)]
            }
            names(catColors) <- unique(V(coGraph)$nodeCat)
            catColors["Genes"] <- geneColor
        }

        netPlot <- appendNodesAndTexts(netPlot,tag,colorize,tagPalette,
                          showLegend,catColors,pal,fontFamily,colorText,
                          scaleRange,
                          useSeed, ret, tagColors=tagPalette,
                          discreteColorWord=discreteColorWord)
        if (autoScale) {
            scaleRange <- c((500 * (1 / numWords))/2.5,
                500 * (1 / numWords))
        }
        netPlot <- netPlot +
            scale_size(range=scaleRange, name="Frequency")+
            scale_edge_width(range=scaleEdgeWidth, name = "Correlation")+
            scale_edge_color_gradient(low=pal[1],high=pal[2],
                name = "Correlation", na.value=naEdgeColor)+
            theme_graph()

        if (!is.null(genePathPlot)) {
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
        returnDf <- data.frame(word = names(matSorted),freq=matSorted) |>
            na.omit()
        wcCol <- NULL
        returnDf$wcColor <- "black"
                
        if (genePlot) {
            ## We need DTM with gene ID with row.names
            DTM <- t(as.matrix(docs))
            if (!is.null(mergeCorpus)) {
                stop("Cannot perform genePlot when merging corpus")
            }
            revID <- AnnotationDbi::select(orgDb,
                keys = as.character(fil$Gene_ID), 
                columns = c("SYMBOL"),
                keytype = "ENTREZID")$SYMBOL
            row.names(DTM) <- revID
            
            genemap <- NULL
            nodeName <- returnDf$word
            for (rn in nodeName){
                tmp <- DTM[, rn]
                for (nm in names(tmp[tmp!=0])){
                    genemap <- rbind(genemap, c(rn, nm))
                }
            }
            
            gcnt <- table(genemap[,2])
            gcnt <- gcnt[order(gcnt, decreasing=TRUE)]
            if (length(gcnt)!=1) {ret@geneCount <- gcnt}
            ret@geneMap <- genemap
            genemap <- data.frame(genemap) |> `colnames<-`(c("word","gene"))
            collapsed_genemap <- genemap %>%
                group_by(.data$word) %>%
                summarise(gene_name=paste0(.data$gene, collapse=","))
            returnDf <- merge(returnDf, collapsed_genemap, by="word")
        }

        if (tag!="none") {
            freqWords <- names(matSorted)
            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            if (tagWhole){
                pvc <- pvclust(as.matrix(dist(as.matrix(docs))), method.dist=tag, parallel=cl)
            } else {
                pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))),  method.dist=tag, parallel=cl)
            }
            pvcl <- pvpick(pvc)
            ret@pvclust <- pvc
            ret@pvpick <- pvcl

            wcCol <- returnDf$word
            if (is.null(tagPalette)) {
                tagPalette <- colorRampPalette(brewer.pal(11, "RdBu"))(length(pvcl$clusters |> unique()))
                names(tagPalette) <- pvcl$clusters |> unique()
            }
            for (i in seq_along(pvcl$clusters)){
                for (j in pvcl$clusters[[i]])
                    wcCol[wcCol==j] <- tagPalette[i]
            }
            wcCol[!wcCol %in% tagPalette] <- "grey"
            returnDf$wcColor <- wcCol
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
            returnDf$freq <- returnDf$freq*scaleFreq
        } else {
            showFreq <- returnDf$freq
        }
        
        if (!("min.freq" %in% names(argList))) {
            argList[["min.freq"]] <- 3
        }
        ret@freqDf <- returnDf

        returnDf <- returnDf[returnDf$freq > argList[["min.freq"]], ]
        

        if (tag!="none"){
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- returnDf$freq
            argList[["family"]] <- fontFamily
            argList[["colors"]] <- returnDf$wcColor
            argList[["random.order"]] <- FALSE
            argList[["ordered.colors"]] <- TRUE
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                if (genePlot) {
                    argList[["label_content"]] <- 
                    sprintf("%s<span style='font-size:7.5pt'><br>(%s)</span>",
                        returnDf$word, returnDf$gene_name)
                }
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill="transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        } else {
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- returnDf$freq
            argList[["family"]] <- fontFamily
            if ("bg.color" %in% names(argList)) {
                argList[["bg.colour"]] <- argList[["bg.color"]]
            }
            if (useggwordcloud) {
                if (genePlot) {
                    argList[["label_content"]] <- 
                    sprintf("%s<span style='font-size:7.5pt'><br>(%s)</span>",
                        returnDf$word, returnDf$gene_name)
                }
                wc <- do.call(ggwordcloud::ggwordcloud, argList)+
                scale_size_area(max_size = wcScale)+
                theme(plot.background = element_rect(fill = "transparent",
                    colour = NA))
            } else {
                wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
            }
        }
        ret@wc <- wc
    }
    return(ret)
}
