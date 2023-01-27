#' wcGeneSummary
#' 
#' Plot wordcloud of RefSeq description obtained by GeneSummary
#' Testing function for new class
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
#' @param onWholeTDM calculate correlation network
#'                   on whole dataset or top-words specified by numWords
#' @param verb show verb on edges (using udpipe)
#'             peform udpipe_download_model(language = "english") beforehand,
#'             and path udpipe file to udpipeModel
#' @param verbPOS Part of Speech tags in verbs, like "VBN", "VBZ"
#' @param udpipeModel udpipe model file name
#' @param cl for parPvclust, parallel clustering can be performed
#' @param stem whether to use stemming
#' @param nodePal node palette when tag is TRUE
#' @param preserve preserve original characters
#' @param takeMax take max values for each term in term-document matrix
#' @param filterMax use pre-calculated filter based on max-values when excluding TfIdf
#' Otherwise take sum.
#' @param collapse default to FALSE, collapse all the sentences
#' @param argList parameters to pass to wordcloud()
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
#' geneList <- c("DDX41")
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
                            pvclAlpha=0.95, bn=FALSE, R=20, cl=FALSE,
                            ngram=NA, plotType="wc", onlyTDM=FALSE, stem=FALSE,
                            colorText=FALSE, corThresh=0.2, genePlot=FALSE,
                            genePathPlot=NA, genePathPlotSig=0.05, tag=FALSE,
                            layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, 
                            enrich=NULL, topPath=10, ora=FALSE, tagWhole=FALSE,
                            mergeCorpus=NULL, numOnly=TRUE, madeUpperGenes=TRUE,
                            onWholeTDM=FALSE, pre=TRUE, verbPOS=c("VBZ"),
                            nodePal=palette(), collapse=FALSE,
                            verb=FALSE, udpipeModel="english-ewt-ud-2.5-191206.udpipe",
                            argList=list()) {
    ret <- new("osplot")
    ret@query <- geneList
    ret@type <- "refseq"

    if (verb) {
        qqcat("Using verb mode\n")
        if (bn) {stop("verb can be only used with undirected graph")}
        edgeLabel <- TRUE
        udmodel_english <- udpipe::udpipe_load_model(file = udpipeModel)}
    if (madeUpperGenes){
        madeUpper <- c(madeUpper, tolower(keys(orgDb, "SYMBOL")))
    }
    # if (ora & !numOnly) {
    #     stop("ora should be used with numOnly=TRUE, as the background is calculated based on numOnly=TRUE")
    # }
    # if (ora & stem) {
    #     stop("ora should be used with stem=FALSE, as the background is calculated based on stem=FALSE")
    # }

    if (is.null(mergeCorpus)) {
        qqcat("Input genes: @{length(geneList)}\n")
        if (keyType!="ENTREZID"){
            geneList <- AnnotationDbi::select(orgDb,
                keys = geneList, columns = c("ENTREZID"),
                keytype = keyType)$ENTREZID
            geneList <- geneList[!is.na(geneList)]
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

        ## If specified pathway option
        if (!is.null(enrich)) {
            if (genePlot) {stop("genePlot can't be performed in enrichment analysis mode")}
            qqcat("Performing enrichment analysis")
            if (enrich=="reactome"){
                pathRes <- enrichPathway(geneList)
                pathRes@result$Description <- gsub("Homo sapiens\r: ",
                                "",
                                pathRes@result$Description)                
            } else if (enrich=="kegg"){
                pathRes <- enrichKEGG(geneList)
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
        } else {
            ## Load from GeneSummary
            ## Already performed, and automatically loaded
            # load("allFreqGeneSummary.rda") 
            tb <- loadGeneSummary(organism = organism)
            fil <- tb %>% filter(Gene_ID %in% geneList)
            fil <- fil[!duplicated(fil$Gene_ID),]
            ret@rawText <- fil

            if (ora){
                qqcat("Performing ORA\n")
                sig <- textORA(geneList)
                sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>0.05]
                qqcat("Filtered @{length(sigFilter)} words (ORA)\n")
                filterWords <- c(filterWords, sigFilter)
                ret@ora <- sig
                # if (oraPlot) {
                #     ret <- returnORAPlot(ret)
                # }
            }

            if (verb) {
                ## Annotate verbs
                s <- udpipe::udpipe_annotate(udmodel_english, fil$Gene_summary)
                x <- data.frame(s)
                x2 <- x |> dplyr::filter(upos=="VERB")
                x3 <- x2[x2$xpos %in% verbPOS,]
                verbs <- tolower(unique(x3$token))
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

    if (takeMax) {
        perterm <- apply(mat, 1, max, na.rm=TRUE)
    } else {
        perterm <- rowSums(mat)
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
            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
            returnDf[returnDf$word == i,"word"] <- toupper(i)
        }
        ## Do we need to return converted freqDf? 
        ## TODO: Force all functions to return lower freqDf.
        ret@freqDf <- returnDf

        # freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
        ## TODO: before or after?
        freqWordsTDM <- t(as.matrix(docs))
        
        if (tag) {
            ## TODO: tagging based on cluster_walktrap
            ## TODO: tag based on all words (currently top words for computational time)
            if (tagWhole){
                pvc <- pvclust(as.matrix(dist(t(freqWordsTDM))), parallel=cl)
            } else {
            pvc <- pvclust(as.matrix(dist(
                t(
                    freqWordsTDM[,colnames(freqWordsTDM) %in% freqWords]
                    )
                )), parallel=cl)
            }
            pvcl <- pvpick(pvc, alpha=pvclAlpha)
            ret@pvclust <- pvc
            ret@pvpick <- pvcl
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
            row.names(freqWordsTDM) <- revID
        }

        ## genePathPlot: plot associated genes and pathways
        if (!is.na(genePathPlot)) {
            
            if (genePathPlot == "reactome") {
                pathRes <- enrichPathway(geneList)
                pathRes@result$Description <- gsub("Homo sapiens\r: ",
                                "",
                                pathRes@result$Description)
            }
            else if (genePathPlot == "kegg") {
                pathRes <- enrichKEGG(geneList)
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
            qqcat("bn specified, R=@{R}\n")
            # To avoid computaitonal time, subset to numWords
            bnboot <- boot.strength(
                data.frame(freqWordsTDM[,colnames(freqWordsTDM) %in% freqWords]),
                algorithm = "hc", R=R)
            ret@strength <- bnboot
            av <- averaged.network(bnboot)
            avig <- bnlearn::as.igraph(av)
            el <- data.frame(as_edgelist(avig))
            colnames(el) <- c("from","to")
            mgd <- merge(el, bnboot, by=c("from","to"))
            colnames(mgd) <- c("from","to","weight","direction")
            coGraph <- graph_from_data_frame(mgd, directed=TRUE)
        } else {
            ## Check correlation
            ## TODO: speed up calculation using Rcpp 
            if (onWholeTDM){
                corData <- cor(freqWordsTDM)
            } else {
                corData <- cor(freqWordsTDM[,colnames(freqWordsTDM) %in% freqWords])
            }
            ret@corMat <- corData

            ## Set correlation below threshold to zero
            corData[corData<corThresh] <- 0
            coGraph <- graph.adjacency(corData, weighted=TRUE,
                        mode="undirected", diag = FALSE)
        }

        if (verb) {
            eg <- as_data_frame(coGraph)
            newegconc <- c()
            newatt <- c()
            nmls <- list()
            for (j in seq_len(dim(eg)[1])){
                tmp <- as.character(eg[j,1:2])
                if (length(intersect(verbs, tmp))==0){
                    newegconc <- rbind(newegconc,c(tmp,eg[j,3]))
                    newatt <- c(newatt,NA)
                } else if (length(intersect(verbs, tmp))==1){
                    ve <- tmp[tmp%in%verbs]
                    nve <- tmp[!tmp%in%verbs]
                    nmls[[ve]] <- c(nmls[[ve]],nve)
                } else {
                    next
                }
            }
            for (i in names(nmls)){
                if (length(nmls[[i]])!=1) {
                    cmbn <- t(combn(nmls[[i]],2))
                    cmbn <- cbind(cmbn,rep(NA,dim(cmbn)[1]))
                    newegconc <- rbind(newegconc,cmbn)
                    newatt <- c(newatt,rep(i,dim(cmbn)[1]))
                }
            }
            ## Restore cor values
            for (i in seq_len(nrow(newegconc))){
              tmp <- newegconc[i,]
              if (is.na(tmp[3])){
                y <- try(corData[tmp[1],tmp[2]])
                if (class(y)=="try-error"){
                  y <- NA
                }
                tmp[3] <- NA
              }
              newegconc[i,] <- tmp
            }
            ## Merge
            verbDf <- data.frame(cbind(newegconc, newatt))
            colnames(verbDf) <- c("from","to","Weight","Verb")
            coGraph <- graph_from_data_frame(verbDf,
                directed = FALSE)
            coGraph <- simplify(coGraph,
                remove.multiple = TRUE,
                edge.attr.comb = "concat")
            E(coGraph)$verb <- vapply(E(coGraph)$Verb,
                function(x) paste0(x[!is.na(x)],collapse=","),
                FUN.VALUE = "character")
            E(coGraph)$weight <- unlist(sapply(E(coGraph)$Weight,
                function(x) as.numeric(x[!is.na(x)])[1]))
        }

        ## before or after?
        coGraph <- induced.subgraph(coGraph, names(V(coGraph)) %in% freqWords)
        V(coGraph)$Freq <- matSorted[V(coGraph)$name]

        if (deleteZeroDeg){
            coGraph <- induced.subgraph(coGraph, degree(coGraph) > 0)
        }

        nodeName <- V(coGraph)$name
        tdmCol <- colnames(freqWordsTDM)
        for (i in madeUpper) {
            tdmCol[tdmCol == i] <- toupper(i)
            nodeName[nodeName == i] <- toupper(i)
        }
        V(coGraph)$name <- nodeName
        colnames(freqWordsTDM) <- tdmCol

        if (genePlot) {
            genemap <- c()
            for (rn in nodeName){
                tmp <- freqWordsTDM[ ,rn]
                for (nm in names(tmp[tmp!=0])){
                    genemap <- rbind(genemap, c(rn, nm))
                }
            }

            gcnt <- table(genemap[,2])
            gcnt <- gcnt[order(gcnt, decreasing=TRUE)]
            ret@geneCount <- gcnt
            
            incGene <- names(gcnt)[1:genePlotNum]
            genemap <- genemap[genemap[,2] %in% incGene,]
            ret@geneMap <- genemap
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
                    if (verb) {
                        netPlot <- netPlot +
                                    geom_edge_link(
                                        aes(width=weight,
                                        color=weight,
                                        label=verb),
                                        angle_calc = 'along',
                                        label_dodge = unit(2.5, 'mm'),
                                        alpha=0.5,
                                        show.legend = showLegend)
                    } else {
                        netPlot <- netPlot +
                                    geom_edge_link(
                                        aes(width=weight,
                                        color=weight,
                                        label=round(weight,3)),
                                        angle_calc = 'along',
                                        label_dodge = unit(2.5, 'mm'),
                                        alpha=0.5,
                                        show.legend = showLegend)                        
                    }
                } else {
                    netPlot <- netPlot +
                                geom_edge_link(aes(width=weight, color=weight),
                                    alpha=0.5, show.legend = showLegend)
                }
            } else {
                if (edgeLabel){
                    if (verb) {
                        netPlot <- netPlot +
                                    geom_edge_diagonal(
                                        aes(width=weight,
                                        color=weight,
                                        label=verb),
                                        angle_calc = 'along',
                                        label_dodge = unit(2.5, 'mm'),
                                        alpha=0.5,
                                        show.legend = showLegend)
                    } else {
                        netPlot <- netPlot +
                                    geom_edge_diagonal(
                                        aes(width=weight,
                                        color=weight,
                                        label=round(weight,3)),
                                        angle_calc = 'along',
                                        label_dodge = unit(2.5, 'mm'),
                                        alpha=0.5,
                                        show.legend = showLegend)                        
                    }
                } else {
                    netPlot <- netPlot +
                                geom_edge_diagonal(aes(width=weight, color=weight),
                                    alpha=0.5, show.legend = showLegend)                
                }
            }
        }

        if (tag) {
            netPlot <- netPlot + geom_node_point(aes(size=Freq, color=tag),
                                                show.legend = showLegend) +
            scale_color_manual(values=nodePal)
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
                        geom_node_text(aes(label=name, size=Freq),
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
        ret@net <- netPlot
    } else {
        ## WC

        returnDf <- data.frame(word = names(matSorted),freq=matSorted)

        if (tag) {
            freqWords <- names(matSorted)
            freqWordsTDM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
            if (tagWhole){
                pvc <- pvclust(as.matrix(dist(as.matrix(docs))), parallel=cl)
            } else {
                pvc <- pvclust(as.matrix(dist(t(freqWordsTDM))), parallel=cl)
            }
            pvcl <- pvpick(pvc)
            ret@pvclust <- pvc
            ret@pvpick <- pvcl

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
        if (preserve) {
            for (nm in unique(returnDf$word)) {
                if (nm %in% names(pdic)) {
                    returnDf[returnDf$word == nm, "word"] <- pdic[nm]
                }
            }
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
            argList[["words"]] <- returnDf$word
            argList[["freq"]] <- showFreq
            wc <- as.ggplot(as_grob(~do.call("wordcloud", argList)))
        }
        ret@freqDf <- returnDf
        ret@wc <- wc
    }

    return(ret)
}

#' osp
#' 
#' alias for wcGeneSummary
#' 
#' @examples osp("DDX41")
#' @export
osp <- wcGeneSummary