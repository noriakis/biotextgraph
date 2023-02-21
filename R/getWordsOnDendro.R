#' plotEigengeneNetworksWithWords
#' 
#' plot an eigengene dendrogram with word annotation
#' 
#' 
#' @param MEs module eigengene data (data.frame of row as sample and col as gene cluster IDs)
#' @param colors named vector of cluster
#' @param numberOfWords number of words to be included in pyramid plots
#' @param geneNumLimit clusters with gene number above this threshold are ignored
#' @param geneVecType type of IDs in `colors`
#' @param border whether to draw border on pyramid plots
#' @param nboot pvclust bootstrap number
#' @param type "words" or "enrich"
#' @param highlight words to highlight
#' @param showType when "enrich", which labels to show
#' @param argList passed to wcGeneSummary()
#' @param textSize text size in pyramid plots
#' @param candidateNodes if NULL, all nodes are investigated
#' @param takeIntersect take intersection or frequent words
#' @param useWC plot wordcloud instead of pyramid plot
#' @param wcScale max_size of wordcloud
#' @param dendPlot type of dendrogram plot
#' @param dhc user-specified dendrogram
#' @param horiz horizontal plot or not
#' @param wcArg argument list, pass to ggwordcloud
#' @param useFunc if not specified, use RefSeq summary data
#' @param useDf data.frame to subset when manual function is specified
#' "query" column will be used to subset
#' @param useWGCNA pad named color vector with prefix "ME"
#' @param spacer spacing for grob
#' 
#' @export
#' @import grid gridExtra
#' @import ggplot2
#' @return plot of dendrogram between module eigengenes with word annotations
#' @examples
#' mod <- returnExample()
#' plotEigengeneNetworksWithWords(mod$MEs, mod$colors)
#' @importFrom pvclust pvclust
#' @importFrom dendextend hang.dendrogram pvclust_show_signif_gradient as.ggdend
plotEigengeneNetworksWithWords <- function (MEs, colors, nboot=100,
                                            numberOfWords=10, geneNumLimit=1000,
                                            geneVecType="ENSEMBL", useWC=FALSE,
                                            border=TRUE, type="words", wcScale=3,
                                            candidateNodes=NULL, takeIntersect=TRUE,
                                            showType="ID", textSize=3.5, horiz=FALSE,
                                            dendPlot="pvclust", dhc=NULL, wcArg=list(),
                                            useDf=NULL, useWGCNA=TRUE, spacer=0.005,
                                            highlight=NULL, argList=list(), useFunc=NULL) {

    if (is.null(candidateNodes)) {
        candidateNodes <- unique(colors)
        if (useWGCNA) {
            candidateNodes <- paste0("ME",candidateNodes)
        }
    } else {
        # No limit when nodes are specified
        geneNumLimit <- Inf
    }


    
    ## Make dendrogram
    if (dendPlot=="pvclust") {
        ## Perform pvclust on ME data.frame
        result <- pvclust(MEs, method.dist="cor",
                        method.hclust="average", nboot=nboot)    
        dhc <- result |>
            as.dendrogram() |>
            hang.dendrogram()
    }

    ## Make named vector
    if (useWGCNA) {
        geneVec <- paste0("ME", colors)
        names(geneVec) <- names(colors)
    } else {
        geneVec <- colors
    }
    
    ## Get pyramid plot list using the function.
    ## It takes time when geneNumLimit is large.
    grobList <- getWordsOnDendro(dhc, geneVec,
                                 numberOfWords = numberOfWords,
                                 geneNumLimit = geneNumLimit,
                                 geneVecType = geneVecType,
                                 type=type, highlight=highlight,
                                 textSize=textSize, wcArg=wcArg,
                                 candidateNodes=candidateNodes,
                                 showType=showType, takeIntersect=takeIntersect,
                                 argList=argList, useWC=useWC, wcScale=wcScale,
                                 useFunc=useFunc, useDf=useDf)
    
    ## Plot dendrogram ggplot, using the pvclust p-values
    if (dendPlot=="pvclust") {
        dendroPlot <- dhc |> pvclust_show_signif_gradient(result) |> ggplot(horiz=horiz) 
    } else {
        dendroPlot <- dhc |> as.ggdend() |> ggplot(horiz=horiz)
    }
    
    ## Plot the grob on dendrogram using annotation_custom.
    ## If border is TRUE, border line is drawn using grid.rect.
    for (gr in grobList){
        if (border){
            addPlot <- ggplotify::as.grob(function(){
                grid.arrange(gr$plot)
                grid.rect(width = .98, height = .98,
                          gp = gpar(lwd = 1, col = "black", fill = NA))
            })
        } else {
            addPlot <- gr$plot
        }
        if (horiz) {
            dendroPlot <- dendroPlot +
                annotation_custom(addPlot, xmin=gr$xmin, xmax=gr$xmax,
                                  ymin=-1*gr$height-spacer, ymax=-1*gr$heightup+spacer)
        } else {
            dendroPlot <- dendroPlot +      
                annotation_custom(addPlot, xmin=gr$xmin, xmax=gr$xmax,
                                  ymin=gr$height+spacer, ymax=gr$heightup-spacer)
        }
    }
    
    dendroPlot
}




#' getWordsOnDendro
#' 
#' Get grobs to plot on the dendrogram with the position information
#' 
#' @param dhc dendrogram
#' @param geneVec gene-named vector of node names in dendrogram
#' @param geneNumLimit when the gene number is above this threshold,
#'                     pyramid plots are not produced 
#' @param geneVecType type of the name of geneVec (default: ENSEMBL)
#' @param numberOfWords the number of words to plot (default: 25)
#' @param type "words" or "enrich"
#' @param showType when "enrich", which labels to show
#' @param highlight words to highlight
#' @param argList passed to wcGeneSummary
#' @param textSize text size in pyramid plots
#' @param candidateNodes if NULL, all nodes are investigated
#' @param takeIntersect take intersection or frequent words
#' @param useWC use wordcloud
#' @param wcScale max_size of wordcloud
#' @param wcArg argument list for ggwordcloud
#' @param useFunc function to summarize text
#' @param useDf data.frame to subset when manual function is specified
#' @return list of pyramid plot grobs and its positions
#' @import tm
#' @import org.Hs.eg.db
#' @importFrom ggdendro dendro_data
#' @importFrom dplyr select
#' @importFrom ggdendro label
#' @importFrom dendextend get_nodes_attr get_subdendrograms nnodes
#' 
#' @examples
#' mod <- returnExample()
#' result <- hclust(dist(t(mod$MEs)))
#' dhc <- result |>
#'     as.dendrogram() |>
#'     dendextend::hang.dendrogram()
#' geneVec <- paste0("ME", mod$colors)
#' names(geneVec) <- names(mod$colors)
#' getWordsOnDendro(dhc, geneVec, numberOfWords = 2)
#' @export
#' 
getWordsOnDendro <- function(dhc, geneVec, geneNumLimit=1000,
                            geneVecType="ENSEMBL", useWC=FALSE,
                            numberOfWords=25, showType="ID", wcArg=list(),
                            highlight=NULL, textSize=3.5, wcScale=3,
                            candidateNodes=NULL, takeIntersect=TRUE, useDf=NULL,
                            type="words", argList=list(), useFunc=NULL) {
    
    ## Filter high frequency words if needed
    # filterWords <- allFreqGeneSummary[
    #                         allFreqGeneSummary$freq > excludeFreq,]$word
    # if (!is.na(excludeTfIdf)){
    #     filterWords <- c(filterWords,
    #         allTfIdfGeneSummary[
    #                         allTfIdfGeneSummary$tfidf > excludeTfIdf,]$word)
    # }
    # filterWords <- c(filterWords, "pmids", "geneid") ## Excluded by default
    # qqcat("filtered @{length(filterWords)} words (frequency | tfidf)\n")

    
    grobList <- list()
    alHeights <- c(((dhc %>%
        get_subdendrograms(k=1))[[1]] %>%
        get_nodes_attr("height")))
    alNodes <- c(((dhc %>%
        get_subdendrograms(k=1))[[1]] %>%
        get_nodes_attr("label")))
    curHeights <- c(((dhc %>%
        get_subdendrograms(k=1))[[1]] %>%
        get_nodes_attr("height"))[1])
    k <- 2
    grobNum <- 1
    ddata <-dendro_data(dhc)
    labelPos <- ddata$labels

    while (k < length(dhc %>% labels)){
        subdendro <- dhc %>% get_subdendrograms(k=k)
        
        for (i in subdendro){

            NODES <- i %>% get_nodes_attr("label")
            NODES <- NODES[!is.na(NODES)]

            XMIN <- as.numeric(labelPos %>%
                filter(label==NODES[1]) %>% select(.data$x))
            XMAX <- as.numeric(labelPos %>%
                filter(label==NODES[length(NODES)]) %>% select(.data$x))
            HEIGHT <- get_nodes_attr(i, "height")[1]
            
            if (!HEIGHT %in% curHeights){
                if (nnodes(i)!=1){
                    lb <- i %>% dendextend::cutree(k=2) ## not stats::cutree
                    # L <- as.numeric(sapply(strsplit(names(lb[lb==1]),
                    # "ME"), "[", 2)) # WGCNA
                    # R <- as.numeric(sapply(strsplit(names(lb[lb==2]),
                    # "ME"), "[", 2)) # WGCNA
                    if (mean(as.numeric((labelPos %>% filter(
                                    .data$label %in% names(lb[lb==1])) %>%
                                    select(.data$x))[,1])) >
                        mean(as.numeric((labelPos %>% filter(
                                    .data$label %in% names(lb[lb==2])) %>%
                                    select(.data$x))[,1]))) {
                                                        R <- names(lb[lb==1])
                                                        L <- names(lb[lb==2])
                                                    } else {
                                                        L <- names(lb[lb==1])
                                                        R <- names(lb[lb==2])
                    }

                    if (length(names(geneVec)[geneVec %in% L])<geneNumLimit &
                        length(names(geneVec)[geneVec %in% R])<geneNumLimit)
                    {
                        if (sum(L %in% candidateNodes)==length(L) | sum(R %in% candidateNodes)==length(R)) {
                            pyrm <- returnPyramid(L, R, geneVec, geneVecType, highlight=highlight,
                                numberOfWords=numberOfWords, type=type, showType=showType,
                                argList=argList, textSize=textSize, takeIntersect=takeIntersect,
                                useWC=useWC, wcScale=wcScale, wcArg=wcArg, useFunc=useFunc,
                                useDf=useDf)
                            if (!is.null(pyrm)){
                                grobList[[as.character(grobNum)]]$plot <- pyrm
                                grobList[[as.character(grobNum)]]$height <- HEIGHT
                                grobList[[as.character(grobNum)]]$xmin <- XMIN
                                grobList[[as.character(grobNum)]]$xmax <- XMAX

                                hind <- which(alHeights==HEIGHT)
                                revNodes <- rev(alNodes[1:(hind-1)])
                                revHeights <- rev(alHeights[1:(hind-1)])
                                for (tmp in seq_len(length(revHeights))){
                                    if (revHeights[tmp]>HEIGHT){
                                        if (is.na(revNodes[tmp])) {
                                            HEIGHTUP <- revHeights[tmp]
                                            break
                                        }
                                    }
                                }
                                
                                grobList[[as.character(grobNum)]]$heightup <- HEIGHTUP
                                
                                curHeights <- c(curHeights, HEIGHT)
                                grobNum <- grobNum + 1
                            }
                        }
                    }
                    
                }
            }
        }
        k <- k + 1
    }
    return(grobList)
}


#' returnPyramid
#' 
#' Return pyramid plots
#' 
#' 
#' @param L genes in the cluster
#' @param R genes in the other cluster
#' @param geneVec gene-named vector of node names in dendrogram
#' @param geneVecType type of the name of geneVec (default: ENSEMBL)
#' @param numberOfWords the number of words to plot
#' @param widths parameters to pass to patchwork
#' @param lowCol gradient low color
#' @param highCol gradient high color
#' @param type "words" or "enrich"
#' @param showType which labels to show in enrich
#' @param highlight words to highlight
#' @param orgDb organism database to use in enrich
#' @param argList parameters passed to wcGeneSummary()
#' @param wrap wrap the strings
#' @param textSize text size in pyramid plots
#' @param takeIntersect take intersection or frequent words
#' @param useWC return wordcloud
#' @param wcScale if useWC, number of size scaling (max_size)
#' @param wcArg argument list for ggwordcloud
#' @param useFunc function to summarize text
#' @param useDf data.frame to subset when manual function is specified
#' 
#' @return list of pyramid plot grobs and its positions
#' @import tm
#' @import ggplot2
#' @import patchwork
#' @importFrom GeneSummary loadGeneSummary
#' @importFrom dplyr mutate if_else
#' @importFrom ggdendro dendro_data
#' @importFrom dendextend get_nodes_attr get_subdendrograms get_nodes_attr nnodes
#' 
#' 
returnPyramid <- function(L, R, geneVec, geneVecType,
                        numberOfWords=25, widths=c(0.3,0.6,0.3),
                        lowCol="blue", showType="ID",
                        highCol="red", highlight=NULL, wcScale=3,
                        type="words", wrap=15, textSize=3.5,
                        takeIntersect=TRUE, useWC=FALSE, wcArg=list(),
                        orgDb=org.Hs.eg.db, argList=list(), useFunc=NULL, useDf=NULL) {
    ## Convert to ENTREZ ID
    # geneList <- AnnotationDbi::select(orgDb,
    #     keys = names(geneVec)[geneVec %in% L],
    #     columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    # filL <- tb %>% filter(Gene_ID %in% geneList)
    # filL <- filL[!duplicated(filL$Gene_ID),]

    # geneList <- AnnotationDbi::select(orgDb,
    #     keys = names(geneVec)[geneVec %in% R],
    #     columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    # filR <- tb %>% filter(Gene_ID %in% geneList)
    # filR <- filR[!duplicated(filR$Gene_ID),]
    
    # all_L <- paste(filL$Gene_summary, collapse = " ")
    # all_R <- paste(filR$Gene_summary, collapse = " ")
    # all_bet <- c(all_L, all_R)
    
    # all_bet <- VectorSource(all_bet)
    # all_corpus <- VCorpus(all_bet)


    if (type=="words") {

        if (useWC) {

            if (!is.null(useFunc)) {
                if (!is.null(useDf)) {
                    sch <- c(names(geneVec)[geneVec %in% L], names(geneVec)[geneVec %in% R])
                    inputMan <- subset(useDf, useDf$query %in% sch)
                    argList[["df"]] <- inputMan
                    retWC <- do.call(useFunc, argList)
                }
            } else {
                argList[["geneList"]] <- c(names(geneVec)[geneVec %in% L],
                names(geneVec)[geneVec %in% R])
                argList[["keyType"]] <- geneVecType
                retWC <- do.call("wcGeneSummary", argList)
            }

            if (length(wcArg)==0) {
                wcArg[["min.freq"]] <- 1
                wcArg[["max.words"]] <- Inf
                wcArg[["rot.per"]] <- 0.5
                wcArg[["random.order"]] <- FALSE
                wcArg[["colors"]] <- brewer.pal(10, sample(row.names(RColorBrewer::brewer.pal.info), 1))
            }

            wcArg[["words"]] <- retWC@freqDf$word
            wcArg[["freq"]] <- retWC@freqDf$freq



            plt <- do.call(ggwordcloud::ggwordcloud, wcArg)+
                 scale_size_area(max_size = wcScale)+
                 theme(plot.background = element_rect(fill = "white",colour = NA))
            return(plt)
        }

        if (is.null(useFunc)) {
            argList[["geneList"]] <- names(geneVec)[geneVec %in% L]
            argList[["collapse"]] <- TRUE
            argList[["onlyTDM"]] <- TRUE
            argList[["keyType"]] <- geneVecType
            all_L <- as.matrix(do.call("wcGeneSummary",argList))
            
            argList[["geneList"]] <- names(geneVec)[geneVec %in% R]
            all_R <- as.matrix(do.call("wcGeneSummary",argList))
        } else {
            if (!is.null(useDf)) {
                argList[["df"]] <- subset(useDf, useDf$query %in% names(geneVec)[geneVec %in% L])
                argList[["collapse"]] <- TRUE
                argList[["onlyTDM"]] <- TRUE
                all_L <- as.matrix(do.call(useFunc,argList))
                
                argList[["df"]] <- subset(useDf, useDf$query %in% names(geneVec)[geneVec %in% R])
                all_R <- as.matrix(do.call(useFunc,argList))
            }
        }

        if (takeIntersect) {

            common <- intersect(row.names(all_L), row.names(all_R))
            all_m <- cbind(all_L[common, ], all_R[common, ])
            # ## Clean the corpus
            # all_corpus <- makeCorpus(all_corpus, filterWords, additionalRemove, numOnly, stem)
            # if (tfidf){
            #     stop("Use of tfidf on returnPyramid is currently not supported")
            #     all_tdm <- TermDocumentMatrix(all_corpus, list(weighting = weightTfIdf))
            # } else {
            #     all_tdm <- TermDocumentMatrix(all_corpus)
            # }
            # all_m <- as.matrix(all_tdm)
            
            ## Modified from: https://rpubs.com/williamsurles/316682
            common_words <- subset(all_m, all_m[, 1] > 0 & all_m[, 2] > 0)

            difference <- abs(common_words[, 1] - common_words[, 2])
            common_words <- cbind(common_words, difference)
            common_words <- common_words[order(common_words[, 3], decreasing = TRUE), ]
            if (dim(common_words)[1]<numberOfWords){
               numberOfWords <- dim(common_words)[1]
            }
        } else {
            sorted_L <- all_L[order(all_L[,1], decreasing = TRUE),] 
            sorted_R <- all_R[order(all_R[,1], decreasing = TRUE),]
        }
        
        L <- paste0(L, collapse="_")
        R <- paste0(R, collapse="_")
        
        # Referencing:
        # https://stackoverflow.com/questions/54191369/
        # how-to-create-pyramid-bar-chart-in-r-with-y-axis-labels-between-the-bars

        if (numberOfWords==0){
            qqcat("No common words found\n")
            return(NULL)
        } else {
            if (takeIntersect) {
                topDf <- data.frame(
                    ME = factor(c(rep(L,numberOfWords),
                        rep(R,numberOfWords)),
                        levels = c(L,R)),
                    Word = c(common_words[1:numberOfWords, 1],
                        common_words[1:numberOfWords, 2]),
                    Label = c(rownames(common_words[1:numberOfWords, ]),
                        rownames(common_words[1:numberOfWords, ]))
                )
            } else {
               topDf <- data.frame(
                    ME = factor(c(rep(L,numberOfWords),
                                    rep(R,numberOfWords)),
                                    levels = c(L,R)),
                    Word = c(sorted_L[1:numberOfWords],
                             sorted_R[1:numberOfWords]),
                    Label = c(names(sorted_L[1:numberOfWords]),
                            names(sorted_R[1:numberOfWords]))
                    )
               topDf$Label <- factor(topDf$Label, levels = topDf$Label)
            }

            ## Convert to uppercase
            nodeName <- topDf$Label   
            # for (i in madeUpper) {
            #     nodeName[nodeName == i] <- toupper(i)
            # }
            topDf$Label <- nodeName
            
            if (takeIntersect) {
                gg1 <- topDf %>%
                    mutate(Count = if_else(.data$ME == L, .data$Word, 0)) %>%
                    ggplot(aes( .data$Label, .data$Count, fill = .data$Count)) +
                    geom_col(width = 0.6) +
                    coord_flip() +
                    scale_y_reverse()+
                    scale_fill_gradient(low=lowCol,high=highCol)+
                    # scale_fill_manual(values = c("Red", "Blue")) +
                    theme_void()+
                    theme(
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position = "none")
                gg2 <- topDf %>%
                    filter(.data$ME == L) %>%
                    ggplot(aes(.data$Label, 0, label = .data$Label)) +
                    geom_text(size=textSize) +
                    coord_flip() +
                    theme_void()

                gg3 <- topDf %>%
                mutate(Count = if_else(.data$ME == R, .data$Word, 0)) %>%
                ggplot(aes( .data$Label, .data$Count, fill = .data$Count)) +
                geom_col(width = 0.6) +
                coord_flip() +
                theme_void() +
                scale_fill_gradient(low=lowCol,high=highCol)+
                theme(
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    legend.position = "none")

            } else {
                ## TODO: some re-ordering

                gg1 <- topDf %>%
                    mutate(Count = if_else(.data$ME == L, .data$Word, 0)) %>%
                    ggplot(aes( .data$Label, .data$Count, fill = .data$Count)) +
                    geom_col(width = 0.6) +
                    coord_flip() +
                    scale_y_reverse()+
                    scale_fill_gradient(low=lowCol,high=highCol)+
                    # scale_fill_manual(values = c("Red", "Blue")) +
                    theme_void()+
                    theme(
                        axis.title.y=element_blank(),
                        axis.text.y=element_blank(),
                        axis.ticks.y=element_blank(),
                        legend.position = "none")

                gg2 <- topDf %>%
                    ggplot(aes(.data$Label, 0, label = .data$Label)) +
                    geom_text(size=textSize) +
                    coord_flip() +
                    theme_void()

                gg3 <- topDf %>%
                mutate(Count = if_else(.data$ME == R, .data$Word, 0)) %>%
                ggplot(aes( .data$Label, .data$Count, fill = .data$Count)) +
                geom_col(width = 0.6) +
                coord_flip() +
                theme_void() +
                scale_fill_gradient(low=lowCol,high=highCol)+
                theme(
                    axis.title.y=element_blank(),
                    axis.text.y=element_blank(),
                    axis.ticks.y=element_blank(),
                    legend.position = "none")
            }

            pyramid <- gg1+gg2+gg3+plot_layout(widths=widths)
            pyramidGrob <- patchworkGrob(pyramid)
            return(pyramidGrob)
        }
    } else {
        ## Use KEGG
        LI <- names(geneVec)[geneVec %in% L]
        RI <- names(geneVec)[geneVec %in% R]
        if (geneVecType!="ENTREZID") {
            LI <- clusterProfiler::bitr(LI, fromType=geneVecType, toType="ENTREZID", OrgDb=orgDb)$ENTREZID
            RI <- clusterProfiler::bitr(RI, fromType=geneVecType, toType="ENTREZID", OrgDb=orgDb)$ENTREZID
        }
        leftE <- clusterProfiler::enrichKEGG(LI)
        rightE <- clusterProfiler::enrichKEGG(RI)

        lSig <- leftE@result[1:numberOfWords,]
        rSig <- rightE@result[1:numberOfWords,]

        lSig <- lSig[!is.na(lSig$ID),]
        rSig <- rSig[!is.na(rSig$ID),]
        
        if (dim(lSig)[1]==0 & dim(rSig)[1]==0) {
          stop("No enriched term found for these clusters")
          return(NULL)
        }
        lSig$value <- -1 * -log10(lSig$p.adjust)
        rSig$value <- -log10(rSig$p.adjust)

        if (showType=="Description") {
            if (!is.null(wrap)) {
              lSig$plotID <- stringr::str_wrap(lSig$Description, wrap)
              rSig$plotID <- stringr::str_wrap(rSig$Description, wrap)
            }            
        } else if (showType=="ID") {
            lSig$plotID <- lSig$ID
            rSig$plotID <- rSig$ID
        } else {
            stop("showType must be ID or Description")
        }

        if (!is.null(highlight)) {
            sigCol <- NULL
            for (id in lSig$plotID) {
                if (id %in% highlight) {
                    sigCol <- c(sigCol, "highlight")
                } else {
                    sigCol <- c(sigCol, "not_highlight")
                }
            }
            lSig$plotCol <- sigCol

            sigCol <- NULL
            for (id in rSig$plotID) {
                if (id %in% highlight) {
                    sigCol <- c(sigCol, "highlight")
                } else {
                    sigCol <- c(sigCol, "not_highlight")
                }
            }
            rSig$plotCol <- sigCol
            highlightCol <- c("black","red")
            names(highlightCol) <- c("not_highlight","highlight")
        } else {
            lSig$plotCol <- rep("not_highlight",length(lSig$plotID))
            rSig$plotCol <- rep("not_highlight",length(rSig$plotID))
            highlightCol <- c("black")
            names(highlightCol) <- c("not_highlight")

        }

        gg1 <- lSig %>%
          ggplot(aes( .data$plotID, .data$value, fill = .data$value)) +
          geom_col(width = 0.6) +
          coord_flip() +
          scale_fill_gradient(low=lowCol,high=highCol)+
          # scale_fill_manual(values = c("Red", "Blue")) +
          theme_void()+
          theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none")

        gg2 <- lSig %>%
          ggplot(aes(.data$plotID, 0, label = .data$plotID, color=.data$plotCol)) +
          geom_text(size=textSize) +
          scale_color_manual(values=highlightCol, guide="none")+
          coord_flip() +
          theme_void()

        gg3 <- rSig %>%
          ggplot(aes( .data$plotID, .data$value, fill = .data$value)) +
          geom_col(width = 0.6) +
          coord_flip() +
          theme_void() +
          scale_fill_gradient(low=lowCol,high=highCol)+
          theme(
            axis.title.y=element_blank(),
            axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            legend.position = "none")

        gg4 <- rSig %>%
          ggplot(aes(.data$plotID, 0, label = .data$plotID, color=.data$plotCol)) +
          geom_text(size=textSize) +
          scale_color_manual(values=highlightCol, guide="none")+
          coord_flip() +
          theme_void()

        areas <- "
        AA#
        AA#
        #BB
        #BB
        "
        pyramid <- (gg1 + gg2) / (gg4 + gg3) + plot_layout(design=areas)
        pyramidGrob <- patchworkGrob(pyramid)
        return(pyramidGrob)
    }
}
