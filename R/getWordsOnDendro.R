#' plotEigengeneNetworksWithWords
#' 
#' plot an eigengene dendrogram with word annotation
#' 
#' 
#' @param MEs module eigengene data
#' @param colors named vector of cluster
#' @param numberOfWords number of words to be included in pyramid plots
#' @param geneNumLimit clusters with gene number above this threshold are ignored
#' @param geneVecType type of IDs in `colors`
#' @param excludeFreq exclude the words above this threshold
#' @param excludeTfIdf exclude the words above this threshold
#' @param tfidf use tfidf when making TDM
#' @param border whether to draw border on pyramid plots
#' @param madeUpper words with uppercase
#' @param nboot pvclust bootstrap number
#' 
#' @export
#' @import grid gridExtra
#' @import ggplot2
#' @return plot of dendrogram between module eigengenes with word annotations
#' @examples
#' mod <- returnExample()
#' plotEigengeneNetworksWithWords(mod$MEs, mod$colors)
#' @importFrom pvclust pvclust
#' @importFrom dendextend hang.dendrogram pvclust_show_signif_gradient
plotEigengeneNetworksWithWords <- function (MEs, colors, nboot=100,
                                            numberOfWords=10, geneNumLimit=1000,
                                            geneVecType="ENSEMBL",
                                            excludeFreq=5000, excludeTfIdf=NA,
                                            border=TRUE, tfidf=FALSE,
                                            madeUpper=c("rna","dna")) {
    ## Perform pvclust on ME data.frame
    result <- pvclust(MEs, method.dist="cor",
        method.hclust="average", nboot=nboot)
    
    # make dendrogram    
    dhc <- result |>
        as.dendrogram() |>
        hang.dendrogram()

    ## Make named vector
    geneVec <- paste0("ME", colors)
    names(geneVec) <- names(colors)
    
    ## Get pyramid plot list using the function.
    ## It takes time when geneNumLimit is large.
    grobList <- getWordsOnDendro(dhc, geneVec, numberOfWords = numberOfWords,
                                 geneNumLimit = geneNumLimit,
                                 geneVecType = geneVecType,
                                 excludeFreq = excludeFreq, excludeTfIdf=excludeTfIdf,
                                 tfidf = tfidf,
                                 madeUpper = madeUpper)
    
    ## Plot dendrogram ggplot, using the pvclust p-values
    dendroPlot <- dhc |> pvclust_show_signif_gradient(result) |> ggplot() 
    
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
        dendroPlot <- dendroPlot +
            annotation_custom(addPlot, xmin=gr$xmin, xmax=gr$xmax,
                              ymin=gr$height+0.005, ymax=gr$heightup-0.005)
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
#' @param excludeFreq words with the frequency above
#'                    this threshold is excluded beforehand
#' @param excludeTfIdf filter using tfidf
#' @param orgDb organism database
#' @param tfidf use tfidf when making TDM
#' @param madeUpper words with uppercase
#' 
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
                            geneVecType="ENSEMBL", excludeTfIdf=NA,
                            numberOfWords=25, excludeFreq=5000,
                            orgDb=org.Hs.eg.db, tfidf=FALSE,
                            madeUpper=c("rna","dna")){
    
    ## Filter high frequency words if needed
    filterWords <- allFreqGeneSummary[
                            allFreqGeneSummary$freq > excludeFreq,]$word
    if (!is.na(excludeTfIdf)){
        filterWords <- c(filterWords,
            allTfIdfGeneSummary[
                            allTfIdfGeneSummary$tfidf > excludeTfIdf,]$word)
    }
    filterWords <- c(filterWords, "pmids", "geneid") ## Excluded by default

    # qqcat("filtered @{length(filterWords)} words (frequency | tfidf) ...\n")

    
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
                filter(label==NODES[1]) %>% select(x))
            XMAX <- as.numeric(labelPos %>%
                filter(label==NODES[length(NODES)]) %>% select(x))
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
                                    select(x))[,1])) >
                        mean(as.numeric((labelPos %>% filter(
                                    .data$label %in% names(lb[lb==2])) %>%
                                    select(x))[,1]))) {
                                                        R <- names(lb[lb==1])
                                                        L <- names(lb[lb==2])
                                                    } else {
                                                        L <- names(lb[lb==1])
                                                        R <- names(lb[lb==2])
                    }

                    if (length(names(geneVec)[geneVec %in% L])<geneNumLimit &
                        length(names(geneVec)[geneVec %in% R])<geneNumLimit)
                    {

                        pyrm <- returnPyramid(L,R, geneVec, geneVecType,
                                filterWords, numberOfWords, orgDb=orgDb,
                                additionalRemove=NA, madeUpper=madeUpper, tfidf=tfidf)
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
        k <- k + 1
    }
    return(grobList)
}


#' makeCorpus
#' 
#' Clean-up the corpus
#' 
#' @param docs corpus to clean
#' @param filterWords words to filter based on frequency
#' @param additionalRemove words to filter
#' 
#' @return cleaned corpus
#' @import tm
#' 
#' 
makeCorpus <- function (docs, filterWords, additionalRemove) {
    docs <- docs %>%
        tm_map(FUN=content_transformer(tolower)) %>%
        tm_map(FUN=removeNumbers) %>%
        tm_map(removeWords, stopwords::stopwords("english",
            "stopwords-iso")) %>%
        tm_map(removeWords, filterWords) %>% 
        tm_map(FUN=removePunctuation) %>%
        tm_map(FUN=stripWhitespace)
    if (prod(is.na(additionalRemove))!=1){
        docs <- docs %>% tm_map(removeWords, additionalRemove)
    }
    return(docs)
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
#' @param filterWords words to filter based on frequency
#' @param additionalRemove words to filter
#' @param numberOfWords the number of words to plot
#' @param orgDb database to change the ID to ENTREZ ID
#' @param widths parameters to pass to patchwork
#' @param lowCol gradient low color
#' @param highCol gradient high color
#' @param madeUpper words with uppercase
#' @param tfidf use tfidf
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
returnPyramid <- function(L, R, geneVec, geneVecType, filterWords,
                        numberOfWords, orgDb, additionalRemove=NA,
                        widths=c(0.3,0.3,0.3), lowCol="blue", tfidf=FALSE,
                        highCol="red", madeUpper=c("rna","dna")) {
    tb <- loadGeneSummary()
    ## Convert to ENTREZ ID
    geneList <- AnnotationDbi::select(orgDb,
        keys = names(geneVec)[geneVec %in% L],
        columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    filL <- tb %>% filter(Gene_ID %in% geneList)
    filL <- filL[!duplicated(filL$Gene_ID),]

    geneList <- AnnotationDbi::select(orgDb,
        keys = names(geneVec)[geneVec %in% R],
        columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    filR <- tb %>% filter(Gene_ID %in% geneList)
    filR <- filR[!duplicated(filR$Gene_ID),]
    
    all_L <- paste(filL$Gene_summary, collapse = " ")
    all_R <- paste(filR$Gene_summary, collapse = " ")
    all_bet <- c(all_L, all_R)
    
    all_bet <- VectorSource(all_bet)
    all_corpus <- VCorpus(all_bet)
    
    ## Clean the corpus
    all_corpus <- makeCorpus(all_corpus, filterWords, additionalRemove)
    if (tfidf){
        stop("Use of tfidf on returnPyramid is currently not supported ...")
        all_tdm <- TermDocumentMatrix(all_corpus, list(weighting = weightTfIdf))
    } else {
        all_tdm <- TermDocumentMatrix(all_corpus)
    }
    all_m <- as.matrix(all_tdm)
    
    ## Modified from: https://rpubs.com/williamsurles/316682
    common_words <- subset(all_m, all_m[, 1] > 0 & all_m[, 2] > 0)

    difference <- abs(common_words[, 1] - common_words[, 2])
    common_words <- cbind(common_words, difference)
    common_words <- common_words[order(common_words[, 3], decreasing = TRUE), ]
    
    L <- paste0(L, collapse="_")
    R <- paste0(R, collapse="_")
    
    # Referencing:
    # https://stackoverflow.com/questions/54191369/
    # how-to-create-pyramid-bar-chart-in-r-with-y-axis-labels-between-the-bars

    if (dim(common_words)[1]<numberOfWords){
        numberOfWords <- dim(common_words)[1]
    }

    if (numberOfWords==0){
        qqcat("No common words ...\n")
        return(NULL)
    } else {    
        topDf <- data.frame(
            ME = factor(c(rep(L,numberOfWords),
                rep(R,numberOfWords)),
                levels = c(L,R)),
            Word = c(common_words[1:numberOfWords, 1],
                common_words[1:numberOfWords, 2]),
            Label = c(rownames(common_words[1:numberOfWords, ]),
                rownames(common_words[1:numberOfWords, ]))
        )
        
        ## Convert to uppercase
        nodeName <- topDf$Label   
        for (i in madeUpper) {
            nodeName[nodeName == i] <- toupper(i)
        }
        topDf$Label <- nodeName
        
        gg1 <- topDf %>%
            mutate(Count = if_else(ME == L, Word, 0)) %>%
            ggplot(aes( Label, Count, fill = Count)) +
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
            filter(ME == L) %>%
            ggplot(aes(Label, 0, label = Label)) +
            geom_text(size=3.5) +
            coord_flip() +
            theme_void()
        
        gg3 <- topDf %>%
            mutate(Count = if_else(ME == R, Word, 0)) %>%
            ggplot(aes( Label, Count, fill = Count)) +
            geom_col(width = 0.6) +
            coord_flip() +
            theme_void() +
            scale_fill_gradient(low=lowCol,high=highCol)+
        theme(
                axis.title.y=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                legend.position = "none")
        

        pyramid <- gg1+gg2+gg3+plot_layout(widths=widths)
        pyramidGrob <- patchworkGrob(pyramid)
        return(pyramidGrob)
    }
}
