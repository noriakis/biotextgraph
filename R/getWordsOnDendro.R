#' getWordsOnDendro
#' 
#' Get grobs to plot on the dendrogram with the position information
#' 
#' @param dhc dendrogram
#' @param geneVec gene-named vector of node names in dendrogram
#' @param geneNumLimit when the gene number is above this threshold, pyramid plots are not produced 
#' @param geneVecType type of the name of geneVec (default: ENSEMBL)
#' @param numberOfWords the number of words to plot (default: 25)
#' @param excludeFreq words with the frequency above this threshold is excluded beforehand
#' @param orgDb organism database
#' 
#' @return list of pyramid plot grobs and its positions
#' @import tm
#' @importFrom ggdendro dendro_data
#' @importFrom dplyr select
#' @importFrom dendextend get_nodes_attr get_subdendrograms get_nodes_attr nnodes
#' 
#' @examples \dontrun{getWordsOnDendro(dhc, geneVec)}
#' @export
#' 
getWordsOnDendro <- function(dhc, geneVec, geneNumLimit=1000, geneVecType="ENSEMBL", numberOfWords=25, excludeFreq=10000, orgDb=org.Hs.eg.db){
    
    ## Filter high frequency words
    filterWords <- allFreqGeneSummary[allFreqGeneSummary$freq>excludeFreq,]$word
    filterWords <- c(filterWords, "pmids", "geneid") # 'PMIDs' is excluded by default
    
    grobList <- list()
    alHeights <- c(((dhc %>% get_subdendrograms(k=1))[[1]] %>% get_nodes_attr("height")))
    alNodes <- c(((dhc %>% get_subdendrograms(k=1))[[1]] %>% get_nodes_attr("label")))
    curHeights <- c(((dhc %>% get_subdendrograms(k=1))[[1]] %>% get_nodes_attr("height"))[1])
    k <- 2
    grobNum <- 1
    ddata <-dendro_data(dhc)
    labelPos <- ddata$labels
    
    while (k < length(dhc %>% labels)){
        subdendro <- dhc %>% get_subdendrograms(k=k)
        
        for (i in subdendro){
            NODES <- i %>% get_nodes_attr("label")
            NODES <- NODES[!is.na(NODES)]
            XMIN <- as.numeric(labelPos %>% filter(label==NODES[1]) %>% select(x))
            XMAX <- as.numeric(labelPos %>% filter(label==NODES[length(NODES)]) %>% select(x))
            HEIGHT <- get_nodes_attr(i, "height")[1]
            
            if (!HEIGHT %in% curHeights){
                if (nnodes(i)!=1){
                    lb <- i %>% cutree(k=2)
                    # L <- as.numeric(sapply(strsplit(names(lb[lb==1]), "ME"), "[", 2)) # WGCNA
                    # R <- as.numeric(sapply(strsplit(names(lb[lb==2]), "ME"), "[", 2)) # WGCNA
                    if (mean(as.numeric((labelPos %>% filter(label %in% names(lb[lb==1])) %>% select(x))[,1])) > mean(as.numeric((labelPos %>% filter(label %in% names(lb[lb==2])) %>% select(x))[,1]))) {
                        R <- names(lb[lb==1])
                        L <- names(lb[lb==2])
                    } else {
                        L <- names(lb[lb==1])
                        R <- names(lb[lb==2])
                    }

                    if (length(names(geneVec)[geneVec %in% L])<geneNumLimit & length(names(geneVec)[geneVec %in% R])<geneNumLimit)
                    {
                        grobList[[as.character(grobNum)]]$plot <- returnPyramid(L,R, geneVec, geneVecType, filterWords, numberOfWords, orgDb=orgDb, additionalRemove=NA)
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
        tm_map(removeWords, stopwords::stopwords("english", "stopwords-iso")) %>%
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
returnPyramid <- function(L, R, geneVec, geneVecType, filterWords, numberOfWords, orgDb, additionalRemove=NA,
                          widths=c(0.3,0.3,0.3), lowCol="blue", highCol="red") {
    tb <- loadGeneSummary()
    ## Convert to ENTREZ ID
    geneList = AnnotationDbi::select(orgDb, keys = names(geneVec)[geneVec %in% L], columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    filL <- tb %>% filter(Gene_ID %in% geneList)
    geneList = AnnotationDbi::select(orgDb, keys = names(geneVec)[geneVec %in% R], columns = c("ENTREZID"), keytype = geneVecType)$ENTREZID
    filR <- tb %>% filter(Gene_ID %in% geneList)
    
    all_L <- paste(filL$Gene_summary, collapse = "")
    all_R <- paste(filR$Gene_summary, collapse = "")
    all_bet <- c(all_L, all_R)
    
    all_bet <- VectorSource(all_bet)
    all_corpus <- VCorpus(all_bet)
    
    ## Clean the corpus
    all_corpus <- makeCorpus(all_corpus, filterWords, additionalRemove)
    all_tdm <- TermDocumentMatrix(all_corpus)
    all_m <- as.matrix(all_tdm)
    
    ## Modified from: https://rpubs.com/williamsurles/316682
    common_words <- subset(all_m, all_m[, 1] > 0 & all_m[, 2] > 0)
    difference <- abs(common_words[, 1] - common_words[, 2])
    common_words <- cbind(common_words, difference)
    common_words <- common_words[order(common_words[, 3], decreasing = T), ]
    
    L <- paste0(L, collapse="_")
    R <- paste0(R, collapse="_")
    
    # Referencing: https://stackoverflow.com/questions/54191369/how-to-create-pyramid-bar-chart-in-r-with-y-axis-labels-between-the-bars
    if (dim(common_words)[1]<numberOfWords){
        numberOfWords <- dim(common_words)[1]
    }
    
    topDf <- data.frame(
        ME = factor(c(rep(L,numberOfWords), rep(R,numberOfWords)), levels = c(L,R)),
        Word = c(common_words[1:numberOfWords, 1], common_words[1:numberOfWords, 2]),
        Label = c(rownames(common_words[1:numberOfWords, ]), rownames(common_words[1:numberOfWords, ]))
    )
    
    
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
