#' wcAbst
#' 
#' make word cloud or correlation network from PubMed abstracts using easyPubMed
#' 
#' 
#' @param queries gene symbols (max: 5)
#' @param redo if plot in other parameters, input the previous list
#' @param madeUpper make the words uppercase in resulting plot
#' @param palette palette for color gradient in correlation network
#' @param numWords the number of words to be shown
#' @param plotType "wc" or "network"
#' @param scaleRange scale for label and node size in correlation network
#' @param corThresh the correlation threshold
#' @param layout the layout for correlation network, defaul to "nicely"
#' @param edgeLink if FALSE, use geom_edge_diagonal
#' @param deleteZeroDeg delete zero degree node from plot in correlation network
#' @param showLegend whether to show legend in correlation network
#' @param colorText color text label based on frequency in correlation network
#' @param additionalRemove specific words to be excluded
#' @param ... parameters to pass to wordcloud()
#' 
#' @export
#' @examples wcAbst("DDX41")
#' @return list of data frame and ggplot2 object
#' @import tm easyPubMed
#' @import GeneSummary
#' @import wordcloud
#' @import igraph
#' @import ggraph ggplot2
#' @importFrom GetoptLong qqcat
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
wcAbst <- function(queries, redo=NA, madeUpper=c("dna","rna"),
                   palette=c("blue","red"), numWords=30, scaleRange=c(5,10),
                   showLegend=FALSE, plotType="wc", colorText=FALSE, corThresh=0.6,
                   layout="nicely", edgeLink=TRUE, deleteZeroDeg=TRUE, additionalRemove=NA, ...)
        {
        	if (!is.list(redo)) {
	        	fetched <- list() # store results

	        	if (length(queries)>5){stop("please limit the gene number to 5")}
				query <- paste(queries, collapse=" OR ")

				qqcat("querying pubmed for @{query}")

				pubmedIds <- get_pubmed_ids(query)
				pubmedData <- fetch_pubmed_data(pubmedIds)
				allXml <- articles_to_list(pubmedData)

				message("converting to a data frame ...\n")
				dataDf <- do.call(rbind, lapply(allXml, article_to_df, 
				                                max_chars = -1, getAuthors = FALSE))
				fetched[["rawdf"]] <- dataDf

				docs <- VCorpus(VectorSource(dataDf$abstract))
				docs <- makeCorpus(docs, additionalRemove, additionalRemove)

				docs <- TermDocumentMatrix(docs)
				mat <- as.matrix(docs)
				matSorted <- sort(rowSums(mat), decreasing=TRUE)
				fetched[["rawfrequency"]] <- matSorted
				fetched[["TDM"]] <- docs
			} else {
				fetched <- redo
				matSorted <- fetched[["rawfrequency"]]
				docs <- fetched[["TDM"]]
			}

		    if (plotType=="network"){
		        matSorted <- matSorted[1:numWords]
		        freqWords <- names(matSorted)
		        freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
		        
		        ## Check correlation
		        corData <- cor(freqWordsDTM)
		        fetched[["corMat"]] <- corData

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
		        netPlot <- netPlot + geom_node_point(aes(size=Freq, color=Freq), show.legend = showLegend)
		        if (colorText){
		            netPlot <- netPlot + geom_node_text(aes(label=name, size=Freq, color=Freq), check_overlap=TRUE, repel=TRUE,# size = labelSize,
		                           bg.color = "white", segment.color="black",
		                           bg.r = .15, show.legend=showLegend)
		        } else{
		            netPlot <- netPlot <- netPlot + geom_node_text(aes(label=name, size=Freq), check_overlap=TRUE, repel=TRUE,# size = labelSize,
		                           color = "black",
		                           bg.color = "white", segment.color="black",
		                           bg.r = .15, show.legend=showLegend) 
		        }
		        netPlot <- netPlot+
		            scale_size(range=scaleRange, name="Frequency")+
		            scale_edge_width(range=c(1,3), name = "Correlation")+
		            scale_color_gradient(low=palette[1],high=palette[2], name = "Frequency")+
		            scale_edge_color_gradient(low=palette[1],high=palette[2], name = "Correlation")+
		            theme_graph()
		        fetched[["net"]] <- netPlot
		    } else {
		    	matSorted <- matSorted[1:numWords]
		        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
		        for (i in madeUpper) {
		            # returnDf$word <- str_replace(returnDf$word, i, toupper(i))
		            returnDf[returnDf$word == i,"word"] <- toupper(i)
		        }
		        wc <- as.ggplot(as_grob(~wordcloud(words = returnDf$word, freq = returnDf$freq, ...)))
		        fetched[["df"]] <- returnDf
		        fetched[["wc"]] <- wc
		    }
    return(fetched)
}