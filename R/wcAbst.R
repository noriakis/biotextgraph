#' wcAbst
#' 
#' make word cloud or correlation network from PubMed abstracts using easyPubMed
#' 
#' 
#' @param queries gene symbols (max: 5)
#' @param redo if plot in other parameters, input the previous list
#' @param madeUpper make the words uppercase in resulting plot
#' @param pal palette for color gradient in correlation network
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
#' @param additionalRemove specific words to be excluded
#' @param target "abstract" or "title"
#' @param tag cluster the words based on text
#' @param genePlot plot associated genes (default: FALSE)
#' @param usefil filter based on "gstfidf" or "bsdbtfidf"
#' @param filnum specify filter tfidf
#' @param geneUpper make queries uppercase
#' @param apiKey api key for eutilities
#' @param tfidf use TfIdf when making TDM
#' @param pvclAlpha alpha for pvpick()
#' @param onlyCorpus return only corpus
#' @param onlyTDM return only TDM
#' @param preset filter preset words
#' @param ... parameters to pass to wordcloud()
#' 
#' @export
#' @examples wcAbst("DDX41")
#' @return list of data frame and ggplot2 object
#' @import tm easyPubMed
#' @import GeneSummary
#' @import wordcloud
#' @import igraph
#' @import stringr
#' @import ggraph ggplot2
#' @importFrom GetoptLong qqcat
#' @importFrom dplyr filter
#' @importFrom igraph graph.adjacency
#' @importFrom cowplot as_grob
#' @importFrom ggplotify as.ggplot
wcAbst <- function(queries, redo=NA, madeUpper=c("dna","rna"),
				   target="abstract", usefil=NA, filnum=0,
				   pvclAlpha=0.95,
				   geneUpper=FALSE, apiKey=NULL, tfidf=FALSE,
                   pal=c("blue","red"), numWords=30, scaleRange=c(5,10),
                   showLegend=FALSE, plotType="wc", colorText=FALSE,
                   corThresh=0.2, layout="nicely", tag=FALSE,
                   onlyCorpus=FALSE, onlyTDM=FALSE,
                   edgeLabel=FALSE, edgeLink=TRUE, ngram=NA, genePlot=FALSE,
                   deleteZeroDeg=TRUE, additionalRemove=NA,
                   preset=FALSE, ...)
        {
        	if (preset) {
        		additionalRemove <- c(additionalRemove,"genes","gene","patients","hub",
                                 "analysis","cells","cell","expression","doi",
                                 "deg","degs")
        	}
        	if (!is.list(redo)) {
	        	fetched <- list() # store results
	        	if (length(queries)>10){
	        		stop("please limit the gene number to 10")}
				# query <- paste(queries, collapse=" ")
	        	query <- paste(queries, collapse=" OR ")
				# qqcat("querying pubmed for @{query}\n")
				# allDataDf <- c()

				## Query all by OR and later grepl on title and abstract
				qqcat("querying pubmed for @{query}\n")
				pubmedIds <- get_pubmed_ids(query, api_key=apiKey)
				pubmedData <- fetch_pubmed_data(pubmedIds)
				allXml <- articles_to_list(pubmedData)
				qqcat("converting to a data frame ...\n")
				allDataDf <- do.call(rbind, lapply(allXml, article_to_df, 
				                        max_chars = -1, getAuthors = FALSE))
				incs <- c()
				for (i in queries){
				    li <- tolower(i)
				    inc <- grepl(li, tolower(allDataDf$abstract), fixed=TRUE)
				    inc[inc] <- i
				    inct <- grepl(li, tolower(allDataDf$title), fixed=TRUE)
				    inct[inct] <- i
				    incs <- cbind(incs, inc, inct)
				}
				allDataDf$query <- apply(incs, 1,
					function(x) paste(unique(x[x!="FALSE"]), collapse=","))

				# for (q in queries) {
				# 	qqcat("querying pubmed for @{q}\n")
				# 	pubmedIds <- get_pubmed_ids(q, api_key=apiKey)
				# 	pubmedData <- fetch_pubmed_data(pubmedIds)
				# 	allXml <- articles_to_list(pubmedData)

				# 	qqcat("converting to a data frame ...\n")
				# 	dataDf <- do.call(rbind, lapply(allXml, article_to_df, 
				# 	                        max_chars = -1, getAuthors = FALSE))
				# 	dataDf$query <- q
				# 	allDataDf <- rbind(allDataDf, dataDf)
				# }
				fetched[["rawdf"]] <- allDataDf
			} else {
				qqcat("Resuming from the previous results ...\n")
				fetched <- redo
				allDataDf <- fetched[["rawdf"]]
			}
			if (geneUpper){
				aq <- allDataDf$query
				aq <- aq[aq!=""]
				madeUpper <- c(madeUpper, tolower(unique(unlist(strsplit(aq, ",")))))
				# print(madeUpper)
			}
			if (target=="abstract"){
				docs <- VCorpus(VectorSource(allDataDf$abstract))
			} else if (target=="title"){
				docs <- VCorpus(VectorSource(allDataDf$title))
			} else {
				stop("specify target or abstract")
			}

			if (!is.na(usefil)){
				if (usefil=="gstfidf") {
					qqcat("filter based on GeneSummary\n")
					filterWords <- allTfIdfGeneSummary[
                                allTfIdfGeneSummary$tfidf > filnum,]$word
				} else if (usefil=="bsdbtfidf"){
					qqcat("filter based on BugSigDB\n")
					filterWords <- allTfIdfBSDB[
                                allTfIdfBSDB$tfidf > filnum,]$word
				} else {
					stop("please specify gstfidf or bsdbtfidf")
				}
			} else {
				filterWords <- NA
			}

			docs <- makeCorpus(docs, filterWords, additionalRemove)
			if (onlyCorpus){
				return(docs)
			}
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

		    if (onlyTDM){
		    	return(docs)
		    }

			mat <- as.matrix(docs)
			matSorted <- sort(rowSums(mat), decreasing=TRUE)
			fetched[["rawfrequency"]] <- matSorted
			fetched[["TDM"]] <- docs


		    if (plotType=="network"){
		        matSorted <- matSorted[1:numWords]
		        freqWords <- names(matSorted)
		        # freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
		        ## TODO: before or after?
		        freqWordsDTM <- t(as.matrix(docs))

                if (tag) {
                	if (is.list(redo) & "pvcl" %in% names(fetched)){
                		qqcat("Using previous pvclust results ...")
                		pvcl <- fetched[["pvcl"]]
                	} else {
			            pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
			            pvcl <- pvpick(pvc, alpha=pvclAlpha)
			            fetched[["pvcl"]] <- pvcl
			        }
		        }

		        ## Check correlation
		        corData <- cor(freqWordsDTM)
		        fetched[["corMat"]] <- corData

		        ## Set correlation below threshold to zero
		        corData[corData<corThresh] <- 0
		        coGraph <- graph.adjacency(corData, weighted=TRUE,
		        	mode="undirected", diag = FALSE)
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
		                	if (nm!=""){
		                		if (grepl(",",nm,fixed=TRUE)){
		                			for (nm2 in unlist(strsplit(nm, ","))){
		                				genemap <- rbind(genemap, c(rn, paste(nm2, "(Q)")))
		                			}
		                		} else {
				                    genemap <- rbind(genemap, c(rn, paste(nm, "(Q)")))
				                }
		                	}
		                }
		            }
		            genemap <- simplify(igraph::graph_from_edgelist(genemap, directed = FALSE))
		            coGraph <- igraph::union(coGraph, genemap)
		            tmpW <- E(coGraph)$weight
		            if (corThresh < 0.1) {corThreshGenePlot <- 0.01} else {
		                corThreshGenePlot <- corThresh - 0.1}
		            tmpW[is.na(tmpW)] <- corThreshGenePlot
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

		        ## Main plot
		        fetched[["ig"]] <- coGraph
		        netPlot <- ggraph(coGraph, layout=layout)

		        if (edgeLink){
		            if (edgeLabel){
		                netPlot <- netPlot +
		                            geom_edge_link(aes(width=weight,
		                            	color=weight,
		                                label=round(weight,3)),
		                                angle_calc = 'along',
		                                label_dodge = unit(2.5, 'mm'),
		                                alpha=0.5, show.legend = showLegend)
		            } else {
		                netPlot <- netPlot +
		                            geom_edge_link(aes(width=weight,
		                            	color=weight),
		                                alpha=0.5, show.legend = showLegend)
		            }
		        } else {
		            if (edgeLabel){
		                netPlot <- netPlot +
		                            geom_edge_diagonal(aes(width=weight,
		                            	color=weight,
		                                label=round(weight,3)),
		                                angle_calc = 'along',
		                                label_dodge = unit(2.5, 'mm'),
		                                alpha=0.5, show.legend = showLegend)
		            } else {
		                netPlot <- netPlot +
		                            geom_edge_diagonal(aes(width=weight,
		                            	color=weight),
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
		        netPlot <- netPlot+
		            scale_size(range=scaleRange, name="Frequency")+
		            scale_edge_width(range=c(1,3), name = "Correlation")+
		            scale_edge_color_gradient(low=pal[1],high=pal[2],
		            	name = "Correlation")+
		            theme_graph()
		        fetched[["net"]] <- netPlot
		    } else {
		    	matSorted <- matSorted[1:numWords]
		        returnDf <- data.frame(word = names(matSorted),freq=matSorted)
		        if (tag) {

		            freqWords <- names(matSorted)
		            freqWordsDTM <- t(as.matrix(docs[Terms(docs) %in% freqWords, ]))
		            pvc <- pvclust(as.matrix(dist(t(freqWordsDTM))))
		            pvcl <- pvpick(pvc)
		            fetched[["pvcl"]] <- pvcl

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
                fetched[["df"]] <- returnDf
        		fetched[["wc"]] <- wc
		    }
    return(fetched)
}