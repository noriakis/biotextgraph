
#' obtain_bugsigdb
#' 
#' obtain microbe-related information from bugsigdb
#' 
#' @param mb_list microbe list
#' @param target title or abstract
#' @param api_key API key for PubMed
#' @param curate curated articles (TRUE) or search for PubMed (FALSE)
#' if FALSE, same for fetching from PubMed by specified query
#' @export
#' @return biotext class object
obtain_bugsigdb <- function(mb_list,
	target="title",
	curate=TRUE,
	api_key=NULL) {

    ret <- new("biotext")
    ret@query <- mb_list
    ret@curate <- curate
    ret@type <- paste0("microbe_",target)

    qqcat("Input microbes: @{length(mb_list)}\n")
    tb <- importBugSigDB()
    subTb <- NULL
    for (m in mb_list) {
        tmp <- tb[grepl(m, tb$`MetaPhlAn taxon names`,
                 fixed = TRUE),]
        if (dim(tmp)[1]>0) {
            qqcat("  Found @{length(unique(tmp$PMID))} entries for @{m}\n")
            tmp$query <- m
            subTb <- rbind(subTb, tmp)
        }
    }
    qqcat("Including @{dim(subTb)[1]} entries\n")

    titles <- unique(subTb$Title)
    titles <- titles[!is.na(titles)]
    fil <- NULL
    for (title in titles){
        tmp <- subset(subTb, subTb$Title==title)
        if (dim(tmp)[1]>1){
            ## Duplicated entry is deleted based on query and paper title
            tmp$title2 <- paste0(tmp$Title,"_",tmp$query)
            tmp2 <- tmp[!duplicated(tmp$title2),]
            fil <- rbind(fil, c(paste(tmp2$query, collapse=","),
                unique(tmp2$Title), unique(tmp2$PMID)))
        } else {
            fil <- rbind(fil, c(tmp$query, tmp$Title, tmp$PMID))
        }
    }
    fil <- data.frame(fil)
    colnames(fil) <- c("query","text","ID")
    fil <- fil[!is.na(fil$ID),]
    ret@pmids <- fil$ID
    ret@rawTextBSDB <- subTb

	if (curate) {
        if (target=="abstract"){
            qqcat("Target is abstract\n")
            pmids <- unique(fil$ID)
            qqcat("  Querying PubMed for @{length(pmids)} pmids\n")
            queryUrl <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/",
                               "efetch.fcgi?db=pubmed&retmode=XML&id=",
                                   paste(pmids, collapse=","))
            if (!is.null(api_key)){
                queryUrl <- paste0(queryUrl, "&api_key=", api_key)
            } else {
                qqcat("  Querying without API key\n")
            }
            onequery <- url(queryUrl, open = "rb", encoding = "UTF8")
            xmlString <- readLines(onequery, encoding="UTF8")
            parsedXML <- xmlParse(xmlString)
            obtainText <- function(pmid) {xpathSApply(parsedXML,
                paste('//PubmedArticle/MedlineCitation[PMID=',pmid,
                    ']/Article/Abstract/AbstractText'),
                xmlValue)}
            PMIDs <- as.numeric(xpathSApply(parsedXML,
                "//PubmedArticle/MedlineCitation/PMID",
                xmlValue))
            abstDf <- NULL
            for (pmid in PMIDs) {
                if (length(obtainText(pmid))!=0) {
                    for (text in obtainText(pmid)) {
                        tax <- unique(subset(fil, fil$ID==pmid)$query)
                        abstDf <- rbind(abstDf, c(pmid, text, tax))
                    }
                }
            }
            abstDf <- data.frame(abstDf) |> `colnames<-`(c("ID","text","query"))
            ret@rawText <- abstDf
        } else {
        	ret@rawText <- fil
        }
    } else {
        abstArg[["queries"]] <- mb_list
        abstArg[["quote"]] <- TRUE
        abstArg[["target"]] <- target
        abstDf <- do.call(wcAbst, abstArg)
        ret@rawText <- abstDf
    }
    ret
}





#' obtain_refseq
#' 
#' obtain RefSeq text data
#' 
#' @param geneList gene list
#' @param keyType key type of gene list
#' @param organism organism
#' @param org_db used to convert ID
#' @export
#' @return biotext class object
obtain_refseq <- function(geneList, keyType="SYMBOL", organism=9606, org_db=org.Hs.eg.db) {
    ret <- new("biotext")
    ret@query <- geneList
    ret@type <- "refseq"

    qqcat("Input genes: @{length(geneList)}\n")
    if (keyType!="ENTREZID"){
        geneList <- AnnotationDbi::select(org_db,
            keys = geneList, columns = c("ENTREZID"),
            keytype = keyType)$ENTREZID
        geneList <- geneList[!is.na(geneList)] |> unique()
        qqcat("  Converted input genes: @{length(geneList)}\n")
    }

    tb <- loadGeneSummary(organism = organism)
    fil <- tb %>% filter(tb$Gene_ID %in% geneList)
    fil <- fil[!duplicated(fil$Gene_ID),]
    ret@rawText <- fil
    ret
}

#' obtain_pubmed
#' 
#' obtain PubMed text data
#' 
#' @param queries query list
#' @param delim delimiter for query, default to OR
#' @param api_key API key for PubMed
#' @param target title or abstract
#' @param ret_max max articles to get
#' @param sort_order sort by relevance or pubdate
#' @param limit query number limit
#' @param quote whether to quote each query
#' @export
#' @return biotext class object
obtain_pubmed <- function(queries, target="title",
	delim="OR", api_key=NULL, ret_max=10, limit=10,
	sort_order="relevance", quote=FALSE) {

    ret <- new("biotext")
    ret@type <- paste0("pubmed_",target)
    if (length(queries)>limit){
      stop("Number of queries exceeded specified limit number")}

    ret@delim <- delim
    if (quote) {
      query <- paste(dQuote(queries,
      	options(useFancyQuotes = FALSE)),
      collapse=paste0(" ",delim," "))
    } else {
      query <- paste(queries, collapse=paste0(" ",delim," "))
    }
    ret@query <- query
    clearQuery <- gsub('\"', '', queries)
    ret <- getPubMed(ret, query, clearQuery, type=target, apiKey=api_key,
                           retMax=ret_max, sortOrder=sort_order)
    ret@retMax <- ret_max
    ret@sortOrder <- sort_order
    ret
}

#' set_filter_words
#' 
#' set filtering words to class
#' 
#' @param ret biotext object
#' @param exclude GS for using GeneSummary, 
#' @param exclude_by `frequency` or `tfidf`
#' @param exclude_type ">" or "<"
#' @param exclude_number threshold to exclude
#' @param filterMax filter based on pre-computed maximum values 
#' for each word
#' @param additional_remove your customized words to remove
#' @param pre remove pre-defined words
#' @param pre_microbe remove pre-defined words for microbes
#' @export
#' @return biotext class object
set_filter_words <- function(ret, exclude_by="frequency",
	exclude_type=">", exclude="GS", exclude_number=2000, filterMax=FALSE,
	additional_remove=NULL, pre=TRUE, pre_microbe=FALSE) {

	filterWords <- ret@filtered

    if (exclude_by=="frequency") {
        pref = paste0(exclude,"_Freq")
    } else {
        pref = paste0(exclude,"_TfIdf")
    }
    if (filterMax) {
        pref <- paste0(pref, "_Max")
    }
    filterWords <- c(filterWords, retFiltWords(useFil=pref,
    	filType=exclude_type, filNum=exclude_number))
    if (!is.null(additional_remove)) {
    	filterWords <- c(filterWords, additional_remove)
    }
    if (pre) {
    	## [TODO] include easyPubMed PubMed_stopwords
	    filterWords <- c(filterWords,"genes","gene","patients","hub",
                      "analysis","cells","cell","expression","doi",
                      "deg","degs","author","authors","elsevier",
                      "oxford","wiley","pmids", "geneid", "pmid", "geneids")
	}
	if (pre_microbe) {
		filterWords <- c(filterWords, "microbiota",
        "microbiome","relative","abundance","abundances",
        "including","samples","sample","otu","otus","investigated",
        "taxa","taxon")
	}
    qqcat("Filtered @{length(filterWords)} words (frequency and/or tfidf)\n")
    ret@filtered <- filterWords
    ret
}


#' perform_ora
#' @param ret biotext class
#' @param threshold ORA threshold to filter (bonferroni-corrected p-value)
#' @export
#' @return biotext class object
perform_ora <- function(ret, threshold=0.05) {
	if (ret@type!="refseq") {stop("ORA for type other than refseq is not supported")}
    qqcat("Performing ORA\n")
    sig <- textORA(ret@rawText$Gene_ID |> unique())
    sigFilter <- names(sig)[p.adjust(sig, "bonferroni")>threshold]
    qqcat("Filtered @{length(sigFilter)} words (ORA)\n")
    ret@filtered <- c(ret@filtered, sigFilter)
    ret@ora <- sig
    ret
}


#' make_corpus
#' make corpus based on biotext class
#' @param ret biotext class object
#' @param collapse collapse all the sentences to one sentence
#' @param num_only delete number only (not deleting XXX123)
#' @param stem stem the words or not (default to FALSE)
#' @param preserve preserve the original cases
#' @param ngram n-gram
#' @export
#' @return biotext class object
make_corpus <- function(ret, collapse=FALSE,
	num_only=TRUE, stem=FALSE, preserve=TRUE, ngram=1) {
	if (ret@type=="refseq") {
		texts <- ret@rawText$Gene_summary
	# } else if (startsWith(a@type,"pubmed")) {
	} else {
		texts <- ret@rawText$text
	}
    if (collapse) {
        docs <- VCorpus(VectorSource(paste(texts, collapse=" ")))
    } else {
        docs <- VCorpus(VectorSource(texts))
    }
    ret@numOnly <- num_only
    ret@stem <- stem
    ret@ngram <- ngram
    if (preserve) {pdic <- preserveDict(docs, ngram,
    	num_only, stem)
            ret@dic <- pdic}

    docs <- makeCorpus(docs, ret@filtered, NULL, num_only, stem)
    ret@corpus <- docs
    ret@numOnly <- num_only
    ret@stem <- stem
    ret
}

#' make_TDM
#' 
#' @param ret biotext class
#' @param tfidf use TF-IDF
#' @param normalize normalize the values to document number
#' @param takeMean take the mean value for words to rank
#' @param takeMax take the max value for words to rank
#' @return biotext class object
#' @export
make_TDM <- function(ret, tfidf=FALSE,
	normalize=FALSE, takeMean=FALSE, takeMax=FALSE) {
	ngram <- ret@ngram
    docs <- ret@corpus
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
    ret@TDM <- docs
    ret
}

#' tag_words
#' 
#' tag the words based on pvclust and word matrix
#' 
#' @param ret biotext class object
#' @param cl cluster for parallel computing
#' @param pvclAlpha alpha value for `pvpick`
#' @param whole perform clustering on whole matrix (take time)
#' @param num_words if set, subset to this number of words
#' based on ranking
#' @return biotext class object
#' @export
tag_words <- function(ret, cl=FALSE, pvclAlpha=0.95, whole=FALSE,
	num_words=30) {
	DTM <- t(as.matrix(ret@TDM))
	freqWords <- names(ret@wholeFreq[1:num_words])
    if (whole){
        pvc <- pvclust(as.matrix(dist(t(DTM))), parallel=cl)
    } else {
    pvc <- pvclust(as.matrix(dist(
        t(
            DTM[, colnames(DTM) %in% freqWords]
            )
        )), parallel=cl)
    }
    pvcl <- pvpick(pvc, alpha=pvclAlpha)
    ret@pvclust <- pvc
    ret@pvpick <- pvcl
    ret
}

#' make_graph
#' 
#' make correlation or cooccurrence network
#' 
#' @param ret biotext class object
#' @param num_words number of words to include
#' @param cor_threshold correlation threshold
#' @param cooccurrence use cooccurrence
#' @param bn perform Bayesian network inference
#' @param R parameter for boot.strength function
#' @param whole return correlation or cooccurrence for whole words
#' @return biotext class object
#' @export
make_graph <- function(ret, num_words=30, cor_threshold=0.2,
	cooccurrence=FALSE, whole=FALSE, bn=FALSE, R=20) {
    matrixs <- obtainMatrix(ret, bn=bn, R=R, DTM=t(as.matrix(ret@TDM)),
    	freqWords=names(ret@wholeFreq[1:num_words]),
        corThresh=cor_threshold, cooccurrence=cooccurrence,
        onWholeDTM=whole)
    coGraph <- matrixs$coGraph
    ret <- matrixs$ret
    ret@numWords <- num_words
    ret@igraphRaw <- coGraph
    returnDf <- data.frame(word = names(ret@wholeFreq)[1:num_words],
    	freq=ret@wholeFreq[1:num_words])
    ret@freqDf <- returnDf[!is.na(returnDf$word),]
    ret
}

#' graph_cluster
#' @param ret biotext class object
#' @param func community detection algorithm in igraph
#' @return biotext class object
#' @export
graph_cluster <- function(ret, func=igraph::cluster_leiden) {
	ret@communities <- do.call(func, list(graph=ret@igraphRaw))
	ret
}


#' process_network_microbe
#' 
#' perform tasks before plotting
#' like obtaining associated diseases to plot
#' 
#' @param ret biotext class object
#' @param delete_zero_degree exclude zero degree nodes from plotting
#' @param make_upper make these words to uppercase
#' @param disease_plot disease plotting
#' @param ec_plot enzyme plotting
#' @param metab_plot metabolite plotting
#' @param mb_plot microbe plotting
#' @param ec_file enzyme database file
#' @param up_tax_file UniProt taxonomy file
process_network_microbe <- function(ret, delete_zero_degree=TRUE,
    make_upper=NULL, disease_plot=FALSE, ec_plot=FALSE,
    metab=NULL, metab_col=NULL, metab_thresh=0.2,
    mb_plot=FALSE, 
    ec_file=NULL, up_tax_file=NULL) {
    mb_list <- ret@query
    addNet <- list()
    subTb <- ret@rawTextBSDB                


    if (disease_plot) {## This does not need to be deduplicated
        mb_plot<-TRUE
        dis <- NULL
        for (i in seq_len(nrow(subTb))){
            dis <- rbind(dis,
            c(subTb[i, "Condition"], subTb[i, "query"]))
        }
        dis <- dis[!is.na(dis[,1]),]
        dis <- dis[!is.na(dis[,2]),]
        dmg <- simplify(graph_from_data_frame(dis, directed=FALSE))
        addNet[["Diseases"]] <- dmg
    }
    if (ec_plot) {
        mb_plot <- TRUE
        if (is.null(ec_file)) {stop("Please provide EC file")}
        if (is.null(up_tax_file)) {stop("Please provide UniProt taxonomy file")}
        ecDf <- wcEC(file=ec_file, ecnum="all", taxec=TRUE,
            taxFile=up_tax_file, candTax=mb_list)
        if (!is.null(ecDf)) {
           ecg <- simplify(graph_from_data_frame(ecDf[,c("desc","query")], 
                directed=FALSE))
           addNet[["Enzymes"]] <- ecg
        }
    }
    if (!is.null(metab)) {
        mb_plot<-TRUE
        if (is.null(metab_col)) {
            stop("No column names specified")
        }
        qqcat("Checking metabolites\n")
        metabGraph <- NULL
        for (sp in mb_list) {
            tmp <- metab[grepl(sp,metab[[metCol[1]]]),]
            tmp <- tmp[abs(tmp[[metCol[3]]])>metab_thresh,]
            if (dim(tmp)[1]!=0) {
                for (met in tmp[[metCol[2]]]) {
                    metabGraph <- rbind(metabGraph, c(sp, met))
                }
            } else {
                qqcat("  Found no metabolites for @{sp}\n")
            }
        }
        if (!is.null(metabGraph)) {
            metabGraph <- simplify(graph_from_edgelist(metabGraph,
                directed=FALSE))
            addNet[["Metabolites"]] <- metabGraph
        }
    }

    DTM <- t(as.matrix(ret@TDM))
    if (mb_plot) {
        row.names(DTM) <- ret@rawText$query
    }
    matSorted <- ret@wholeFreq

    coGraph <- ret@igraphRaw
    freqWords <- names(matSorted[1:ret@numWords])

    coGraph <- induced.subgraph(coGraph,
    names(V(coGraph)) %in% freqWords)
    V(coGraph)$Freq <- matSorted[V(coGraph)$name]

    if (delete_zero_degree){
        coGraph <- induced.subgraph(coGraph,
            degree(coGraph) > 0)
    }

    nodeName <- V(coGraph)$name
    dtmCol <- colnames(DTM)
    for (i in make_upper) {
        dtmCol[dtmCol == i] <- toupper(i)
        nodeName[nodeName == i] <- toupper(i)
    }
    V(coGraph)$name <- nodeName
    colnames(DTM) <- dtmCol

    nodeN <- NULL
    if (mb_plot) {
        mbmap <- c()
        for (rn in nodeName){
            tmp <- DTM[ ,rn]
            for (nm in names(tmp[tmp!=0])){
                if (grepl(",", nm, fixed = TRUE)) {
                    for (microbe in unlist(strsplit(nm, ","))){
                        mbmap <- rbind(mbmap, c(rn, microbe))
                    }
                } else {
                    mbmap <- rbind(mbmap, c(rn, nm))
                }
            }
        }
        mbmap <- mbmap[mbmap[,2]!="",]
        gcnt <- table(mbmap[,2])
        if (length(mb_list)!=1) {
            gcnt <- gcnt[order(gcnt, decreasing=TRUE)]
        }
        if (is.table(gcnt)) {
            ret@geneCount <- gcnt
        }
        
        inc_microbes <- names(gcnt)[1:length(mb_list)]
        mbmap <- mbmap[mbmap[,2] %in% inc_microbes,]
        ret@geneMap <- mbmap

        candMb <- unique(mb_list)
        nodeN <- rep("Microbes",length(candMb))
        names(nodeN) <- candMb

        mbmap <- simplify(graph_from_edgelist(mbmap,
            directed = FALSE))
        coGraph <- igraph::union(coGraph, mbmap)
        
        ## If present, add additional edges
        if (length(addNet)!=0) {
            for (netName in names(addNet)) {
                tmpAdd <- addNet[[netName]]
                tmpNN <- names(V(tmpAdd))
                tmpNN <- tmpNN[!tmpNN %in% names(nodeN)]

                newNN <- rep(netName, length(tmpNN))
                names(newNN) <- tmpNN
                nodeN <- c(nodeN, newNN)

                coGraph <- igraph::union(coGraph, tmpAdd)
            }
        }
        ## Set pseudo-edge weight
        ## Probably set to NA would be better.

        E(coGraph)$edgeColor <- E(coGraph)$weight
        tmpW <- E(coGraph)$weight
        tmpW[is.na(tmpW)] <- ret@corThresh
        E(coGraph)$weight <- tmpW

    } else {
        E(coGraph)$edgeColor <- E(coGraph)$weight
    }

    ## Node attributes
    if (length(ret@pvpick)!=0) { ## If tag
        netCol <- tolower(names(V(coGraph)))
        for (i in seq_along(ret@pvpick$clusters)){
            for (j in ret@pvpick$clusters[[i]])
                netCol[netCol==j] <- paste0("cluster",i)
        }
        netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
        V(coGraph)$tag <- netCol

        ## Add disease and other labs
        if (!is.null(nodeN)) {
            addC <- V(coGraph)$tag
            for (nn in seq_along(names(V(coGraph)))) {
                if (names(V(coGraph))[nn] %in% names(nodeN)) {
                    addC[nn] <- nodeN[names(V(coGraph))[nn]]
                } else {
                    next
                }
            }
            V(coGraph)$tag <- addC
        }
    }

    if (!is.null(nodeN)) {
        nodeCat <- NULL
        for (nn in seq_along(names(V(coGraph)))) {
            if (names(V(coGraph))[nn] %in% names(nodeN)) {
                nodeCat[nn] <- nodeN[names(V(coGraph))[nn]]
            } else {
                nodeCat[nn] <- "Words"
            }
        }
    } else {
        nodeCat <- rep("Words",length(V(coGraph)))
    }    
    V(coGraph)$nodeCat <- nodeCat

    if (!identical(ret@dic, logical(0))) {
        pdic <- ret@dic
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
    ret
}

#' process_network_gene
#' 
#' perform tasks before plotting
#' like obtaining associated genes to plot
#' 
#' @param ret biotext class object
#' @param delete_zero_degree exclude zero degree nodes from plotting
#' @param make_upper make these words to uppercase
#' @param make_upper_gene make gene symbol to uppercase
#' @param gene_plot plot and return associated genes
#' @param gene_plot_num how many genes are to be plotted
#' @param gene_path_plot plot associated pathways related to genes
#' @param gene_path_plot_threshold threshold for pathway enrichment
#' @param distinguish_query distinguish query with words in obtained text
#' @param org_db organism database to convert IDs
#' @return biotext class object
#' @export
process_network_gene <- function(ret, delete_zero_degree=TRUE,
	make_upper=NULL, make_upper_gene=TRUE, org_db=org.Hs.eg.db,
	gene_plot=FALSE, gene_plot_num=5, gene_path_plot=NULL,
	gene_path_plot_threshold=0.05, distinguish_query=FALSE) {


	DTM <- t(as.matrix(ret@TDM))
	matSorted <- ret@wholeFreq

    if (!is.null(gene_path_plot)) {gene_plot <- TRUE}
    if (gene_plot) {

		if (startsWith(ret@type,"pubmed")) {
			row.names(DTM) <- ret@rawText$query
		} else {
	        revID <- AnnotationDbi::select(org_db,
	            keys = as.character(ret@rawText$Gene_ID), 
	            columns = c("SYMBOL"),
	            keytype = "ENTREZID")$SYMBOL
	        row.names(DTM) <- revID
	    }
    }

	if (make_upper_gene) {
        make_upper <- c(make_upper,
        	tolower(keys(org_db, "SYMBOL")))
	}

	coGraph <- ret@igraphRaw
	freqWords <- names(matSorted[1:ret@numWords])

    coGraph <- induced.subgraph(coGraph,
    names(V(coGraph)) %in% freqWords)
    V(coGraph)$Freq <- matSorted[V(coGraph)$name]

    if (delete_zero_degree){
        coGraph <- induced.subgraph(coGraph,
            degree(coGraph) > 0)
    }

    nodeName <- V(coGraph)$name
    dtmCol <- colnames(DTM)
    for (i in make_upper) {
        dtmCol[dtmCol == i] <- toupper(i)
        nodeName[nodeName == i] <- toupper(i)
    }
    V(coGraph)$name <- nodeName
    colnames(DTM) <- dtmCol

    if (gene_plot) {
        genemap <- c()
        for (rn in nodeName){
	        tmp <- DTM[ ,rn]
	        for (nm in names(tmp[tmp!=0])){
	          if (nm!=""){
	            if (grepl(",",nm,fixed=TRUE)){
	              for (nm2 in unlist(strsplit(nm, ","))){
	                if (distinguish_query) {
	                  genemap <- rbind(genemap, c(rn, paste(nm2, "(Q)")))
	                } else {
	                  genemap <- rbind(genemap, c(rn, nm2))
	                }
	              }
	            } else {
	              if (distinguish_query) {
	                genemap <- rbind(genemap, c(rn, paste(nm, "(Q)")))
	              } else {
	                genemap <- rbind(genemap, c(rn, nm))
	              }
	            }
	          }
	        }
	    }
        gcnt <- table(genemap[,2])
        gcnt <- gcnt[order(gcnt, decreasing=TRUE)]

        if (!is.integer(gcnt)) {
	        ret@geneCount <- gcnt
	    }
        
        incGene <- names(gcnt)[1:gene_plot_num]
        genemap <- genemap[genemap[,2] %in% incGene,]
        ret@geneMap <- genemap
        genemap <- simplify(igraph::graph_from_edgelist(genemap,
            directed = FALSE))
        coGraph <- igraph::union(coGraph, genemap)
        E(coGraph)$edgeColor <- E(coGraph)$weight

        tmpW <- E(coGraph)$weight
        tmpW[is.na(tmpW)] <- ret@corThresh
        E(coGraph)$weight <- tmpW
    } else {
        E(coGraph)$edgeColor <- E(coGraph)$weight
    }

    if (!is.null(gene_path_plot)) {
    	pathGraph <- return_gene_path_graph(ret, gene_path_plot, org_db,
    		threshold=gene_path_plot_threshold)
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
    if (!identical(ret@dic, logical(0))) {
    	pdic <- ret@dic
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
    ret@igraph <- coGraph
    ret
}

#' return_gene_path_graph
#' @noRd
return_gene_path_graph <- function(ret, gene_path_plot="kegg",
	org_db=org.Hs.eg.db,
	threshold=0.05) {
    if (gene_path_plot == "reactome") {
        pathRes <- ReactomePA::enrichPathway(ret@rawText$Gene_ID)
        pathRes@result$Description <- gsub("Homo sapiens\r: ",
                        "",
                        pathRes@result$Description)
    }
    else if (gene_path_plot == "kegg") {
        pathRes <- clusterProfiler::enrichKEGG(ret@rawText$Gene_ID)
    }
    else {
        stop("Please specify 'reactome' or 'kegg'")
    }
    
    sigPath <- subset(pathRes@result, p.adjust<threshold)
    if (gene_path_plot=="kegg"){pathRes@keytype <- "ENTREZID"}
    ret@enrichResults <- pathRes@result
    sigPath <- subset(clusterProfiler::setReadable(pathRes,
        org_db)@result, p.adjust<threshold)

    if (dim(sigPath)[1]==0) {
        stop("No enriched term found.")
    } else {
        qqcat("Found @{dim(sigPath)[1]} enriched term\n")
    }

    pathGraph <- NULL
    for (i in 1:nrow(sigPath)){
        pa <- sigPath[i, "Description"]
        for (j in unlist(strsplit(sigPath[i, "geneID"], "/"))){
            pathGraph <- rbind(pathGraph, c(pa, j))
        }
    }
    pathGraph
}

#' plot_biotextgraph
#' 
#' plot biotextgraph after processing graph
#' 
#' @param ret biotext object
#' @param edge_link whether to use edge link or edge diagonal
#' @param edge_label whether to show edge label
#' @param show_legend whether to show legend
#' @param cat_color named vector specifying color for each category of node
#' @param query_color color for queried nodes
#' @param font_family font family
#' @param colorize colorize the word frequency and query separately
#' if FALSE, no node is plotted for query nodes. it overrides color_by_tag.
#' @param color_by_tag color by tagging information (color scale will be discrete)
#' @param pal palette for node color specified as (low, high)
#' @param tag_colors palette for node color in discrete mode
#' @param color_text whether to color the text
#' @param layout graph layout
#' @param na_edge_color edge color for NA value
#' @param use_seed seed value for ggrepel
#' @param scale_range scale range for node and text size
#' @export
#' @return biotext class object
plot_biotextgraph <- function(ret,
	edge_link=TRUE,
	edge_label=FALSE,
	show_legend=FALSE,
    cat_colors=NULL,
	query_color="grey",
	font_family="sans",
	colorize=FALSE,
	color_by_tag=FALSE,
	tag_colors=NULL,
	color_text=TRUE,
	scale_range=c(5,10),
	pal=c("blue","red"),
	layout="nicely",
	na_edge_color="grey",
	use_seed=42) {

	coGraph <- ret@igraph

    if (colorize) {
        fre <- V(coGraph)$Freq
        fre[is.na(fre)] <- min(fre, na.rm=TRUE)
        V(coGraph)$Freq <- fre
    }

    E(coGraph)$weightLabel <- round(E(coGraph)$weight, 3)

    if (color_by_tag) {
        netCol <- tolower(names(V(coGraph)))
        for (i in seq_along(ret@pvpick$clusters)){
            for (j in ret@pvpick$clusters[[i]])
                netCol[netCol==j] <- paste0("cluster",i)
        }
        netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
        for (nn in seq_along(V(ret@igraph)$nodeCat)) {
            if (V(ret@igraph)$nodeCat[nn]!="Words") {
                netCol[nn] <- V(ret@igraph)$nodeCat[nn]
            }
        }
        V(coGraph)$tag <- netCol
    }

    ## Make ggraph from here ##

    netPlot <- ggraph(coGraph, layout=layout)

    netPlot <- appendEdges(netPlot,
    	bn=is.directed(coGraph),
    	edgeLink=edge_link,
        edgeLabel=edge_label,
        showLegend=show_legend,
        fontFamily=font_family)

    cols <- V(ret@igraph)$nodeCat |> unique()
    if (is.null(cat_colors)) {
        cat_colors <- RColorBrewer::brewer.pal(length(cols), "Dark2")
        names(cat_colors) <- cols
        cat_colors["Genes"] <- query_color
        cat_colors["Microbes"] <- query_color
    }

    if (!is.null(V(coGraph)$tag)) {
        cols <- V(coGraph)$tag |> unique()
        if (is.null(tag_colors)) {
            tag_colors <- RColorBrewer::brewer.pal(length(cols), "Dark2")
            names(tag_colors) <- cols
            tag_colors["Genes"] <- query_color
            tag_colors["Microbes"] <- query_color
        }
    }

    netPlot <- appendNodesAndTexts(netPlot,
        tag=color_by_tag,
        colorize=colorize,
        nodePal=NULL,
        showLegend=show_legend,
        catColors=cat_colors,
        pal=pal,
        fontFamily=font_family,
        colorText=color_text,
        scaleRange=scale_range,
        useSeed=use_seed,
        ret=ret,
        tagColors=tag_colors)

    netPlot <- netPlot +
        scale_edge_width(range=c(1,3), name = "Correlation")+
        scale_edge_color_gradient(low=pal[1],high=pal[2],
            name = "Correlation", na.value=na_edge_color)+
        theme_graph()

    if (!is.null(V(coGraph)$grp)) {
        netPlot <- netPlot + ggforce::geom_mark_hull(
            aes(netPlot$data$x,
                netPlot$data$y,
                group = .data$grp,
                label=.data$grp, fill=.data$grp,
                filter = !is.na(.data$grp)),
            concavity = 4,
            expand = unit(2, "mm"),
            alpha = 0.25,
            na.rm = TRUE,
            show.legend = FALSE
        )
    }
    ret@net <- netPlot
    ret
}
#' plot_wordcloud
#' 
#' @param ret biotext object
#' @param num_words number of words to include
#' @param pal palette for colorizing words
#' @param make_upper make uppercase for these words
#' @param scale_freq scale the frequency for visualization
#' @param color_by_tag colorize based on tagging information
#' @param font_family font family
#' @param use_ggwordcloud use ggwordcloud or wordcloud function
#' @param arg_list arguments to pass to wordcloud functions
#' @param wc_scale scale factor for ggwordcloud
#' @export
#' @return biotext class object
plot_wordcloud <- function(ret, num_words=30, pal=palette(),
	make_upper=NULL, scale_freq=NULL, color_by_tag=FALSE,
	font_family="sans", use_ggwordcloud=TRUE, arg_list=list(),
	wc_scale=10) {
	matSorted <- ret@wholeFreq
	if (num_words > length(matSorted)) {num_words <- length(matSorted)}
    returnDf <- data.frame(word = names(matSorted)[1:num_words],freq=matSorted[1:num_words])

    if (length(ret@pvpick)!=0) {
        wcCol <- returnDf$word
        for (i in seq_along(ret@pvpick$clusters)){
            for (j in ret@pvpick$clusters[[i]])
                wcCol[wcCol==j] <- pal[i]
        }
        wcCol[!wcCol %in% pal] <- "grey"
    }

    for (i in make_upper) {
        returnDf[returnDf$word == i,"word"] <- toupper(i)
    }
    if (!identical(ret@dic, logical(0))) {
    	pdic <- ret@dic
        for (nm in unique(returnDf$word)) {
            if (nm %in% names(pdic)) {
                returnDf[returnDf$word == nm, "word"] <- pdic[nm]
            }
        }
    }

    if (!is.null(scale_freq)) {
        showFreq <- returnDf$freq*scale_freq
    } else {
        showFreq <- returnDf$freq
    }

    if (color_by_tag){
        arg_list[["words"]] <- returnDf$word
        arg_list[["freq"]] <- showFreq
        arg_list[["family"]] <- font_family
        arg_list[["colors"]] <- wcCol
        arg_list[["random.order"]] <- FALSE
        arg_list[["ordered.colors"]] <- TRUE
        if ("bg.color" %in% names(arg_list)) {
            arg_list[["bg.colour"]] <- arg_list[["bg.color"]]
        }
        if (use_ggwordcloud) {
            wc <- do.call(ggwordcloud::ggwordcloud, arg_list)+
            scale_size_area(max_size = wc_scale)+
            theme(plot.background = element_rect(fill="transparent",
                colour = NA))
        } else {
            wc <- as.ggplot(as_grob(~do.call("wordcloud", arg_list)))
        }
    } else {
        arg_list[["words"]] <- returnDf$word
        arg_list[["freq"]] <- showFreq
        arg_list[["family"]] <- font_family
        if ("bg.color" %in% names(arg_list)) {
            arg_list[["bg.colour"]] <- arg_list[["bg.color"]]
        }
        if (use_ggwordcloud) {
            wc <- do.call(ggwordcloud::ggwordcloud, arg_list)+
            scale_size_area(max_size = wc_scale)+
            theme(plot.background = element_rect(fill = "transparent",
                colour = NA))
        } else {
            wc <- as.ggplot(as_grob(~do.call("wordcloud", arg_list)))
        }
    }
    ret@freqDf <- returnDf[!is.na(returnDf$word),]
    ret@wc <- wc
    return(ret)
}