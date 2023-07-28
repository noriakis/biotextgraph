#' obtain_manual
#' 
#' obtain biotext class object from manually specified data
#' 
#' @param df data frame
#' @export
#' @examples obtain_manual(data.frame(text=c("test")))
#' @return biotext class object
#' 
obtain_manual <- function(df) {
    if (!is.data.frame(df)) {
      if (is.vector(df)) {
        df <- data.frame(df) |> `colnames<-`(c("text"))
      }
    }
    ret <- new("biotext")
    ret@type <- paste0("manual")
    ret@rawText <- df
    ret
}



#' obtain_alliance
#' 
#' obtain gene description from 
#' The Alliance of Genome Resources
#' https://www.alliancegenome.org/downloads
#' 
#' @param geneList gene ID list
#' @param file filename, default to "GENE-DESCRIPTION-TSV_HUMAN.tsv"
#' @param keyType key type of gene list
#' @param org_db used to convert ID
#' @export
#' @importFrom data.table fread
#' @examples \dontrun{obtain_alliance(c("PNKP"))}
#' @return biotext class object
obtain_alliance <- function(geneList, file="GENE-DESCRIPTION-TSV_HUMAN.tsv",
    keyType="SYMBOL", org_db=org.Hs.eg.db) {
    ret <- new("biotext")
    ret@query <- geneList
    ret@type <- "alliance_genome_resources"

    qqcat("Input genes: @{length(geneList)}\n")
    if (keyType!="SYMBOL"){
        geneList <- AnnotationDbi::select(org_db,
            keys = geneList, columns = c("SYMBOL"),
            keytype = keyType)$SYMBOL
        geneList <- geneList[!is.na(geneList)] |> unique()
        qqcat("  Converted input genes: @{length(geneList)}\n")
    }

    df <- fread(file, header=FALSE)
    df <- df |> dplyr::filter(.data$V2 %in% geneList)
    colnames(df) <- c("HGNC", "Gene_ID", "text")
    ret@rawText <- data.frame(df)
    ret
}


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
#' @examples obtain_bugsigdb("Veillonella dispar")
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
        abstDf <- do.call(pubmed, abstArg)
        ret@rawText <- abstDf
    }
    ret
}

#' obtain_enzyme
#' 
#' obtain EC-related text data from PubMed
#' 
#' @param file file downloaded from expasy
#' @param ec_num candidate ecnum, like those obtained from eggNOG-mapper
#' @param only_term only return quoted queries to pubmed
#' @param only_df only return ec description data.frame
#' if onlyTerm and onlyDf are both specified, onlyTerm have priority
#' @param tax_ec link taxonomy to EC using UniProt Taxonomy ID file
#' If this is TRUE, data.frame is returned
#' @param tax_file UniProt organism ID file path
#' @param cand_tax when taxec=TRUE, search only for these species.
#' @param arg_list passed to obtain_pubmed()
#' @param target abstract or title
#' @param api_key api key for PubMed
#' @export
#' @return biotext class object
obtain_enzyme <- function(file, ec_num,
    only_term=FALSE, only_df=FALSE, target="abstract",
    tax_ec=FALSE, tax_file=NULL, cand_tax=NULL, arg_list=list(),
    api_key=NULL) {
  flg <- FALSE
  candecs <- NULL
  allFlag <- FALSE
  if (length(ec_num)==1) {
    if (ec_num=="all") {
      allFlag <- TRUE
    }
  }
  qqcat("Processing EC file\n")
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (startsWith(line,"ID")) {
      ccs <- NULL
      ecs <- NULL
      drs <- NULL
      ec <- gsub("ID   ","",line)
      if (allFlag) {
        flg <- TRUE
      } else {
        if (ec %in% ec_num) {
          flg <- TRUE
        } else {
          flg <- FALSE
        }
      }
    }
    if (flg) {
      if (startsWith(line, "DE")) {
        de <- gsub("\\.","",gsub("DE   ","",line))
      }
      if (startsWith(line, "CC")) {
        cc <- gsub("\\.","",gsub("CC   ","",line))
        cc <- gsub("-!- ", "", cc)
        ccs <- c(ccs, cc)
      }
      if (startsWith(line, "DR")) {
        stLine <- gsub(" ","",gsub("DR", "", line))
        drs <- c(drs, unlist(strsplit(stLine,";")))
      }
      if (startsWith(line, "//")) {
        flg <- FALSE
        ecs <- c(ec, de, paste0(ccs, collapse=" "),
                 paste0(drs, collapse=";"))
        candecs <- rbind(candecs, ecs)
        }
    }
  }
  close(con)
  if (is.null(candecs)) {return(NULL)}
  candecs <- data.frame(candecs) |>
    `colnames<-`(c("number","desc","comment","DRs"))
  if (tax_ec) {
    qqcat("  Linking taxonomy to EC\n")
    retTaxEC <- NULL
    if (is.null(tax_file)) {stop("Please provide UniProt Taxonomy file")}
    if (!is.null(cand_tax)) {
      taxo <- getUPtax(tax_file, candUP="all", candTax=cand_tax)
      for (num in candecs$number) {
        desc <- subset(candecs, candecs$number==num)$desc
        allCharIDs <- as.character(unlist(vapply(unlist(strsplit(subset(candecs, candecs$number==num)$DRs,";")),
                                                 function(x) unlist(strsplit(x, "_"))[2], "character")))
        if (length(intersect(allCharIDs, taxo$UPID))>=1) {
          for (ta in intersect(allCharIDs, taxo$UPID)) {
            for (cta in subset(taxo, taxo$UPID==ta)$Taxonomy) {
              retTaxEC <- rbind(retTaxEC, c(num, desc, ta, cta))
            }
          }
        }
      }
    } else {
      stop("Please specify cand_tax when tax_ec=TRUE")
    }
    if (is.null(retTaxEC)) {stop("No EC could be found for query")}
    retTaxEC <- data.frame(retTaxEC) |> 
      `colnames<-`(c("number","desc","taxonomy","scientificName"))
    queryCheck <- NULL
    for (ct in cand_tax) {
      queryCheck <- cbind(queryCheck,
                          grepl(ct, retTaxEC$scientificName))
    }
    retTaxEC$query <- apply(queryCheck, 1, function(x) {
      if (length(cand_tax[x])==1) {
        cand_tax[x]
      } else {
        paste0(cand_tax[x],",")
      }})
    return(retTaxEC)
  }
  
  quoted <- dQuote(candecs$desc,options(useFancyQuotes = FALSE))
  if (only_term) {return(quoted)}
  if (only_df) {return(candecs)}
  arg_list[["target"]] <- "abstract"
  arg_list[["queries"]] <- quoted
  arg_list[["api_key"]] <- api_key
  abst <- do.call("obtain_pubmed", arg_list)
  abst@ec <- candecs
  abst@type <- "EC"
  return(abst)
}


#' obtain_enrich
#' 
#' obtain enrichment analysis description
#' 
#' @param geneList gene list
#' @param keyType key type of gene list
#' @param enrich `kegg` or `reactome`
#' @param org_db used to convert ID
#' @param top_path the number of pathways to be obtained
#' sorted by p-values
#' @export
#' @examples obtain_enrich(c("PNKP","DDX41"))
#' @return biotext class object
obtain_enrich <- function(geneList, keyType="SYMBOL", enrich="reactome",
    org_db=org.Hs.eg.db, top_path=30) {
    ret <- new("biotext")
    ret@query <- geneList
    ret@type <- "enrich"

    qqcat("Input genes: @{length(geneList)}\n")
    if (keyType!="ENTREZID"){
        geneList <- AnnotationDbi::select(org_db,
            keys = geneList, columns = c("ENTREZID"),
            keytype = keyType)$ENTREZID
        geneList <- geneList[!is.na(geneList)] |> unique()
        qqcat("  Converted input genes: @{length(geneList)}\n")
    }

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
    ret@rawText <- pathRes@result[1:top_path,"Description"] |> 
    data.frame() |> `colnames<-`(c("text"))
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
#' @examples obtain_refseq(c("PNKP"))
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
#' If filtered words are already present, add words to them
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
#' @examples obtain_refseq(c("PNKP")) |> set_filter_words()
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
    qqcat("Filtered @{length(unique(filterWords))} words (frequency and/or tfidf)\n")
    ret@filtered <- unique(filterWords)
    ret
}


#' perform_ora
#' @param ret biotext class
#' @param threshold ORA threshold to filter (bonferroni-corrected p-value)
#' @export
#' @examples obtain_refseq(c("DDX41")) |> perform_ora()
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
#' @examples obtain_refseq("PNKP") |> make_corpus()
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
#' @examples obtain_refseq("PNKP") |> make_corpus() |> make_TDM()
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
#' @param corMat tagging based on correlation matrix
#' based on ranking
#' @param mat correlation matrix
#' @return biotext class object
#' @examples obtain_refseq(c("IRF3","PNKP")) |> 
#' make_corpus() |> make_TDM() |> tag_words()
#' @export
tag_words <- function(ret, cl=FALSE, pvclAlpha=0.95, whole=FALSE,
	num_words=30, corMat=FALSE, mat=NULL) {
	
	if (corMat) {
		if (is.null(mat)) {stop("Please provide matrix")}
		pvc <- pvclust(as.matrix(as.dist(mat,diag=TRUE,upper=TRUE)))
	    pvcl <- pvpick(pvc, alpha=pvclAlpha)
	    ret@pvclust <- pvc
    	ret@pvpick <- pvcl
    	return(ret)
	}
	
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
#' @examples obtain_refseq(c("DDX41", "PNKP")) |> 
#' make_corpus() |> make_TDM() |> make_graph()
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
#' @param factorize convert to factor upon assigning
#' @return biotext class object
#' @examples refseq(c("PNKP","DDX41")) |> graph_cluster()
#' @export
graph_cluster <- function(ret, func=igraph::cluster_leiden, factorize=TRUE) {
	ret@communities <- do.call(func, list(graph=ret@igraphRaw))
    if (factorize) {
        V(ret@igraphRaw)$community <- factor(ret@communities$membership)
    } else {
        V(ret@igraphRaw)$community <- ret@communities$membership
    }
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
#' @param metab metabolite data.frame
#' @param metab_thresh threshold of association in metabolites - taxa relationship
#' (e.g., correlation coefficient)
#' @param metab_col metabolite data frame column name in the order of
#' "candidate taxon", "metabolite", "quantitative values for thresholding"
#' @param mb_plot microbe plotting
#' @param ec_file enzyme database file
#' @param up_tax_file UniProt taxonomy file
#' @export
#' @examples bugsigdb("Veillonella dispar") |> process_network_microbe()
#' @return biotext class object
#' 
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
        vtx <- data.frame(cbind(c(dis[,1], dis[,2]), c(rep("Diseases",length(dis[,1])),
            rep("Microbes",length(dis[,2]))))) |> `colnames<-`(c("name","type"))
        vtx <- vtx[!duplicated(vtx),]
        dmg <- tbl_graph(nodes=vtx, edges=data.frame(dis), directed=FALSE)
        addNet[["Diseases"]] <- dmg
    }

    if (ec_plot) {
        mb_plot <- TRUE
        if (is.null(ec_file)) {stop("Please provide EC file")}
        if (is.null(up_tax_file)) {stop("Please provide UniProt taxonomy file")}
        ecDf <- enzyme(file=ec_file, ecnum="all", taxec=TRUE,
            taxFile=up_tax_file, candTax=ret@query)
        if (!is.null(ecDf)) {
            ecDf <- ecDf[,c("desc","query")]
            vtx <- data.frame(
                cbind(c(ecDf[,1], ecDf[,2]),
                      c(rep("Enzymes",length(ecDf[,1])),
                        rep("Microbes",length(ecDf[,2]))))) |> 
            `colnames<-`(c("name","type"))
            vtx <- vtx[!duplicated(vtx),]
            ecGraph <- tbl_graph(nodes=vtx,
                edges=data.frame(ecDf),
                directed=FALSE)
            addNet[["Enzymes"]] <- ecGraph
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
            tmp <- metab[grepl(sp,metab[[metab_col[1]]]),]
            tmp <- tmp[abs(tmp[[metab_col[3]]])>metab_thresh,]
            if (dim(tmp)[1]!=0) {
                for (met in tmp[[metab_col[2]]]) {
                    metabGraph <- rbind(metabGraph, c(sp, met))
                }
            } else {
                qqcat("  Found no metabolites for @{sp}\n")
            }
        }
        if (!is.null(metabGraph)) {
            vtx <- data.frame(
                cbind(c(metabGraph[,1], metabGraph[,2]),
                      c(rep("Microbes",length(metabGraph[,1])),
                        rep("Metabolites",length(metabGraph[,2]))))) |> `colnames<-`(c("name","type"))
            vtx <- vtx[!duplicated(vtx),]
            metabGraph <- tbl_graph(nodes=vtx, edges=data.frame(metabGraph), directed=FALSE)
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
        
        incGene <- names(gcnt)[1:length(mb_list)]
        mbmap <- mbmap[mbmap[,2] %in% incGene,]
        ret@geneMap <- mbmap

        vtx <- data.frame(cbind(c(mbmap[,1], mbmap[,2]), c(rep("Words",length(mbmap[,1])),
            rep("Microbes",length(mbmap[,2]))))) |> `colnames<-`(c("name","type"))
        vtx <- vtx[!duplicated(vtx),]
        
        ####
        ## It is possible not to distinguish between query microbe name 
        ## and words as in `pubmed()`, as the words are lowercases here, 
        ## and altered to titlecase after
        ####

        mbmap <- tbl_graph(nodes=vtx,
            edges=data.frame(mbmap), directed=FALSE)
        V(coGraph)$type <- "Words"

        coGraph <- graph_join(as_tbl_graph(coGraph),
            as_tbl_graph(mbmap))
        coGraph <- coGraph |> activate("nodes") |>
            mutate(type=ifelse(is.na(.data$Freq),"Microbes","Words"))

        ## If present, add additional graphs
        if (length(addNet)!=0) {
            for (netName in names(addNet)) {
                tmpAdd <- addNet[[netName]]
                coGraph <- graph_join(as_tbl_graph(coGraph),
                    as_tbl_graph(tmpAdd))
                coGraph <- coGraph |> activate("nodes") |>
                    mutate(type=ifelse(is.na(.data$type),netName,.data$type))
            }
        }
        ## Set edge weight
        ## Probably set to NA would be better.


        E(coGraph)$edgeColor <- E(coGraph)$weight
        tmpW <- E(coGraph)$weight
        tmpW[is.na(tmpW)] <- ret@corThresh
        E(coGraph)$weight <- tmpW

    } else {
        V(coGraph)$type <- "Words"
        E(coGraph)$edgeColor <- E(coGraph)$weight
        coGraph <- as_tbl_graph(coGraph)
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
        V(coGraph)$tag_raw <- netCol

        ## Add disease and other labs
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

    nodeN <- (coGraph |> activate("nodes") |> data.frame())$type
    V(coGraph)$nodeCat <- nodeN

    if (!identical(ret@dic, logical(0))) {
        pdic <- ret@dic
        nodeDf <- coGraph |> activate("nodes") |> data.frame()
        V(coGraph)$name <- apply(nodeDf,
              1,
              function(x) {ifelse(x["type"]=="Words", pdic[x["name"]],
                x["name"])})
    }
    coGraph <- assign_community(ret, coGraph)  
      if (!is.tbl_graph(coGraph)) {
          ret@igraph <- coGraph
      } else {
          ret@igraph <- as.igraph(coGraph)
      }
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
#' @examples refseq(c("DDX41","PNKP")) |>
#' process_network_gene()
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
            if (ret@type=="refseq") {
                revID <- AnnotationDbi::select(org_db,
                    keys = as.character(ret@rawText$Gene_ID), 
                    columns = c("SYMBOL"),
                    keytype = "ENTREZID")$SYMBOL                
            } else {
                revID <- ret@rawText$Gene_ID
            }
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

        if (is.table(gcnt)) {
	        ret@geneCount <- gcnt
	    }
        
        incGene <- names(gcnt)[1:gene_plot_num]
        genemap <- genemap[genemap[,2] %in% incGene,]
        ret@geneMap <- genemap
        genemap <- simplify(igraph::graph_from_edgelist(genemap,
            directed = FALSE))
        coGraph <- tidygraph::graph_join(as_tbl_graph(coGraph),
            as_tbl_graph(genemap))
        coGraph <- coGraph |> activate("nodes") |>
            mutate(type=ifelse(is.na(.data$Freq),"Genes","Words"))
        E(coGraph)$edgeColor <- E(coGraph)$weight

        tmpW <- E(coGraph)$weight
        tmpW[is.na(tmpW)] <- ret@corThresh
        E(coGraph)$weight <- tmpW
    } else {
        E(coGraph)$edgeColor <- E(coGraph)$weight
        V(coGraph)$type <- "Words"
    }

    if (!is.null(gene_path_plot)) {
    	pathGraph <- return_gene_path_graph(ret, gene_path_plot, org_db,
    		threshold=gene_path_plot_threshold)
        pathGraph <- pathGraph |> data.frame()
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
    ## Assign node category
    nodeN <- (as_tbl_graph(coGraph) |> activate("nodes") |> data.frame())$type
    V(coGraph)$nodeCat <- nodeN
    names(nodeN) <- V(coGraph)$name

    ## Node attributes
    if (!identical(ret@dic, logical(0))) {
    	pdic <- ret@dic
        nodeDf <- as_tbl_graph(coGraph) |> activate("nodes") |> data.frame()
        V(coGraph)$name <- apply(nodeDf,
              1,
              function(x) {ifelse(x["type"]=="Words", pdic[x["name"]],
                x["name"])})
    }

    if (length(ret@pvpick)!=0) { ## If tag
        netCol <- tolower(names(V(coGraph)))
        for (i in seq_along(ret@pvpick$clusters)){
          for (j in ret@pvpick$clusters[[i]])
            netCol[netCol==j] <- paste0("cluster",i)
        }
        netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
        V(coGraph)$tag <- netCol
        V(coGraph)$tag_raw <- netCol

        
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
    coGraph <- assign_community(ret, coGraph)
      if (!is.tbl_graph(coGraph)) {
          ret@igraph <- coGraph
      } else {
          ret@igraph <- as.igraph(coGraph)
      }
    ret
}

#' return_gene_path_graph
#' @noRd
return_gene_path_graph <- function(ret, gene_path_plot="kegg",
	org_db=org.Hs.eg.db,
	threshold=0.05) {

    if (ret@type=="alliance_genome_resources") {
        gene_list <- AnnotationDbi::select(org_db,
            keys = ret@rawText$Gene_ID, columns = c("ENTREZID"),
            keytype = "SYMBOL")$ENTREZID
        gene_list <- gene_list[!is.na(gene_list)] |> unique()
    } else {
        gene_list <- ret@rawText$Gene_ID |> unique()
    }


    if (gene_path_plot == "reactome") {
        pathRes <- ReactomePA::enrichPathway(gene_list)
        pathRes@result$Description <- gsub("Homo sapiens\r: ",
                        "",
                        pathRes@result$Description)
    }
    else if (gene_path_plot == "kegg") {
        pathRes <- clusterProfiler::enrichKEGG(gene_list)
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


#' process_network_manual
#' 
#' process graph made from user-specified df
#' 
#' @param ret biotext class object
#' @param delete_zero_degree delete zero degree nodes
#' @param make_upper make these words uppercase
#' @param query_plot plot the column other than text
#' @param drop_ID drop `ID` column
#' @export
#' @return biotext class object
#' @examples 
#' obtain_manual(data.frame(text=c("test4 test4",
#'                                 "test3 test3 testi",
#'                                 "test2 test2"))) |>
#' make_corpus(num_only=TRUE) |>
#' make_TDM() |>
#' make_graph() |>
#' process_network_manual()
#' 

process_network_manual <- function(ret, delete_zero_degree=TRUE,
                       make_upper=NULL, query_plot=FALSE, drop_ID=TRUE) {
  DTM <- t(as.matrix(ret@TDM))
  
  if (query_plot) {
    if (!"query" %in% colnames(ret@rawText)) {stop("No query words specified")}
    row.names(DTM) <- ret@rawText$query
    
  }
  if (drop_ID) {
    ret@rawText$ID <- NULL
  }

  matSorted <- ret@wholeFreq
  coGraph <- as_tbl_graph(ret@igraphRaw)
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
  
  df <- ret@rawText
  incCols <- colnames(df)
  incCols <- incCols[!incCols %in% c("query","text")]
  V(coGraph)$type <- "Words"
  if (length(incCols)!=0) {
    qqcat("Including columns @{paste(incCols, collapse=' and ')} to link with query\n")
    for (ic in incCols) {
          querymap <- df[,c(ic, "query")]
          vtx <- data.frame(cbind(c(querymap[,2], querymap[,1]),
            c(rep("Queries",length(querymap[,2])),
              rep(ic,length(querymap[,1]))))) |> 
          `colnames<-`(c("name","type"))
          vtx <- vtx[!duplicated(vtx),]
          vtx <- vtx |> `rownames<-`(1:nrow(vtx))
          eds <- data.frame(querymap)
          words <- vtx |> subset(.data$type==ic)
          queriesDf <- vtx |> subset(.data$type=="Queries")
          row.names(words)[which(words$name %in% eds[,1])]
          row.names(queriesDf)[which(queriesDf$name %in% eds[,2])]
          eds[,1] <- sapply(eds[,1], function(x) {
            as.integer(row.names(words)[which(words$name %in% x)])
          })
          eds[,2] <- sapply(eds[,2], function(x) {
            as.integer(row.names(queriesDf)[which(queriesDf$name %in% x)])
          })
          qmap <- tbl_graph(nodes=vtx,edges=eds,directed=FALSE)
          coGraph <- graph_join(as_tbl_graph(coGraph),
                                qmap)
    }
  }
  
  nodeN <- NULL
  for (coln in c(incCols, "query")) {
    tmpn <- df[[coln]]
    tmpnn <- rep(coln, length(tmpn))
    names(tmpnn) <- tmpn
    nodeN <- c(nodeN, tmpnn)
  }

  if (query_plot) {
    genemap <- c()
    for (rn in nodeName){
      tmp <- DTM[ ,rn]
      for (nm in names(tmp[tmp!=0])){
        if (nm!=""){
          if (grepl(",",nm,fixed=TRUE)){
            for (nm2 in unlist(strsplit(nm, ","))){
              genemap <- rbind(genemap, c(rn, paste(nm2)))
            }
          } else {
            genemap <- rbind(genemap, c(rn, paste(nm)))
          }
        }
      }
    }
    vtx <- data.frame(cbind(c(genemap[,2], genemap[,1]),
      c(rep("Queries",length(genemap[,2])),
        rep("Words",length(genemap[,1]))))) |> 
    `colnames<-`(c("name","type"))
    vtx <- vtx[!duplicated(vtx),]
    vtx <- vtx |> `rownames<-`(1:nrow(vtx))
    eds <- data.frame(genemap)
    words <- vtx |> subset(.data$type=="Words")
    queriesDf <- vtx |> subset(.data$type=="Queries")
    row.names(words)[which(words$name %in% eds[,1])]
    row.names(queriesDf)[which(queriesDf$name %in% eds[,2])]
    eds[,1] <- sapply(eds[,1], function(x) {
      as.integer(row.names(words)[which(words$name %in% x)])
    })
    eds[,2] <- sapply(eds[,2], function(x) {
      as.integer(row.names(queriesDf)[which(queriesDf$name %in% x)])
    })
    qmap <- tbl_graph(nodes=vtx,edges=eds,directed=FALSE)
    coGraph <- graph_join(as_tbl_graph(coGraph),
                          qmap)
  }
  
  if (length(E(coGraph))==0) {stop("No edge present, stopping.")}
  E(coGraph)$edgeColor <- E(coGraph)$weight
  tmpW <- E(coGraph)$weight
  tmpW[is.na(tmpW)] <- ret@corThresh
  E(coGraph)$weight <- tmpW
  
  ## Node attributes
  if (length(ret@pvpick)!=0) { ## If tag
    netCol <- tolower(names(V(coGraph)))
    for (i in seq_along(ret@pvpick$clusters)){
      for (j in ret@pvpick$clusters[[i]])
        netCol[netCol==j] <- paste0("cluster",i)
    }
    netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
    V(coGraph)$tag <- netCol
    V(coGraph)$tag_raw <- netCol

    
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
  
  ## Assign node category, not tag
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
  
  coGraph <- assign_community(ret, coGraph)  

  if (!identical(ret@dic, logical(0))) {
    pdic <- ret@dic
    nodeDf <- as_tbl_graph(coGraph) |> activate("nodes") |> data.frame()
    V(coGraph)$name <- apply(nodeDf,
          1,
          function(x) {ifelse(x["type"]=="Words", pdic[x["name"]],
            x["name"])})
  }
  if (!is.tbl_graph(coGraph)) {
      ret@igraph <- coGraph
  } else {
      ret@igraph <- as.igraph(coGraph)
  }
  ret
}

#' assign_community
#' assign community and assign remaining queries the other category
#' If all words, just use like V(g)$community <- community$membership
#' @noRd
assign_community <- function(ret, coGraph) {
    if (!is.null(attributes(ret@communities)$names)) {
        memb <- ret@communities$membership
        uniq_memb <- memb |> unique()
        netCol <- tolower(names(V(coGraph)))
        for (i in seq_along(uniq_memb)){
            for (j in ret@communities$names[memb==uniq_memb[i]]) {
                netCol[netCol==j] <- paste0(uniq_memb[i])
            }
        }
        # netCol[!startsWith(netCol, "cluster")] <- "not_assigned"
        for (nn in seq_along(V(coGraph)$nodeCat)) {
            if (V(coGraph)$nodeCat[nn]!="Words") {
                netCol[nn] <- V(coGraph)$nodeCat[nn]
            }
        }
        V(coGraph)$community <- netCol
    }
    coGraph
}

#' plot_biotextgraph
#' 
#' plot biotextgraph after processing graph
#' 
#' @param ret biotext object
#' @param edge_link whether to use edge link or edge diagonal
#' @param edge_label whether to show edge label
#' @param show_legend whether to show legend
#' @param cat_colors named vector specifying color for each category of node
#' @param query_color color for queried nodes
#' @param font_family font family
#' @param colorize colorize the word frequency and query separately
#' if FALSE, no node is plotted for query nodes. it overrides color_by_tag.
#' @param color_by_tag color by tagging information (color scale will be discrete)
#' @param color_by_community color by community information in graph
#' @param pal palette for node color specified as (low, high)
#' @param tag_colors palette for node color in discrete mode
#' @param color_text whether to color the text
#' @param layout graph layout
#' @param na_edge_color edge color for NA value
#' @param use_seed seed value for ggrepel
#' @param scale_range scale range for node and text size
#' @param discrete_color_word colorize words by "Words" category, not frequency.
#' @param add_pseudo_freq add pseudo value for nodes other than words
#' @export
#' @examples refseq(c("PNKP","DDX41")) |> plot_biotextgraph()
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
    color_by_community=FALSE,
	tag_colors=NULL,
	color_text=TRUE,
	scale_range=c(5,10),
	pal=c("blue","red"),
	layout="nicely",
	na_edge_color="grey",
    add_pseudo_freq=FALSE,
    discrete_color_word=FALSE,
	use_seed=42) {

    if (color_by_community & color_by_tag) {
        stop("Cannot specify both of community or tag for coloring")
    }

	coGraph <- ret@igraph

    if (colorize | add_pseudo_freq) {
        fre <- V(coGraph)$Freq
        fre[is.na(fre)] <- min(fre, na.rm=TRUE)
        V(coGraph)$Freq <- fre
    }

    E(coGraph)$weightLabel <- round(E(coGraph)$weight, 3)

    if (color_by_community) {V(coGraph)$tag <- V(coGraph)$community}

    ## Make ggraph from here ##

    netPlot <- ggraph(coGraph, layout=layout)

    netPlot <- appendEdges(netPlot,
    	bn=is.directed(coGraph),
    	edgeLink=edge_link,
        edgeLabel=edge_label,
        showLegend=show_legend,
        fontFamily=font_family)

    cols <- V(coGraph)$nodeCat |> unique()
    if (is.null(cat_colors)) {
        cat_colors <- RColorBrewer::brewer.pal(length(cols), "Dark2")
        names(cat_colors) <- cols
        cat_colors["Genes"] <- query_color
        cat_colors["Microbes"] <- query_color
        cat_colors["query"] <- query_color
    }

    if (!is.null(V(coGraph)$tag) | color_by_community) {
        cols <- V(coGraph)$tag |> unique()
        if (is.null(tag_colors)) {
            tag_colors <- RColorBrewer::brewer.pal(length(cols), "Dark2")
            if (length(tag_colors) < length(cols)) {
                tag_colors <- colorRampPalette(RColorBrewer::brewer.pal(length(cols), "Dark2"))(length(cols))
            }
            names(tag_colors) <- cols
            tag_colors["Genes"] <- query_color
            tag_colors["query"] <- query_color
            tag_colors["Microbes"] <- query_color
        }
    }


    netPlot <- appendNodesAndTexts(netPlot,
        tag=color_by_tag | color_by_community,
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
        tagColors=tag_colors,
        discreteColorWord=discrete_color_word)

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
#' @examples refseq("DDX41", plotType="wc") |>
#' plot_wordcloud()
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