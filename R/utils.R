
#'
#' returnQuanteda
#' 
#' simple function to generate dfm and export to wc* functions.
#' It did not utilize full functionality of quanteda, and 
#' will add better support.
#' @param ret osplot object
#' @param quantedaArgs args to passed to tokens function
#' @param numWords number of words
#' @param ngram tokens_ngrams
#' @param tfidf use dfm_tfidf or not
#' @param filterWords filtered words
#' @param additionalRemove user-specified filter words
#' @return osplot object after filtering using quanteda
returnQuanteda <- function(ret, quantedaArgs,numWords,ngram,
                           filterWords,additionalRemove, tfidf) {
    
    if (length(quantedaArgList)==0) {
        quantedaArgList <- list(
            remove_punct=TRUE,
            remove_symbols=TRUE,
            remove_numbers=TRUE,
            remove_url=TRUE,
            remove_separators=TRUE,
            split_hyphens=FALSE
        )
    }

    docs <- quanteda::corpus(ret@rawText$text)
    ret@corpusQuanteda <- docs
    
    ## case_insensitive=TRUE
    
    quantedaArgList[["x"]] <- docs
    tokens <- do.call("tokens", quantedaArgList) |>
        quanteda::tokens_remove(quanteda::stopwords("english"))

    if (length(additionalRemove[!is.na(additionalRemove)])!=0) {
        tokens <- quanteda::tokens_remove(tokens, additionalRemove)
    }    
    if (length(filterWords[!is.na(filterWords)])!=0) {
        tokens <- quanteda::tokens_remove(tokens, filterWords)
    }
    if (!is.na(ngram)) {
        tokens <- quanteda::tokens_ngrams(tokens, n=ngram)
    }

    freqWordsDFM <- quanteda::dfm(tokens)
    if (tfidf) {
      freqWordsDFM <- quanteda::dfm_tfidf(freqWordsDFM)
    }
    ret@dfm <- freqWordsDFM
    # freqWordsAll <- sort(quanteda::featfreq(freqWordsDFM),
    #                   decreasing=TRUE)
    # freqWords <- names(freqWordsAll[1:numWords])
    return(ret)
}


#' convertMetaCyc
#'
#' this function needs taxonomizr package to be installed,
#' and needs download of databases by prepareDatabase(getAccessions=FALSE)
#' 
#' @param ids tax ids from metacyc
#' @param onlySpecies parse only species
#' @return coverted species name
#' @examples
#' convertMetaCyc("TAX-9606")
#' @export
convertMetaCyc <- function (ids, onlySpecies=FALSE) {
  convIds <- sapply(ids, function(x) if (grepl("TAX-",x)) unlist(strsplit(x,"-")[[1]])[2] else x)
  if (onlySpecies) {
    parsed <- apply(taxonomizr::getTaxonomy(convIds), 1, function(x) x[7])
  } else {
    parsed <- apply(taxonomizr::getTaxonomy(convIds), 1, function(x) paste(x, collapse=";"))
  }
  as.character(parsed)
}


#' clearPath
#' 
#' clear MetaCyc pathway strings
#' 
#' @param ex vector of text
#' @noRd
clearPath <- function (ex) {
    ex <- gsub("|FRAME: ", "", ex)
    ex <- gsub("FRAME: ", "", ex)
    ex <- gsub("|CITS: ", "", ex)
    ex <- gsub("CITS: ", "", ex)
    ex <- gsub("\\[[^][]*]", "", ex)
    ex <- gsub("\"", "", ex)
    ## Clean HTML tags
    ex <- gsub("<.*?>", "", ex)
    ## lastly remove \\|
    ex <- gsub("\\|","",ex)
    ex
}


#'
#' parseMetaCycPathway
#' 
#' parse MetaCyc "pathways.dat"
#' 
#' 
#' @param file path to pathways.dat
#' @param candSp species to grepl
#' @param withTax parse taxonomy information
#' @param noComma no comma separated when taxonomy parsing
#' @param clear delete HTML tags and some symbols
#' @return data.frame of MetaCyc pathway information
#' @examples
#' file <- "pathways.dat"
#' \donttest{parseMetaCycPathway(file, candSp="all")}
#' @export
#' 
parseMetaCycPathway <- function(file, candSp, withTax=FALSE, noComma=FALSE, clear=FALSE) {
  flg <- FALSE
  allFlg <- FALSE
  allmeta <- NULL
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (startsWith(line,"UNIQUE-ID - ")) {
      com <- NA
      commn <- NA
      spec <- NULL
      taxr <- NULL
      pwy <- gsub("UNIQUE-ID - ","",line)
      flg <- TRUE
    }
    if (flg) {
      if (startsWith(line, "/")) {
        com <- c(com, line)
      }
      if (startsWith(line, "COMMON-NAME")) {
        commn <- gsub("COMMON-NAME - ","",line)
      }
      if (startsWith(line, "SPECIES - ")) {
        spec <- c(spec, gsub("SPECIES - ","",line))
      }
      if (startsWith(line, "TAXONOMIC-RANGE - ")) {
        taxr <- c(taxr, gsub("TAXONOMIC-RANGE - ","",line))
      }
      if (startsWith(line,"//")) {
        coms <- paste(com[!is.na(com)], collapse=" ")
        coms <- gsub("/","",coms)

        if (!noComma) {
          if (length(spec)!=0) {spec <- paste0(spec, collapse=",")} else {spec <- ""}
          if (length(taxr)!=0) {taxr <- paste0(taxr, collapse=",")} else {taxr <- ""}          
        } else {
          if (length(spec)!=0) {} else {spec <- ""}
          if (length(taxr)!=0) {} else {taxr <- ""}
        }
        if (length(candSp)==1) {
          if (candSp=="all") {
             allFlg <- TRUE
          }
        }

        if (!allFlg) {
          if (grepl(paste(candSp,collapse="|"),coms)) {
            if (withTax) {
              if (noComma) {
                for (sp in spec) {
                  for (tax in taxr) {
                    allmeta <- rbind(allmeta, c(pwy, coms, commn, sp, tax))
                  }
                }
              } else {
                allmeta <- rbind(allmeta, c(pwy, coms, commn, spec, taxr))
              }
            } else {
              allmeta <- rbind(allmeta, c(pwy, coms, commn))
            }
          }
        } else {
            if (withTax) {
              if (noComma) {
                for (sp in spec) {
                  for (tax in taxr) {
                    allmeta <- rbind(allmeta, c(pwy, coms, commn, sp, tax))
                  }
                }
              } else {
                allmeta <- rbind(allmeta, c(pwy, coms, commn, spec, taxr))
              }
            } else {
              allmeta <- rbind(allmeta, c(pwy, coms, commn))
            }
        }
        flg <- FALSE
      }
    }
  }
  close(con)
  if (withTax) {
    allmeta <- data.frame(allmeta) |> `colnames<-`(c("pathwayID","text","commonName","species","taxonomicRange"))
  } else {
    allmeta <- data.frame(allmeta) |> `colnames<-`(c("pathwayID","text","commonName"))
  }
  # allmeta |> dim()
  if (!allFlg) {
    queries <- NULL
    for (q in candSp) {
      queries <- cbind(queries, grepl(q, allmeta$text))
    }
    allmeta$query <- apply(queries, 1, function(x) paste0(candSp[x], collapse=","))
  }

  if (clear) {
    allmeta$text <- clearPath(allmeta$text)
  }
  return(allmeta)
}

#' retFiltWords
#' return filtered words
#' @param useFil use filter
#' @param filType "above" or "below"
#' @param filNum number
#' @noRd

retFiltWords <- function(useFil, filType, filNum) {
    if (filType=="above" | filType==">") {
        if (useFil=="GS_TfIdf") {
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummary[
            allTfIdfGeneSummary$tfidf > filNum,]$word
        } else if (useFil=="BSDB_TfIdf"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDB[
            allTfIdfBSDB$tfidf > filNum,]$word
        } else if (useFil=="GS_TfIdf_Max"){
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummaryMax[
            allTfIdfGeneSummaryMax$tfidf > filNum,]$word
        } else if (useFil=="BSDB_TfIdf_Max"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDBMax[
            allTfIdfBSDBMax$tfidf > filNum,]$word
        } else if (useFil=="GS_Freq"){
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allFreqGeneSummary[
            allFreqGeneSummary$freq > filNum,]$word
        } else if (useFil=="BSDB_Freq"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allFreqBSDB[
            allFreqBSDB$freq > filNum,]$word
        } else {
          stop("Please specify useFil")
        }
    } else {
        if (useFil=="GS_TfIdf") {
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummary[
            allTfIdfGeneSummary$tfidf < filNum,]$word
        } else if (useFil=="BSDB_TfIdf"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDB[
            allTfIdfBSDB$tfidf < filNum,]$word
        } else if (useFil=="GS_TfIdf_Max"){
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allTfIdfGeneSummaryMax[
            allTfIdfGeneSummaryMax$tfidf < filNum,]$word
        } else if (useFil=="BSDB_TfIdf_Max"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allTfIdfBSDBMax[
            allTfIdfBSDBMax$tfidf < filNum,]$word
        } else if (useFil=="GS_Freq"){
          qqcat("Filter based on GeneSummary\n")
          filterWords <- allFreqGeneSummary[
            allFreqGeneSummary$freq < filNum,]$word
        } else if (useFil=="BSDB_Freq"){
          qqcat("Filter based on BugSigDB\n")
          filterWords <- allFreqBSDB[
            allFreqBSDB$freq < filNum,]$word
        } else {
          stop("Please specify useFil")
        }
    }
    return(filterWords)
}


#' removeAloneNumbers
#' 
#' @noRd
#' 
removeAloneNumbers <- 
    function (x) PlainTextDocument(
        gsub('\\s*(?<!\\B|-)\\d+(?!\\B|-)\\s*', " ",
                      x, perl=TRUE))

#' preserveDict
#' 
#' 
#' @param docs corpus
#' @param ngram ngram params
#' @param numOnly remove number only
#' @param stem use stemming
#' 
#' @noRd
#' 
preserveDict <- function(docs, ngram, numOnly, stem) {
    NgramTokenizer <- function(x)
    unlist(lapply(ngrams(words(x), ngram),
        paste, collapse = " "),
        use.names = FALSE)

  ldocs <- docs %>%
      tm_map(FUN=content_transformer(tolower))
  if (numOnly) {
    ldocs <- ldocs %>% tm_map(FUN=removeAloneNumbers)
    docs <- docs %>% tm_map(FUN=removeAloneNumbers)
  } else {
    ldocs <- ldocs %>% tm_map(FUN=removeNumbers)
    docs <- docs %>% tm_map(FUN=removeNumbers)
  }
  docs <- docs %>%
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
  ldocs <- ldocs %>%
    tm_map(FUN=removePunctuation) %>%
    tm_map(FUN=stripWhitespace)
  if (stem) {
    docs <- docs %>% tm_map(stemDocument) 
    ldocs <- ldocs %>% tm_map(stemDocument)
  }
  if (!is.na(ngram)){
      docs <- TermDocumentMatrix(docs,
                                 control = list(tokenize = NgramTokenizer,
                                                tolower=FALSE))
      ldocs <- TermDocumentMatrix(ldocs,
                                    control = list(tokenize = NgramTokenizer,
                                                   tolower=TRUE))
  } else {
      docs <- TermDocumentMatrix(docs, control=list(tolower=FALSE))
      ldocs <- TermDocumentMatrix(ldocs, control=list(tolower=TRUE))
  }
  rawrn <- row.names(docs)
  lrn <- row.names(ldocs)
  dic <- list()
  for (i in rawrn) {
    if (tolower(i) %in% lrn) {
      dic[[tolower(i)]] <- i      
    }
  }
  rawWords <- names(dic)
  newWords <- unlist(dic)
  names(newWords) <- rawWords
  return(newWords)
}

#' getPubMed
#' 
#' obtain pubmed information
#' 
#' @param ret osplot object
#' @param searchQuery search query
#' @param rawQuery raw query
#' @param type abstract or title
#' @param apiKey api key
#' @param retMax retmax
#' @import rentrez
#' 
#' @noRd

getPubMed <- function(ret, searchQuery, rawQuery,
    type="abstract", apiKey=NULL, retMax=10, sortOrder="relevance") {
    if (is.null(apiKey)){
      qqcat("Proceeding without API key\n")
    } else {
      set_entrez_key(apiKey)
    }
    pubmedSearch <- entrez_search("pubmed",
                                  term = searchQuery, 
                                  retmax = retMax,
                                  sort = sortOrder)
    ret@pmids <- pubmedSearch$ids
    searchResults <- entrez_fetch(db="pubmed",
                                  pubmedSearch$ids, rettype="xml", 
                                  parsed=FALSE)
    parsedXML <- xmlTreeParse(as.character(searchResults))

    if (type=="abstract") {
      obt <- "Abstract"
    } else {
      obt <- "ArticleTitle"
    }

    charset <- xmlElementsByTagName(parsedXML$doc$children$PubmedArticleSet,
                                    obt,
                                    recursive = TRUE)

    obtainedText <- as.character(xmlValue(charset))

    incs <- c()
    for (i in rawQuery){
      li <- tolower(i)
      inc <- grepl(li, tolower(obtainedText), fixed=TRUE)
      inc[inc] <- i
      # inct <- grepl(li, tolower(titletext), fixed=TRUE)
      # inct[inct] <- i
      incs <- cbind(incs, inc)
    }
    obtainedDf <- data.frame(cbind(
      seq_len(length(obtainedText)),
      obtainedText,
      apply(incs, 1, function(x) paste(unique(x[x!="FALSE"]), collapse=","))
    )) |> `colnames<-`(c("ID","text","query"))
    ret@rawText <- obtainedDf
    return(ret)
}


#' findTerm
#' 
#' find queried terms in list of gene clusters and return frequency
#' 
#' @param query query words
#' @param listOfGenes named list of genes
#' @param split split the query by space, default to FALSE
#' @param ngram use ngram or not
#' @param tfidf use tfidf when creating TDM or not
#' @param keyType key type of listOfGenes
#' @param calc "sum", "mean" or "highest"
#' @param argList passed to wcGeneSummary
#' @return named list of frequency
#' @export
#' @examples
#' query <- "DNA repair"
#' lg <- list()
#' lg[["sample"]] <- c("ERCC1","ERCC2")
#' findTerm(query, lg)
#' 
findTerm <- function (query, listOfGenes, split=FALSE, ngram=NA,
                      tfidf=TRUE, calc="sum", keyType="SYMBOL", argList=list()) {
    qqcat("Finding query in @{length(listOfGenes)} clusters ...\n")
    if (split) {
        querySplit <- tolower(unlist(strsplit(query, " ")))
    } else {
        querySplit <- tolower(query)
    }
    frq <- list()
    for (clus in names(listOfGenes)) {
        argList[["geneList"]] <- listOfGenes[[clus]]
        argList[["keyType"]] <- keyType
        argList[["ngram"]] <- ngram
        argList[["tfidf"]] <- tfidf
        argList[["onlyTDM"]] <- TRUE
        tmptdm <- do.call("wcGeneSummary", argList)
        # tmptdm <- wcGeneSummary(listOfGenes[[clus]],
        #                         keyType = keyType, ngram=ngram,
        #                         tfidf=tfidf, onlyTDM=TRUE)
        querytdm <- t(as.matrix(tmptdm[Terms(tmptdm) %in% querySplit, ]))
        tmp <- rep(0, length(querySplit))
        names(tmp) <- querySplit
        if (calc=="sum"){
            tmpfrq <- apply(querytdm, 2, sum)# / dim(querytdm)[1]
        } else if (calc=="mean") {
            tmpfrq <- apply(querytdm, 2, mean)# / dim(querytdm)[1]
        } else if (calc=="highest") {
            tmpfrq <- apply(querytdm, 2, max)
        } else {
            stop("please specify sum, mean or highest")
        }
        for (i in names(tmpfrq)){
            tmp[names(tmp)==i] <- tmpfrq[i]
        }
        frq[[clus]] <- tmp
    }
    frq
}

#' returnSim
#' 
#' return similarity matrix of cluster
#' based on intersection over union of words
#' 
#' @param cllist cluster list (named vector or list)
#' @param keyType keytype
#' @param numLimit threshold for gene number limit
#'                 default to 5000
#' @param target target to query
#' @param argList parameters to pass to wcGeneSummary()
#' @return similarity matrix
#' @export
#' @examples
#' ex <- returnExample()
#' returnSim(ex$color, keyType="ENSEMBL")
#' @importFrom GetoptLong qqcat
returnSim <- function (cllist, keyType="ENTREZID", numLimit=5000,
  target="refseq", argList=list()) {
    store <- list()
    if (!is.list(cllist)){
        converted <- list()
        clname <- unique(cllist)
        for (c in clname) {
            converted[[as.character(c)]] <-
                names(cllist)[cllist==c]
        }
    } else {
        converted <- cllist
    }
    qqcat("Number of clusters: @{length(converted)}\n")
    for (i in names(converted)) {
        ## Time consuming cluster
        if (length(converted[[i]])<numLimit){
            qqcat("@{i}\n")
            argList[["geneList"]] <- converted[[i]]
            argList[["keyType"]] <- keyType
            store[[i]] <-
                do.call("wcGeneSummary", argList)@freqDf
        }
    }
    sim <- sapply(store, function(x) sapply(store,
                    function(y)
                        length(intersect(x$word, y$word)) /
                        length(union(x$word, y$word))))
    return(as.matrix(sim))
}

#' makeBar
#' 
#' Makeing a barplot of word frequency from queried genes
#' 
#' @examples
#' geneList <- c("DDX41")
#' makeBar(geneList)
#' @param queries gene IDs
#' @param keyType default to SYMBOL
#' @param top how many numbers of words to be shown
#' @param grad use gradient of frequency
#' @param pal palette used in barplot
#' @param textSize text size in barplot
#' @param reord order by frequency or not
#' @param flip flip the barplot (gene name in y-axis)
#' @param orgDb orgDb
#' @param retList return result of wcGeneSummary
#' @param argList passed to wcGeneSummary()
#' @import org.Hs.eg.db
#' @return barplot of word frequency
#' @importFrom stats reorder
#' @export
#' 
makeBar <- function(queries, top=10, keyType="SYMBOL",
                    pal=NULL, textSize=20, reord=TRUE, orgDb=org.Hs.eg.db,
                    flip=FALSE, grad=FALSE, retList=FALSE, argList=list()) {
  if (is.null(pal)) {
    # palNum <- sample(1:151,1)
    # pal <- pokepal(palNum)
    pal <- palette()
    if (length(pal)<top){
      pal <- rep(pal, ceiling(top/length(pal)))
    }
  }

  argList[["geneList"]] <- queries
  argList[["keyType"]] <- keyType
  argList[["orgDb"]] <- orgDb
  argList[["madeUpper"]] <- c("dna","rna",
                                  tolower(AnnotationDbi::keys(orgDb,
                                                              keytype="SYMBOL")))
  wc <- do.call("wcGeneSummary",argList)
  barp <- utils::head(wc@freqDf, n=top)
  ## Need rewrite
  if (reord){
    if (grad) {
      plt <- ggplot(barp, aes(x=reorder(barp$word, barp$freq),
                              y=barp$freq, fill=barp$freq))+ scale_fill_gradient(low="blue",high="red", guide="none")
    } else {
      plt <- ggplot(barp, aes(x=reorder(barp$word, barp$freq),
                              y=barp$freq, fill=barp$word))+ scale_fill_manual(values=pal, guide="none")
    }
  } else {
    if (grad) {
      plt <- ggplot(barp, aes(x=barp$word,
                              y=barp$freq, fill=barp$freq))+ scale_fill_gradient(low="blue",high="red", guide="none")
    } else {
      plt <- ggplot(barp, aes(x=barp$word,
                              y=barp$freq, fill=barp$word))+ scale_fill_manual(values=pal, guide="none")
    }
  }     
  ## Need rewrite
  if (flip) {
    plt <- plt +    
      geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
      theme_minimal() + 
      theme(axis.text = element_text(size = textSize))+coord_flip()
  } else {
    plt <- plt +    
      geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
      theme_minimal() + 
      theme(axis.text = element_text(size = textSize, angle=90))
  }
  if (retList){
    return(list(wc=wc, plot=plt))
  } else {
    return(plt)
  }
}


#' exportCyjs
#' 
#' Export Cytoscape.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @import jsonlite
#' @importFrom cyjShiny dataFramesToJSON
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \donttest{exportCyjs(g, "./", "net")}
#' @export
#' 
exportCyjs <- function(g, rootDir, netDir) {
    
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    
    nodes <- data.frame(
        id=names(V(g)),
        label=names(V(g)),
        image=V(g)$image,
        size=V(g)$size,
        shape=V(g)$shape
    )
    
    edgeList <- as_edgelist(g)
    edges <- data.frame(source=edgeList[,1],
        target=edgeList[,2], interaction=NA)
    edges$strength <- E(g)$strength
    
    pret <- prettify(dataFramesToJSON(edges, nodes))
    pret <- substr(pret, 18, nchar(pret)-3)
    pret
    
    js <- paste0("
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
          style: cytoscape.stylesheet()
            .selector('node')
            .css({
                      'content': 'data(label)',
                      'shape' : 'data(shape)',
                      'background-image': 'data(image)',
                      'text-valign': 'bottom',
                      'background-color': '#FFF',
                      'background-fit': 'cover',
                      'width': 'data(size)',
                      'height': 'data(size)',
                      'font-size' : 'mapData(size, 0, 100, 1, 20)',
                      'text-outline-width': 1,
                      'text-outline-color': '#FFF',
                      'border-color' : '#555',
                      'border-width': 1
                  })
            .selector('edge')
            .css({
                    'width' : '4',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    'width' : 'mapData(strength, 0.5, 1, 0, 5)'
                  }),
        'elements':
    ", pret, ",
        layout:{
              name: 'cola',
              padding: 0.5,
              avoidOverlap: true, 
              nodeSpacing: function( node ){ return 0.1; },
              nodeDimensionsIncludeLabels: true
          }
        });
    ")
    
    ## Using cola layout by default.
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="cy"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- "
    body {
        font-family: helvetica, sans-serif;
        font-size: 14px;
    }
    
    #cy {
        position: absolute;
        left: 0;
        top: 0;
        right: 0;
        bottom: 0;
        z-index: 999;
    }
    
    h1 {
        opacity: 0.5;
        font-size: 1em;
    }"
    
    
    write(js, file = paste0(rootDir, netDir, "/script.js"))
    write(style, file = paste0(rootDir, netDir, "/style.css"))
    write(html, file = paste0(rootDir, netDir, "/index.html"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' exportVisjs
#' 
#' Export vis.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @import jsonlite
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \donttest{exportVisjs(g, "./", "net")}
#' @export
#' 
exportVisjs <- function(g, rootDir, netDir){
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    if (unique(V(g)$shape)=="rectangle"){
        visjsShape <- "image"
    } else {
        visjsShape <- "circularImage"
    }
    
    nodejson <- toJSON(data.frame(
            id=names(V(g)),
            label=names(V(g)),
            image=V(g)$image,
            shape=visjsShape,
            size=V(g)$size
        ))
    
    edgeList <- as_edgelist(g)
    edgejson <- toJSON(data.frame(from=edgeList[,1], to=edgeList[,2],
        width=E(g)$strength))
    
    # Make JS
    js <- paste0("
    var nodes = null;
    var edges = null;
    var network = null;
    
    function draw() {
    nodes = ", nodejson, ";
    edges = ", edgejson, ";
      var container = document.getElementById('mynetwork');
      var data = {
        nodes: nodes,
        edges: edges,
      };
      var options = {
        nodes: {
          borderWidth: 2,
          size: 30,
          color: {
            border: '#222222',
    background: 'white',
    },
    font: { color: 'black' },
    },
    edges: {
        length: 200,
        color: 'lightgray',
        arrows: { to: {enabled: true} }
    },
    layout: {
        improvedLayout: true
    }
    };
    network = new vis.Network(container, data, options);
    }
    
    window.addEventListener('load', () => {
        draw();
    });
    ")
    
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="mynetwork"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- '
    body {
      font: 10pt arial;
    }
    #mynetwork {
      width: 1000px;
      height: 1000px;
      border: 1px solid lightgray;
      background-color: white;
    }
    '
    
    write(js, file = paste0(rootDir, netDir, "/script.js"))
    write(html, file = paste0(rootDir, netDir, "/index.html"))
    write(style, file = paste0(rootDir, netDir, "/style.css"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' returnExample
#' 
#' return an example dataset used in the analysis
#' 
#' @import org.Hs.eg.db
#' @return return example MEs and colors
#' @examples returnExample()
#' @export
#' 
returnExample <- function() {
    ## Simulate WGCNA results (three modules)
    ccls <- c()
    for (i in c(1,2,3,4,5,6,7,8,9)){
        ccls <- c(ccls, paste0("CCL",i))
    }
    cxcls <- c()
    for (i in c(1,2,3,5,6,8,9,10,11,12,13,14,16)){
        cxcls <- c(cxcls, paste0("CXCL",i))
    }
    erccs <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")

    CCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=ccls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    CXCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=cxcls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    ERCCensg <- AnnotationDbi::select(org.Hs.eg.db, keys=erccs,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL

    CCLensg <- CCLensg[!is.na(CCLensg)]
    CXCLensg <- CXCLensg[!is.na(CXCLensg)]
    ERCCensg <- ERCCensg[!is.na(ERCCensg)]

    CCLcol <- rep(1, length(CCLensg))
    names(CCLcol) <- CCLensg
    CXCLcol <- rep(2, length(CXCLensg))
    names(CXCLcol) <- CXCLensg
    ERCCcol <- rep(3, length(ERCCensg))
    names(ERCCcol) <- ERCCensg
    modColors <- c(CCLcol, CXCLcol, ERCCcol)
    ensg <- names(modColors)

    MEs <- data.frame(
        ME1 = 1:10,
        ME2 = 2:11,
        ME3 = 10:1
    )

    mod <- list()
    mod[["MEs"]] <- MEs
    mod[["colors"]] <- modColors
    mod
}

#' makeCorpus
#' 
#' Clean-up the corpus
#' 
#' @param docs corpus to clean
#' @param filterWords words to filter based on frequency
#' @param additionalRemove words to filter
#' @param numOnly delete number only
#' @param stem use stem or not
#' @param lower transform lower for corpus
#' 
#' @return cleaned corpus
#' @import tm
#' 
#' 
makeCorpus <- function (docs, filterWords, additionalRemove, numOnly, stem, lower=TRUE) {
    if (lower) {
        docs <- docs %>%
            tm_map(FUN=content_transformer(tolower))
    }
    if (numOnly) {
        docs <- docs %>% tm_map(FUN=removeAloneNumbers)
    } else {
        docs <- docs %>% tm_map(FUN=removeNumbers)
    }
    docs <- docs %>%
        tm_map(removeWords, stopwords::stopwords("english",
            "stopwords-iso")) %>%
        tm_map(removeWords, filterWords) %>% 
        tm_map(FUN=removePunctuation) %>%
        tm_map(FUN=stripWhitespace)
    if (prod(is.na(additionalRemove))!=1){
        docs <- docs %>% tm_map(removeWords, additionalRemove)
    }
    if (stem) {
        docs <- docs %>% tm_map(stemDocument)   
    }
    return(docs)
}


#' getUPtax
#' 
#' Obtain the list of UniProt organism identification codes, by querying taxonomy or UniProt codes.
#' https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/docs/speclist.txt
#' The file above must be downloaded and specified to file argument.
#' 
#' @param file downloaded file
#' @param candUP candidate UniProt organism identification codes
#' @param candTax candidate taxonomy name
#' @return data.frame consisting of taxonomy name and UniProt IDs
#' @examples
#' file <- "speclist.txt"
#' \donttest{getUPtax(file, candUP="all")}
#' @export
#' 
getUPtax <- function(file, candUP, candTax=NULL) {
  tmp <- NULL
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (grepl("N=",line)) {
      code <- strsplit(line, " ")[[1]][1]
      if (code=="Code") {next}
      sn <- unlist(strsplit(line, "N="))[2]
      if (length(candUP)==1) {
        if (candUP=="all") {
          if (!is.null(candTax)) {
            for (ct in candTax) {
              if (grepl(ct, sn)){
                tmp <- rbind(tmp, c(code, sn))
              }
            }
          } else {
            tmp <- rbind(tmp, c(code, sn))
          }
        }
      }
      if (code %in% candUP) {
        tmp <- rbind(tmp, c(code, sn))
      }
    }
  }
  close(con)
  if (!is.null(tmp)) {
    tmp <- data.frame(tmp) |> `colnames<-`(c("UPID","Taxonomy"))
    return(tmp)
  } else {
    return(NULL)
  }
}




#' exportCyjsWithoutImage
#' 
#' Export Cytoscape.js script, HTML and stylesheet for the graph without image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @param edgeWidth attribute name for edgeWidth
#' @param nodeColor attribute name for node color
#' @param nodeSize attribute name for node size
#' @param nodeColorDiscretePal color mapping palette for discrete variables
#' @param sizeMin minimum size for scaling node size
#' @param sizeMax maximum size for scaling node size
#' @import jsonlite
#' @importFrom cyjShiny dataFramesToJSON
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$size <- c(1,1)
#' \donttest{exportCyjsWithoutImage(g, "./", "net")}
#' @export
#' 
exportCyjsWithoutImage <- function(g, rootDir, netDir,
  edgeWidth="weight",nodeColor="tag",nodeSize="Freq",nodeColorDiscretePal="RdBu",
  sizeMin=10, sizeMax=50) {
    
    if (is.null(V(g)$shape)){qqcat("No node shape specified, set to 'circle'\n")
      V(g)$shape <- rep("circle", length(V(g)))}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}


    ## Make node color
    nc <- get.vertex.attribute(g, nodeColor)
    nodeColors <- RColorBrewer::brewer.pal(length(unique(nc)), nodeColorDiscretePal)
    names(nodeColors) <- unique(nc)

    ## Make node size
    V(g)$size <- get.vertex.attribute(g, nodeSize)
    rawMin <- min(V(g)$size)
    rawMax <- max(V(g)$size)
    scf <- (sizeMax-sizeMin)/(rawMax-rawMin)
    V(g)$size <- scf * V(g)$size + sizeMin - scf * rawMin

    nodes <- data.frame(
        id=names(V(g)),
        label=names(V(g)),
        size=V(g)$size,
        color=nodeColors[nc],
        shape=V(g)$shape
    )
    
    edgeList <- as_data_frame(g)
    edgeListRename <- colnames(edgeList)
    edgeListRename[1] <- "source";edgeListRename[2] <- "target"
    colnames(edgeList) <- edgeListRename
    edgeList$interaction <- rep(NA, nrow(edgeList))
    edgeList$width <- edgeList[[edgeWidth]]
    # edges <- data.frame(source=edgeList[,1],
    #     target=edgeList[,2], interaction=NA)
    # edges$strength <- E(g)$strength
    
    pret <- prettify(dataFramesToJSON(edgeList, nodes))
    pret <- substr(pret, 18, nchar(pret)-3)
    pret
    
    js <- paste0("
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
          style: cytoscape.stylesheet()
            .selector('node')
            .css({
                      'content': 'data(label)',
                      'shape' : 'data(shape)',
                      'text-valign': 'bottom',
                      'background-color': 'data(color)',
                      'background-fit': 'cover',
                      'width': 'data(size)',
                      'height': 'data(size)',
                      'font-size' : 'mapData(size, 0, 100, 10, 20)',
                      'text-outline-width': 1,
                      'text-outline-color': '#FFF',
                      'border-color' : '#555',
                      'border-width': 1
                  })
            .selector('edge')
            .css({
                    'curve-style': 'bezier',
                    'width' : 'mapData(width, 0.5, 1, 1, 5)'
                  }),
        'elements':
    ", pret, ",
        layout:{
              name: 'cola',
              padding: 0.5,
              avoidOverlap: true, 
              nodeSpacing: function( node ){ return 0.1; },
              nodeDimensionsIncludeLabels: true
          }
        });
    ")
    
    ## Using cola layout by default.
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="cy"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- "
    body {
        font-family: helvetica, sans-serif;
        font-size: 14px;
    }
    
    #cy {
        position: absolute;
        left: 0;
        top: 0;
        right: 0;
        bottom: 0;
        z-index: 999;
    }
    
    h1 {
        opacity: 0.5;
        font-size: 1em;
    }"
    
    
    write(js, file = paste0(rootDir, "/",netDir, "/script.js"))
    write(style, file = paste0(rootDir, "/",netDir, "/style.css"))
    write(html, file = paste0(rootDir, "/",netDir, "/index.html"))
    # message(paste0("Exported to ",rootDir,netDir))
}