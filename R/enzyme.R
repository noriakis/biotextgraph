#' enzyme
#' 
#' @description Text mining the enzyme information.
#' 
#' @details Query the Enzyme Comission number and obtain description,
#' and search pubmed for these enzymes and make word cloud and
#' correlation network. Need to specify the path to "enzyme.dat"
#' downloaded from from ExPASy (https://enzyme.expasy.org/).
#' 
#' @param file file downloaded from expasy
#' @param ecnum candidate ecnum, like those obtained from eggNOG-mapper
#' @param onlyTerm only return quoted queries to wcAbst
#' @param onlyDf only return ec description data.frame
#' if onlyTerm and onlyDf are both specified, onlyTerm have priority
#' @param taxec link taxonomy to EC using UniProt Taxonomy ID file
#' If this is TRUE, data.frame is returned
#' @param taxFile UniProt organism ID file path
#' @param candTax when taxec=TRUE, search only for these species.
#' @param argList passed to osplot(target="pubmed")
#' @param target abstract or title
#' @param apiKey api key for PubMed
#' @return object consisting of data frame and ggplot2 object
#' @examples
#' file <- "enzyme.dat"
#' \dontrun{enzyme(file, ecnum="1.2.1.1")}
#' @export
#' @seealso generalf
enzyme <- function(file, ecnum, onlyTerm=FALSE, onlyDf=FALSE, target="abstract",
                 taxec=FALSE, taxFile=NULL, candTax=NULL, argList=list(),
                 apiKey=NULL) {
  flg <- FALSE
  candecs <- NULL
  allFlag <- FALSE
  if (length(ecnum)==1) {
    if (ecnum=="all") {
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
        if (ec %in% ecnum) {
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
  if (taxec) {
    qqcat("  Linking taxonomy to EC\n")
    retTaxEC <- NULL
    if (is.null(taxFile)) {stop("Please provide UniProt Taxonomy file")}
    if (!is.null(candTax)) {
      taxo <- getUPtax(taxFile, candUP="all", candTax=candTax)
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
      stop("Please specify candTax when taxec=TRUE")
    }
    if (is.null(retTaxEC)) {stop("No EC could be found for query")}
    retTaxEC <- data.frame(retTaxEC) |> 
      `colnames<-`(c("number","desc","taxonomy","scientificName"))
    queryCheck <- NULL
    for (ct in candTax) {
      queryCheck <- cbind(queryCheck,
                          grepl(ct, retTaxEC$scientificName))
    }
    retTaxEC$query <- apply(queryCheck, 1, function(x) {
      if (length(candTax[x])==1) {
        candTax[x]
      } else {
        paste0(candTax[x],",")
      }})
    return(retTaxEC)
  }
  
  quoted <- dQuote(candecs$desc,options(useFancyQuotes = FALSE))
  if (onlyTerm) {return(quoted)}
  if (onlyDf) {return(candecs)}
  argList[["target"]] <- "abstract"
  argList[["queries"]] <- quoted
  argList[["apiKey"]] <- apiKey
  abst <- do.call("pubmed", argList)
  abst@ec <- candecs
  abst@type <- "EC"
  return(abst)
}
