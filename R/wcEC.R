#' wcEC
#' 
#' Query the Enzyme Comission number and obtain description,
#' and search pubmed for these enzymes and make word cloud and
#' correlation network. Need "enzyme.dat" from ExPASy
#' (https://enzyme.expasy.org/).
#' 
#' @param file file downloaded from expasy
#' @param ecnum candidate ecnum, like those obtained from eggNOG-mapper
#' @param onlyTerm only return quoted queries
#' @param ... passed to osplot(target="pubmed")
#' @export
#' 

wcEC <- function(file, ecnum, onlyTerm=FALSE, ...) {
  flg <- FALSE
  candecs <- NULL
  con = file(file, "r")
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    if (startsWith(line,"ID")) {
      ccs <- NULL
      ecs <- NULL
      ec <- gsub("ID   ","",line)
      if (ec %in% ecnum) {
        flg <- TRUE
      } else {
        flg <- FALSE
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
      if (startsWith(line, "//")) {
        flg <- FALSE
        ecs <- c(ec, de, paste0(ccs, collapse=" "))
        candecs <- rbind(candecs, ecs)
        }
    }
  }
  close(con)
  candecs <- data.frame(candecs) |>
    `colnames<-`(c("number","desc","comment"))
  quoted <- dQuote(candecs$desc,options(useFancyQuotes = FALSE))
  if (onlyTerm) {return(quoted)}
  abst <- osplot(target="pubmed",
                 quoted, ...)
  return(abst)
}
