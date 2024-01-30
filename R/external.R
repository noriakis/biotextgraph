
#' refseqWGCNA
#' 
#' Return the list of biotext class object per cluster for the blockwise module results in WGCNA
#' 
#' @param wgcna results of blockwiseModules()
#' @param keyType key type of gene
#' @param argList argument list to pass to refseq()
#' @export
#' @return list of biotext class object
refseqWGCNA <- function(wgcna, keyType="ENSEMBL", argList=list()) {
    all_cols <- unique(wgcna$colors)
    qqcat("Processing a total of @{length(all_cols)} clusters\n")
    lapply(all_cols, function(col) {
        input <- wgcna$colors[wgcna$colors == col] |> names()
        argList[["geneList"]] <- input
        if (!("keyType" %in% names(argList))) {
            argList[["keyType"]] <- keyType
        }
        res <- try(do.call(refseq, argList))
        if ("try-error" %in% class(res)) {
            return(NA)
        } else { return(res) }
    })
}

#' refseqDESeq2
#' 
#' Return the biotext class object by specified filter in DESeq2 results object
#' 
#' @param res results of DESeq2::results()
#' @param condition filtering condition
#' @param keyType key type of gene
#' @param argList argument list to pass to refseq()
#' @export
#' @return list of biotext class object
refseqDESeq2 <- function(res, condition, keyType="ENSEMBL",
    argList=list()) {
    argList[["geneList"]] <- res %>%
        data.frame() %>%
        filter(!!enquo(condition)) %>%
        row.names()
    if (!("keyType" %in% names(argList))) {
        argList[["keyType"]] <- keyType
    }
    do.call(refseq, argList)
}
