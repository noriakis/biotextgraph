#' TextMarkers
#' 
#' Take results of Seurat::FindAllMarkers as input,
#' and plot wordcloud or network for all the clusters.
#' 
#' @param df result of FindAllMarkers()
#' @param keyType keytype
#' @param type wc or network
#' @param genePlot whether to plot relevant genes
#' @param genePlotNum number of genes to plot
#' @param wcArgs parameters to pass to ggwordcloud
#' @param args parameters to pass to wcGeneSummary
#' @param raw obtain raw results of wcGeneSummary instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud,
#' named list of clusters
#' @param pvalThresh adjusted p-value threshold for markers
#' @param withTitle plot title on the plot
#' @export
#' @return list of plots on textual information in the gene cluster
TextMarkers <- function(df, keyType="SYMBOL",type="wc", genePlot=TRUE,
         genePlotNum=5, colorText=TRUE, args=list(), wcArgs=NULL,
         raw=FALSE, col=NULL, pvalThresh=0.05, withTitle=TRUE) {

    plotList <- list()
    wresList <- list()
    args[["plotType"]] <- type
    args[["colorText"]] <- colorText
    args[["genePlotNum"]] <- genePlotNum
    args[["genePlot"]] <- genePlot
    args[["keyType"]] <- keyType

    if (is.null(wcArgs)) {
        wcArgs <- list(
            rot.per=0.4,
            random.order=FALSE
        )
    }
    for (i in unique(df$cluster)) {
        qqcat("@{i}\n")
        if (is.null(col)){
            tmpcol <- RColorBrewer::brewer.pal(8, sample(
                row.names(RColorBrewer::brewer.pal.info),1))
        } else {
            tmpcol <- col[[i]]
        }
        wcArgs[["colors"]] <- tmpcol
        args[["argList"]] <- wcArgs
        candidate <- subset(df, df$p_val_adj<pvalThresh & df$cluster==i)
        cand_genes <- candidate |> row.names()
        args[["geneList"]] <- cand_genes

        wres <- do.call(wcGeneSummary, args)

        if (type=="wc"){
            if (withTitle) {
                plotList[[i]] <- wres@wc + ggtitle(i)
            } else {
                plotList[[i]] <- wres@wc
            }

        } else {
            if (withTitle) {
                plotList[[i]] <- wres@net + ggtitle(i)
            } else {
                plotList[[i]] <- wres@net
            }
        }
                wresList[[i]] <- wres
        col <- NULL
    }
    if (raw) {
        return(wresList)
    } else {
        return(plotList)
    }
}


#' TextMarkersScran
#' 
#' Take results of scran::findMarkers as input,
#' and plot wordcloud or network for all the clusters.
#' 
#' @param res result of findMarkers()
#' @param keyType keytype
#' @param type wc or network
#' @param genePlot whether to plot relevant genes
#' @param genePlotNum number of genes to plot
#' @param wcArgs parameters to pass to ggwordcloud
#' @param args parameters to pass to wcGeneSummary
#' @param raw obtain raw results of wcGeneSummary instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud
#' @param top Top-{top} genes for each cluster will be included
#' @param withTitle ggtitle({cluster name}) will be added
#' @export
TextMarkersScran <- function(res,
    keyType="SYMBOL",
    type="wc",
    genePlot=TRUE,
    genePlotNum=5,
    colorText=TRUE,
    args=list(),
    wcArgs=NULL,
    raw=FALSE,
    col=NULL,
    withTitle=TRUE,
    top=10) {

    plotList <- list()
    wresList <- list()
    args[["plotType"]] <- type
    args[["colorText"]] <- colorText
    args[["genePlotNum"]] <- genePlotNum
    args[["genePlot"]] <- genePlot
    args[["keyType"]] <- keyType

    if (is.null(wcArgs)) {
        wcArgs <- list(
            rot.per=0.4,
            random.order=FALSE
        )
    }
    for (i in names(res)) {
        qqcat("@{i}\n")
        if (is.null(col)){
            tmpcol <- RColorBrewer::brewer.pal(8, sample(
                row.names(RColorBrewer::brewer.pal.info),1))
        } else {
            tmpcol <- col[i]
        }

        wcArgs[["colors"]] <- tmpcol
        args[["argList"]] <- wcArgs

        tmp <- res[[i]]
        candidate <- subset(tmp, tmp$Top < top)
        cand_genes <- candidate |> row.names()
        args[["geneList"]] <- cand_genes

        wres <- do.call(wcGeneSummary, args)

        if (type=="wc"){
            if (withTitle) {
                plotList[[i]] <- wres@wc + ggtitle(i)
            } else {
                plotList[[i]] <- wres@wc
            }

        } else {
            if (withTitle) {
                plotList[[i]] <- wres@net + ggtitle(i)
            } else {
                plotList[[i]] <- wres@net
            }
        }
        wresList[[i]] <- wres
    }
    if (raw) {
        return(wresList)
    } else {
        return(plotList)
    }
}