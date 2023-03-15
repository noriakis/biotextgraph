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
#' @param wcArg parameters to pass to ggwordcloud
#' @param raw obtain raw results of wcGeneSummary instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud
#' @export
TextMarkers <- function(df, keyType="SYMBOL",type="wc", genePlot=TRUE,
         genePlotNum=5, colorText=TRUE, wcArg=NULL, raw=FALSE, col=NULL) {
    plotList <- list()
    wresList <- list()
    if (is.null(wcArg)) {
        wcArg <- list(
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
            tmpcol <- col[i]
        }
        wcArg[["colors"]] <- tmpcol
        candidate <- subset(df, df$p_val_adj<0.05 & df$cluster==i)
        cand_genes <- candidate |> row.names()
        wres <- wcGeneSummary(cand_genes, keyType=keyType, useggwordcloud = TRUE,
                              genePlot=genePlot, genePlotNum=genePlotNum,
                              colorText=colorText, plotType=type,
                              argList=wcArg)
        if (type=="wc"){
            plotList[[i]] <- wres@wc + ggtitle(i)
        } else {
            plotList[[i]] <- wres@net + ggtitle(i)
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
#' @param wcArg parameters to pass to ggwordcloud
#' @param raw obtain raw results of wcGeneSummary instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud
#' @param top Top-{top} genes will be included
#' @export
TextMarkersScran <- function(res, keyType="SYMBOL",type="wc", genePlot=TRUE,
         genePlotNum=5, colorText=TRUE, wcArg=NULL, raw=FALSE, col=NULL, top=10) {
    plotList <- list()
    wresList <- list()
    if (is.null(wcArg)) {
        wcArg <- list(
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

        wcArg[["colors"]] <- tmpcol
        tmp <- res[[i]]
        candidate <- subset(tmp, tmp$Top < top)
        cand_genes <- candidate |> row.names()
        wres <- wcGeneSummary(cand_genes, keyType=keyType, useggwordcloud = TRUE,
                              genePlot=genePlot, genePlotNum=genePlotNum,
                              colorText=colorText, plotType=type,
                              argList=wcArg)
        if (type=="wc"){
            plotList[[i]] <- wres@wc + ggtitle(i)
        } else {
            plotList[[i]] <- wres@net + ggtitle(i)
        }
        wresList[[i]] <- wres
    }
    if (raw) {
        return(wresList)
    } else {
        return(plotList)
    }
}