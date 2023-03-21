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
        i <- as.character(i)
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
        cand_genes <- candidate$gene
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
        i <- as.character(i)
        qqcat("@{i}\n")
        if (is.null(col)){
            tmpcol <- RColorBrewer::brewer.pal(8, sample(
                row.names(RColorBrewer::brewer.pal.info),1))
        } else {
            tmpcol <- col[[i]]
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



#' plotReducedDimWithTexts
#' 
#' @param sce sce object
#' @param marker.info results of findMarkers()
#' @param colour_by colorize by this label
#' @param point_alpha point alpha
#' @param use_shadowtext use shadowtext for wordcloud
#' @param bg.colour shadowtext background color
#' @param which.label which label to plot text
#' @param wc_alpha alpha value for wordcloud
#' @param wcScale scaling value for wordcloud
#' @param withTitle whether to append title on wordcloud
#' @param args parameters to passed to wcGeneSummary
#' @param rot.per ggwordcloud parameter
#' @param dimred dimension reduction method
#' @param random.order ggwordcloud parameter
#' @param r named vector of size of each cluster
#' @param top Top-{top} genes are included
#' @export
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @return single-cell plot with text annotation
plotReducedDimWithTexts <- function(sce, marker.info,
         colour_by="label", point_alpha=0.4, use_shadowtext=TRUE,
         bg.colour="white", which.label=NULL, wc_alpha=1, wcScale=5,
         rot.per=0.4, r=NULL, top=10,
         random.order=FALSE, dimred="PCA",
         withTitle=FALSE, args=list()) {
    args[["wcScale"]] <- wcScale
    if (requireNamespace("scater", quietly = TRUE)) {## pass the plot itself
        rawPlot <- scater::plotReducedDim(sce, dimred=dimred,
                                  colour_by=colour_by,
                                  point_alpha=point_alpha)
    } else {
        stop("Please install scater")
    }

    ## Obtain color and generate colors for wc
    ## Name as character
    g <- ggplot_build(rawPlot)
    colmap <- g$data[[1]][,c("colour","group")]
    colmap <- colmap[!duplicated(colmap),]
    row.names(colmap) <- colmap$group
    cols <- list()
    vec <- NULL
    for (i in seq_len(nrow(colmap))) {
      cols[[as.character(colmap[i,"group"])]] <- 
          colorRampPalette(c("grey",colmap[i,"colour"]))(10)
      vec[as.character(colmap[i,"group"])] <- colmap[i,"colour"]
    }

    if (is.null(which.label)) {
        which.label <- names(marker.info)
    }
    texts <- marker.info[which.label] |> TextMarkersScran(wcArgs=list(alpha=wc_alpha,
                                                                      rot.per=rot.per,
                                                                      random.order=random.order,
                                                                      use_shadowtext=use_shadowtext,
                                                                      bg.colour=bg.colour),
                                                          col=cols,top=top,
                                                          genePlot=FALSE,
                                                          args=args,
                                                          withTitle=withTitle)


    new_points <- rawPlot$data |>
      group_by(.data$colour_by) |>
      summarise(XMi=min(.data$X),
                YMi=min(.data$Y),
                XMa=max(.data$X),
                YMa=max(.data$Y),
                XMe=mean(.data$X),
                YMe=mean(.data$Y))

    if (!is.null(r)) {
        new_points <- data.frame(t(apply(new_points,1,function(x){
            xme <- as.numeric(x["XMe"])
            yme <- as.numeric(x["YMe"])
            c(x["colour_by"],
              xme - r[x["colour_by"]],
              yme - r[x["colour_by"]],
              xme + r[x["colour_by"]],
              yme + r[x["colour_by"]])
        }))) |> `colnames<-`(colnames(new_points)[1:5])
    }

    for (i in names(texts)) {
      tmp <- subset(new_points,
                    new_points$colour_by==i)
      tmpXMi <- as.numeric(tmp$XMi);
      tmpYMi <- as.numeric(tmp$YMi);
      tmpXMa <- as.numeric(tmp$XMa);
      tmpYMa <- as.numeric(tmp$YMa);
      rawPlot <- rawPlot + annotation_custom(ggplotify::as.grob(texts[[i]]),
                                       xmin=tmpXMi, xmax=tmpXMa,
                                       ymin=tmpYMi, ymax=tmpYMa)
    }
    rawPlot
}





#' DimPlotWithTexts
#' 
#' @param seu Seurat object
#' @param markers results of FindAllMarkers()
#' @param label plot label or not
#' @param pt.size point size in plot
#' @param reduction reduction method
#' @param point_alpha point alpha
#' @param use_shadowtext use shadowtext for wordcloud
#' @param bg.colour shadowtext background color
#' @param which.label which label to plot text
#' @param wc_alpha alpha value for wordcloud
#' @param wcScale scaling value for wordcloud
#' @param withTitle whether to append title on wordcloud
#' @param args parameters to passed to wcGeneSummary
#' @param rot.per ggwordcloud parameter
#' @param random.order ggwordcloud parameter
#' @param r named vector of size of each cluster
#' @export
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @return single-cell plot with text annotation
DimPlotWithTexts <- function(seu, markers,
         label=TRUE, pt.size=0.5, reduction="umap",
         point_alpha=0.2, use_shadowtext=TRUE,
         bg.colour="white", which.label=NULL,
         wc_alpha=1, wcScale=5,
         rot.per=0.4, r=NULL,
         random.order=FALSE,
         withTitle=FALSE, args=list()) {
    args[["wcScale"]] <- wcScale
    
    if (requireNamespace("Seurat", quietly = TRUE)) {## pass the plot itself
        plt <- Seurat::DimPlot(seu, reduction = reduction,
                label = label, pt.size = pt.size) 
        # https://github.com/satijalab/seurat/issues/2835
        plt[[1]]$layers[[1]]$aes_params$alpha <- point_alpha
    } else {
        stop("Please install Seurat")
    }

    ## Obtain color and generate colors for wc
    ## Name as character
    g <- ggplot_build(plt)
    g$data[[1]]$group <- plt$data$ident
    colmap <- g$data[[1]][,c("colour","group")]
    colmap <- colmap[!duplicated(colmap),]
    row.names(colmap) <- colmap$group

    cols <- list()
    vec <- NULL
    for (i in seq_len(nrow(colmap))) {
      cols[[as.character(colmap[i,"group"])]] <- 
          colorRampPalette(c("grey",colmap[i,"colour"]))(10)
      vec[as.character(colmap[i,"group"])] <- colmap[i,"colour"]
    }

    if (is.null(which.label)) {
        which.label <- unique(markers$cluster)
    }

    texts <- subset(markers, markers$cluster %in% which.label) |> TextMarkers(
                                        wcArgs=list(alpha=wc_alpha,
                                        rot.per=rot.per,
                                        random.order=random.order,
                                        use_shadowtext=use_shadowtext,
                                        bg.colour=bg.colour),
                                          col=cols,
                                          args=args,
                                          genePlot=FALSE,
                                          withTitle=withTitle)

    plt$data$X <- plt$data[,1]
    plt$data$Y <- plt$data[,2]
    new_points <- plt$data |>
      group_by(.data$ident) |>
      summarise(XMi=min(.data$X),
                YMi=min(.data$Y),
                XMa=max(.data$X),
                YMa=max(.data$Y),
                XMe=mean(.data$X),
                YMe=mean(.data$Y))
    
    if (!is.null(r)) {
        new_points <- data.frame(t(apply(new_points,1,function(x){
            xme <- as.numeric(x["XMe"])
            yme <- as.numeric(x["YMe"])
            c(x["ident"],
              xme - r[x["ident"]],
              yme - r[x["ident"]],
              xme + r[x["ident"]],
              yme + r[x["ident"]])
        }))) |> `colnames<-`(colnames(new_points)[1:5])
    }


    for (i in names(texts)) {
      tmp <- subset(new_points,
                    new_points$ident==i)
      tmpXMi <- as.numeric(tmp$XMi);
      tmpYMi <- as.numeric(tmp$YMi);
      tmpXMa <- as.numeric(tmp$XMa);
      tmpYMa <- as.numeric(tmp$YMa);
      plt <- plt + annotation_custom(ggplotify::as.grob(texts[[i]]),
                                       xmin=tmpXMi, xmax=tmpXMa,
                                       ymin=tmpYMi, ymax=tmpYMa)
    }
    plt
}