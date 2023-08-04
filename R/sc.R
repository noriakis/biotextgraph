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
#' @param args parameters to pass to refseq
#' @param raw obtain raw results of refseq instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud,
#' named list of clusters
#' @param pvalThresh adjusted p-value threshold for markers
#' @param withTitle plot title on the plot
#' @param withggfx applying ggfx filters
#' @param ggfxParams parameter list for ggfx
#' @export
#' @examples 
#' markers <- data.frame(p_val_adj=c(0.01, 0.01, 0.04),
#' gene=c("PNKP","DDX41","IRF3"),cluster=c("1","1","1"))
#' colors <- list("1"="red")
#' TextMarkers(markers, col=colors, type="wc")
#' @return list of plots on textual information in the gene cluster
TextMarkers <- function(df, keyType="SYMBOL",type="wc", genePlot=TRUE,
         genePlotNum=5, colorText=TRUE, args=list(), wcArgs=NULL,
         raw=FALSE, col=NULL, pvalThresh=0.05, withTitle=TRUE,withggfx=NULL,
         ggfxParams=list()) {

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

        wres <- do.call(refseq, args)

        if (type=="wc"){
            if (!is.null(withggfx)) {
              wc <- do.call(eval(parse(text=withggfx)),
                c(list(
                  x = wres@wc
                  ),
                  ggfxParams)
              )
            } else {
              wc <- wres@wc
            }
            if (withTitle) {
                plotList[[i]] <- wc + ggtitle(i)
            } else {
                plotList[[i]] <- wc
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
#' @param args parameters to pass to refseq
#' @param raw obtain raw results of refseq instead of plot
#' @param colorText colorlize text or not, default to TRUE
#' @param col color to be used in wordcloud
#' @param top Top-{top} genes for each cluster will be included
#' @param withTitle ggtitle({cluster name}) will be added
#' @param withggfx applying ggfx filters
#' @param ggfxParams parameter list for ggfx
#' @param FDRThresh FDR threshold
#' @examples 
#' df <- data.frame(
#'   p.value=c(0.01, 0.01),gene=c("PNKP","DDX41"),
#'   Top=c(1,2)
#' )
#' row.names(df) <- df$gene
#' markers <- list("1"=df)
#' colors <- list("1"="blue")
#' TextMarkersScran(markers, col=colors)
#'
#' @return list of ggplot
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
    top=10,
    withggfx=NULL,
    FDRThresh=0.05,
    ggfxParams=list()
    ) {

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
        candidate <- subset(candidate, candidate$FDR < FDRThresh)
        cand_genes <- candidate |> row.names()
        args[["geneList"]] <- cand_genes

        wres <- do.call(refseq, args)

        if (type=="wc"){
            if (!is.null(withggfx)) {
              wc <- do.call(eval(parse(text=withggfx)),
                c(list(
                  x = wres@wc
                  ),
                  ggfxParams)
              )
            } else {
              wc <- wres@wc
            }
            if (withTitle) {
                plotList[[i]] <- wc + ggtitle(i)
            } else {
                plotList[[i]] <- wc
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
#' @param args parameters to passed to refseq
#' @param rot.per ggwordcloud parameter
#' @param dimred dimension reduction method
#' @param random.order ggwordcloud parameter
#' @param rad named vector of size of each cluster
#' @param top Top-{top} genes are included
#' @param sortBy default to avg_log2FC, "log10p" can be specified.
#' @param scaleNumber scale the frequency of words by this number
#' in `gene_name`
#' @param decreasing sort by decreasing order or not
#' @param geneNum number of genes to be included in wordclouds
#' @param gene_name show gene names instead of textual information
#' @param base_ellipse if TRUE, wordclouds are placed based on \code{stat_ellipse}.
#' @param base_dens if TRUE, wordclouds are placed based on density
#' @param withggfx applying ggfx filters
#' @param ggfxParams parameter list for ggfx
#' @export
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @return single-cell plot with text annotation
plotReducedDimWithTexts <- function(sce, marker.info,
         colour_by="label", point_alpha=0.4, use_shadowtext=TRUE,
         bg.colour="white", which.label=NULL, wc_alpha=1, wcScale=5,
         rot.per=0.4, rad=NULL, top=10, gene_name=FALSE, withggfx=NULL, ggfxParams=list(),
         sortBy="summary.logFC", scaleNumber=2, decreasing=TRUE, geneNum=50,
         random.order=FALSE, dimred="PCA", base_ellipse=FALSE, base_dens=FALSE,
         withTitle=FALSE, args=list()) {
    if (!use_shadowtext) {
      bg.colour <- NULL
    }

    args[["wcScale"]] <- wcScale
    if (base_dens) {base_ellipse <- TRUE}
    if (requireNamespace("scater", quietly = TRUE)) {## pass the plot itself
        rawPlot <- scater::plotReducedDim(sce, dimred=dimred,
                                  colour_by=colour_by,
                                  point_alpha=point_alpha)
    } else {
        stop("Please install scater")
    }
    g <- ggplot_build(rawPlot)

    ## Map the group
    map_group <- as.character(rawPlot$data$colour_by)
    names(map_group) <- as.character(g$data[[1]]$group)
    map_group <- map_group[!duplicated(map_group)]

    ## Obtain color and generate colors for wc
    ## Name as character
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

    wcArgs <- list(alpha=wc_alpha,rot.per=rot.per,random.order=random.order,
        bg.colour=bg.colour)
    if (gene_name) {
        subset.marker.info <- marker.info[which.label]
        texts <-  obtainMarkersWCScran(subset.marker.info,
                                    cols=cols,
                                    wcArgs=wcArgs,
                                    wcScale=wcScale,
                                    scaleNumber=scaleNumber,
                                    sortBy=sortBy,
                                    decreasing=decreasing,
                                    geneNum=geneNum,
                                    withggfx=withggfx,
                                    ggfxParams=ggfxParams)
    } else {

        texts <- marker.info[which.label] |> TextMarkersScran(wcArgs=wcArgs,
                                                              col=cols,top=top,
                                                              genePlot=FALSE,
                                                              args=args,
                                                              withTitle=withTitle,
                                                              withggfx=withggfx,
                                                              ggfxParams=ggfxParams)
    }

    if (base_ellipse) {
        el <- ggplot_build(rawPlot + 
                stat_ellipse(aes(x=.data$X,
                    y=.data$Y, group=.data$colour_by)))
        pl <- el$data[[1]]
        el <- el$data[[2]]

        new_points <- NULL
        for (i in unique(el$group)) {## el
          i <- as.character(i)
          tmp_el <- subset(el, el$group==i)[,c("x","y")]
          
          ctr = MASS::cov.trob(tmp_el)$center
          dist2center <- sqrt(rowSums((t(t(tmp_el)-ctr))^2))
          if (is.null(rad)) {
            ar <- pi*min(dist2center)*max(dist2center)
            r <- sqrt(ar / pi)
          } else {
            r <- rad[i]
          }
          
          if (base_dens) {## pl
            tmp_pl <- subset(pl, pl$group==i)[,c("x","y")]
            kde <- MASS::kde2d(tmp_pl$x, tmp_pl$y, n=100)
            ix <- findInterval(tmp_pl$x, kde$x)
            iy <- findInterval(tmp_pl$y, kde$y)
            ii <- cbind(ix, iy)
            tmp_pl$dens <- kde$z[ii]
            dens_max <- tmp_pl[order(tmp_pl$dens, decreasing=TRUE),][1,]
            XMe <- dens_max$x
            YMe <- dens_max$y
            
            new_points <- rbind(new_points,
                                c(map_group[i],
                                  XMe - r,
                                  YMe - r,
                                  XMe + r,
                                  YMe + r))
          } else {
            new_points <- rbind(new_points,
                                c(map_group[i],
                                  ctr["x"] - min(dist2center),
                                  ctr["y"] - min(dist2center),
                                  ctr["x"] + max(dist2center),
                                  ctr["y"] + max(dist2center)))
          }
          new_points <- data.frame(new_points) |>
              `colnames<-`(c("colour_by","XMi","YMi","XMa","YMa"))
        }
      } else {
        new_points <- rawPlot$data |>
          group_by(.data$colour_by) |>
          summarise(XMi=min(.data$X),
                    YMi=min(.data$Y),
                    XMa=max(.data$X),
                    YMa=max(.data$Y),
                    XMe=mean(.data$X),
                    YMe=mean(.data$Y))

        if (!is.null(rad)) {
            new_points <- data.frame(t(apply(new_points,1,function(x){
                xme <- as.numeric(x["XMe"])
                yme <- as.numeric(x["YMe"])
                c(x["colour_by"],
                  xme - rad[x["colour_by"]],
                  yme - rad[x["colour_by"]],
                  xme + rad[x["colour_by"]],
                  yme + rad[x["colour_by"]])
            }))) |> `colnames<-`(colnames(new_points)[1:5])
        }
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
#' @param args parameters to passed to refseq
#' @param rot.per ggwordcloud parameter
#' @param random.order ggwordcloud parameter
#' @param rad named vector of size of each cluster
#' @param sortBy default to avg_log2FC, "log10p" can be specified.
#' @param scaleNumber scale the frequency of words by this number
#' in `gene_name`
#' @param decreasing sort by decreasing order or not
#' @param geneNum number of genes to be included in wordclouds
#' @param gene_name show gene names instead of textual information
#' @param base_ellipse if TRUE, wordclouds are placed based on \code{stat_ellipse}.
#' @param base_dens if TRUE, wordclouds are placed based on density
#' @param withggfx applying ggfx filters
#' @param ggfxParams parameter list for ggfx
#' @export
#' @importFrom dplyr summarise
#' @importFrom dplyr group_by
#' @return single-cell plot with text annotation
DimPlotWithTexts <- function(seu, markers,
         label=TRUE, pt.size=0.5, reduction="umap",
         point_alpha=0.2, use_shadowtext=TRUE,
         bg.colour="white", which.label=NULL,
         wc_alpha=1, wcScale=5,
         rot.per=0.4, rad=NULL, sortBy="avg_log2FC", scaleNumber=2,
         decreasing=TRUE, geneNum=50, base_ellipse=FALSE, base_dens=FALSE,
         random.order=FALSE, gene_name=FALSE, withggfx=NULL, ggfxParams=list(),
         withTitle=FALSE, args=list()) {
    if (!use_shadowtext) {
      bg.colour <- NULL
    }
    args[["wcScale"]] <- wcScale
    if (base_dens) {base_ellipse <- TRUE}
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

    map_group <- as.character(plt$data$ident)
    names(map_group) <- as.character(g$data[[1]]$group)
    map_group <- map_group[!duplicated(map_group)]

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
    wcArgs <- list(alpha=wc_alpha,rot.per=rot.per,random.order=random.order,
        bg.colour=bg.colour)
    if (gene_name) {
        subset.markers <- subset(markers, markers$cluster %in% which.label)
        texts <- obtainMarkersWC(subset.markers,
                                cols=cols,
                                wcArgs=wcArgs,
                                wcScale=wcScale,
                                scaleNumber=scaleNumber,
                                sortBy=sortBy,
                                decreasing=decreasing,
                                geneNum=geneNum,
                                withggfx=withggfx,
                                ggfxParams=ggfxParams
                                )
    } else {
        texts <- subset(markers, markers$cluster %in% which.label) |> TextMarkers(
                                            wcArgs=wcArgs,
                                              col=cols,
                                              args=args,
                                              genePlot=FALSE,
                                              withTitle=withTitle,
                                              withggfx=withggfx,
                                              ggfxParams=ggfxParams
                                              )        
    }

  if (base_ellipse) {
    plt$data$x <- plt$data[,1]
    plt$data$y <- plt$data[,2]

    el <- ggplot_build(plt + 
            stat_ellipse(aes(x=.data$x, y=.data$y, group=.data$ident)))
    pl <- el$data[[1]]
    el <- el$data[[2]]

    new_points <- NULL
    for (i in unique(el$group)) {## el
      i <- as.character(i)
      tmp_el <- subset(el, el$group==i)[,c("x","y")]
      
      ctr = MASS::cov.trob(tmp_el)$center
      dist2center <- sqrt(rowSums((t(t(tmp_el)-ctr))^2))
      if (is.null(rad)) {
        ar <- pi*min(dist2center)*max(dist2center)
        r <- sqrt(ar / pi)
      } else {
        r <- rad[i]
      }
      
      if (base_dens) {## pl
        tmp_pl <- subset(pl, pl$group==i)[,c("x","y")]
        kde <- MASS::kde2d(tmp_pl$x, tmp_pl$y, n=100)
        ix <- findInterval(tmp_pl$x, kde$x)
        iy <- findInterval(tmp_pl$y, kde$y)
        ii <- cbind(ix, iy)
        tmp_pl$dens <- kde$z[ii]
        dens_max <- tmp_pl[order(tmp_pl$dens, decreasing=TRUE),][1,]
        XMe <- dens_max$x
        YMe <- dens_max$y
        
        new_points <- rbind(new_points,
                            c(map_group[i],
                              XMe - r,
                              YMe - r,
                              XMe + r,
                              YMe + r))
      } else {
        new_points <- rbind(new_points,
                            c(map_group[i],
                              ctr["x"] - min(dist2center),
                              ctr["y"] - min(dist2center),
                              ctr["x"] + max(dist2center),
                              ctr["y"] + max(dist2center)))
      }
      new_points <- data.frame(new_points) |>
          `colnames<-`(c("ident","XMi","YMi","XMa","YMa"))
    }
  } else {
        
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
        
        if (!is.null(rad)) {
            new_points <- data.frame(t(apply(new_points,1,function(x){
                xme <- as.numeric(x["XMe"])
                yme <- as.numeric(x["YMe"])
                c(x["ident"],
                  xme - rad[x["ident"]],
                  yme - rad[x["ident"]],
                  xme + rad[x["ident"]],
                  yme + rad[x["ident"]])
            }))) |> `colnames<-`(colnames(new_points)[1:5])
        }
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