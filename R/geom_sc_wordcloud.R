
#' obtainMarkersWC
#' @param markers marker data frame
#' @param cols list of colors
#' @param wcArgs arguments for ggwordcloud
#' @param sort_by default to avg_log2FC, "log10p" can be specified.
#' @param scale_number scale the frequency by this number
#' @param wcScale scale for scale_size_area()
#' @param decreasing sort by decreasing order or not
#' @param geneNum number of genes to be included in wordclouds
#' @param eps when taking log of p-values, this value will be added
#' @export
obtainMarkersWC <- function(markers,
                            cols,
                            wcArgs,
                            eps=1e-10,
                            sort_by="avg_log2FC",
                            decreasing=TRUE,
                            scale_number=10,
                            wcScale=5,
                            geneNum=50) {
    markers$log10p <- -1*log10(markers$p_val+eps)
    wcList <- list()
    for (i in unique(markers$cluster)) {
        tmp_markers <- subset(markers, markers$cluster==i)
        tmp_markers <- tmp_markers[order(tmp_markers[[sort_by]],
                                         decreasing=decreasing),]
        tmp_markers <- tmp_markers[1:geneNum,]
        wcArgs[["colors"]] <- cols[[as.character(i)]]
        wcArgs[["word"]] <- tmp_markers$gene
        wcArgs[["freq"]] <- tmp_markers[[sort_by]]*scale_number
        wc <- do.call(ggwordcloud::ggwordcloud, wcArgs)
        wc <- wc + scale_size_area(max_size = wcScale)
        wcList[[as.character(i)]] <- wc
    }
    wcList
}


#' obtainMarkersWCScran
#' make gene wordcloud from scran::findMarkers() results
#' @param markers marker list
#' @param cols list of colors
#' @param wcArgs arguments for ggwordcloud
#' @param sort_by default to summary.logFC, "log10p" can be specified.
#' @param scale_number scale the frequency by this number
#' @param wcScale scale for scale_size_area()
#' @param decreasing sort by decreasing order or not
#' @param geneNum number of genes to be included in wordclouds
#' @param eps when taking log of p-values, this value will be added
#' @export
obtainMarkersWCScran <- function(markers,
                            cols,
                            wcArgs,
                            eps=1e-10,
                            sort_by="summary.logFC",
                            decreasing=TRUE,
                            scale_number=10,
                            wcScale=5,
                            geneNum=50) {
    wcList <- list()
    for (i in names(markers)) {
        tmp_markers <- markers[[i]]
        tmp_markers$log10p <- -1*log10(tmp_markers$p.value+eps)
        tmp_markers <- tmp_markers[order(tmp_markers[[sort_by]],
                                         decreasing=decreasing),]

        tmp_markers <- tmp_markers[1:geneNum,]
        wcArgs[["colors"]] <- cols[[as.character(i)]]
        wcArgs[["word"]] <- tmp_markers |> row.names()
        wcArgs[["freq"]] <- tmp_markers[[sort_by]]*scale_number
        wc <- do.call(ggwordcloud::ggwordcloud, wcArgs)
        wc <- wc + scale_size_area(max_size = wcScale)
        wcList[[as.character(i)]] <- wc
    }
    wcList
}

#' ggplot_add.geom_sc_wordcloud
#' use ggplot_add to populate single-cell plot with textual information
#' @param object An object to add to the plot
#' @param plot The ggplot object to add object to
#' @param object_name The name of the object to add
#' @export ggplot_add.geom_sc_wordcloud
#' @export
ggplot_add.geom_sc_wordcloud <- function(object, plot, object_name) {
  
  if (is.null(object$show_markers)) {
    show_markers <- unique(object$markers$cluster)
  } else {
    show_markers <- object$show_markers
  }
  
  
  g <- ggplot_build(plot)
  
  map_group <- as.character(plot$data$ident)
  names(map_group) <- as.character(g$data[[1]]$group)
  map_group <- map_group[!duplicated(map_group)]
  
  g$data[[1]]$group <- plot$data$ident
  colmap <- g$data[[1]][,c("colour","group")]
  colmap <- colmap[!duplicated(colmap),]
  row.names(colmap) <- colmap$group
  
  cols <- list()
  for (i in seq_len(nrow(colmap))) {
    cols[[as.character(colmap[i,"group"])]] <- 
      colorRampPalette(c("grey",colmap[i,"colour"]))(10)
  }
  
  
  if (object$base_ellipse) {
    ## TODO: Better option for directly plotting textGrob()
    ## is needed, like calculating density of points and 
    ## show the high-frequent words in high-density area.
    el <- ggplot_build(plot + 
            stat_ellipse(aes(group=.data$ident)))
    pl <- el$data[[1]]
    el <- el$data[[2]]

    new_points <- NULL
    for (i in unique(el$group)) {## el
      i <- as.character(i)
      tmp_el <- subset(el, el$group==i)[,c("x","y")]
      
      ctr = MASS::cov.trob(tmp_el)$center
      dist2center <- sqrt(rowSums((t(t(tmp_el)-ctr))^2))
      if (is.null(object$rad)) {
        ar <- pi*min(dist2center)*max(dist2center)
        r <- sqrt(ar / pi)
      } else {
        r <- object$rad[i]
      }
      
      if (object$base_dens) {## pl
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
    }
    new_points <- data.frame(new_points) |>
      `colnames<-`(c("ident","XMi","YMi","XMa","YMa"))

  } else {
    
    dim1 <- colnames(plot$data)[1]
    dim2 <- colnames(plot$data)[2]
    new_points <- plot$data |>
      group_by(.data$ident) |>
      summarise(XMi=min(.data[[dim1]]),
                YMi=min(.data[[dim2]]),
                XMa=max(.data[[dim1]]),
                YMa=max(.data[[dim2]]),
                XMe=mean(.data[[dim1]]),
                YMe=mean(.data[[dim2]]))
    
    if (!is.null(object$rad)) {
      r <- object$rad
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
  }

  markers <- markers |> 
    dplyr::filter(.data$cluster %in% show_markers) |>
    dplyr::filter(.data$p_val_adj < object$p_val_adj_threshold)
  
  args <- c(object$args,wcScale=object$wcScale)
  wcArgs <- c(object$wcArgs, alpha=object$wc_alpha,
              rot.per=object$rot.per,
              random.order=object$random.order,
              use_shadowtext=object$use_shadowtext,
              bg.colour=object$bg.colour)
  

  if (object$gene_name) {
      wcMarkers <- suppressMessages(obtainMarkersWC(markers,
                                                cols=cols,
                                                wcArgs=wcArgs,
                                                wcScale=object$wcScale,
                                                scale_number=object$scale_number,
                                                sort_by=object$sort_by,
                                                decreasing=object$decreasing,
                                                geneNum=object$geneNum))
  } else {
      wcMarkers <- suppressMessages(TextMarkers(markers,
                                                keyType=object$keyType,
                                                type="wc",
                                                genePlot=FALSE,
                                                col=cols,
                                                withTitle=object$withTitle,
                                                args=args,
                                                wcArgs=wcArgs))      
  }
  
  for (i in names(wcMarkers)) {
    tmp <- subset(new_points,
                  new_points$ident==i)
    tmpXMi <- as.numeric(tmp$XMi);
    tmpYMi <- as.numeric(tmp$YMi);
    tmpXMa <- as.numeric(tmp$XMa);
    tmpYMa <- as.numeric(tmp$YMa);
    plot <- plot + annotation_custom(ggplotify::as.grob(wcMarkers[[i]]),
                                     xmin=tmpXMi, xmax=tmpXMa,
                                     ymin=tmpYMi, ymax=tmpYMa)
  }
  plot
}

#' geom_sc_wordcloud
#' @param markers FindAllMarkers() results
#' @param show_markers candidate clusters to be appear in plot, default to NULL,
#' which means all the clusters are plotted
#' @param gene_name show gene names instead of textual information
#' @param p_val_adj_threshold default to 0.05, used in subsetting marker data frame
#' using \code{p_val_adj}
#' @param keyType key type of those listed in the markers. default to SYMBOL.
#' @param rad named vector of radius for wordclouds. if specified, this value is 
#' used to define positions of annotation_custom.
#' @param wc_alpha alpha value for word cloud
#' @param withTitle whether to show title in each grob of wordcloud
#' @param rot.per ggwordcloud parameter
#' @param random.order ggwordcloud parameter
#' @param use_shadowtext use shadowtext for wordcloud
#' @param bg.colour the background color of shadowtext
#' @param wcScale word cloud scaling factor
#' @param args passed to wcGeneSummary()
#' @param wcArgs passed to ggwordcloud()
#' @param base_ellipse if TRUE, wordclouds are placed based on \code{stat_ellipse}.
#' @param base_dens if TRUE, wordclouds are placed based on density
#' @param sort_by default to avg_log2FC, "log10p" can be specified.
#' @param scale_number scale the frequency of words by this number
#' in `gene_name`
#' @param decreasing sort by decreasing order or not
#' @param geneNum number of genes to be included in wordclouds
#' @importFrom MASS cov.trob
#' @importFrom ggplot2 ggplot_add
#' @export
geom_sc_wordcloud <- function(markers,
                              show_markers=NULL,
                              gene_name=FALSE,
                              p_val_adj_threshold = 0.05,
                              keyType = "SYMBOL",
                              rad=NULL, wc_alpha=1, withTitle=FALSE,
                              rot.per=0.4, random.order=FALSE,
                              use_shadowtext=TRUE, bg.colour="white",
                              wcScale=4, args=list(), wcArgs=list(),
                              sort_by="avg_log2FC", scale_number=2,
                              decreasing=TRUE, geneNum=50,
                              base_ellipse=FALSE, base_dens=FALSE) {
  if (base_dens) {
    base_ellipse <- TRUE
  }
  structure(list(markers = markers,
                 show_markers = show_markers,
                 gene_name = gene_name,
                 p_val_adj_threshold = p_val_adj_threshold,
                 keyType = keyType,
                 rad=rad,
                 wc_alpha=wc_alpha,
                 rot.per=rot.per,
                 use_shadowtext=use_shadowtext,
                 bg.colour=bg.colour,
                 wcScale=wcScale,
                 withTitle=withTitle,
                 args=args,
                 wcArgs=wcArgs,
                 sort_by=sort_by,
                 decreasing=decreasing,
                 geneNum=geneNum,
                 scale_number=scale_number,
                 random.order=random.order,
                 base_ellipse=base_ellipse,
                 base_dens=base_dens),
            class = "geom_sc_wordcloud")
}

