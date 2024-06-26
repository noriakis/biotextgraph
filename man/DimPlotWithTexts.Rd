% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc.R
\name{DimPlotWithTexts}
\alias{DimPlotWithTexts}
\title{DimPlotWithTexts}
\usage{
DimPlotWithTexts(
  seu,
  markers,
  label = TRUE,
  pt.size = 0.5,
  reduction = "pca",
  point_alpha = 0.2,
  use_shadowtext = TRUE,
  bg.colour = "white",
  which.label = NULL,
  wc_alpha = 1,
  wcScale = 5,
  rot.per = 0.4,
  rad = NULL,
  sortBy = "avg_log2FC",
  scaleNumber = 2,
  decreasing = TRUE,
  geneNum = 50,
  base_ellipse = FALSE,
  base_dens = FALSE,
  random.order = FALSE,
  gene_name = FALSE,
  withggfx = NULL,
  ggfxParams = list(),
  withTitle = FALSE,
  args = list()
)
}
\arguments{
\item{seu}{Seurat object}

\item{markers}{results of FindAllMarkers()}

\item{label}{plot label or not}

\item{pt.size}{point size in plot}

\item{reduction}{reduction method}

\item{point_alpha}{point alpha}

\item{use_shadowtext}{use shadowtext for wordcloud}

\item{bg.colour}{shadowtext background color}

\item{which.label}{which label to plot text}

\item{wc_alpha}{alpha value for wordcloud}

\item{wcScale}{scaling value for wordcloud}

\item{rot.per}{ggwordcloud parameter}

\item{rad}{named vector of size of each cluster}

\item{sortBy}{default to avg_log2FC, "log10p" can be specified.}

\item{scaleNumber}{scale the frequency of words by this number
in `gene_name`}

\item{decreasing}{sort by decreasing order or not}

\item{geneNum}{number of genes to be included in wordclouds}

\item{base_ellipse}{if TRUE, wordclouds are placed based on \code{stat_ellipse}.}

\item{base_dens}{if TRUE, wordclouds are placed based on density}

\item{random.order}{ggwordcloud parameter}

\item{gene_name}{show gene names instead of textual information}

\item{withggfx}{applying ggfx filters}

\item{ggfxParams}{parameter list for ggfx}

\item{withTitle}{whether to append title on wordcloud}

\item{args}{parameters to passed to refseq}
}
\value{
single-cell plot with text annotation
}
\description{
Directly output the dimension reduction plot with textual information.
}
\details{
The function accepts the Seurat object and marker gene information
and output the reduced dimension plot with the textual information.
}
\seealso{
plotReducedDimWithTexts
}
