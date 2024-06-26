% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sc.R
\name{TextMarkers}
\alias{TextMarkers}
\title{TextMarkers}
\usage{
TextMarkers(
  df,
  keyType = "SYMBOL",
  type = "wc",
  genePlot = TRUE,
  genePlotNum = 5,
  colorText = TRUE,
  args = list(),
  wcArgs = NULL,
  raw = FALSE,
  col = NULL,
  pvalThresh = 0.05,
  withTitle = TRUE,
  withggfx = NULL,
  ggfxParams = list()
)
}
\arguments{
\item{df}{result of FindAllMarkers()}

\item{keyType}{keytype}

\item{type}{wc or network}

\item{genePlot}{whether to plot relevant genes}

\item{genePlotNum}{number of genes to plot}

\item{colorText}{colorlize text or not, default to TRUE}

\item{args}{parameters to pass to refseq}

\item{wcArgs}{parameters to pass to ggwordcloud}

\item{raw}{obtain raw results of refseq instead of plot}

\item{col}{color to be used in wordcloud,
named list of clusters}

\item{pvalThresh}{adjusted p-value threshold for markers}

\item{withTitle}{plot title on the plot}

\item{withggfx}{applying ggfx filters}

\item{ggfxParams}{parameter list for ggfx}
}
\value{
list of plots on textual information in the gene cluster
}
\description{
Take results of `Seurat::FindAllMarkers` as input and plot wordcloud or network for all the clusters.
}
\details{
Using the results of marker gene identification such as `FindAllMarkers` from Seurat,
Recursively summarize the textual information of markers and output the visualizations.
}
\examples{
markers <- data.frame(p_val_adj=c(0.01, 0.01, 0.04),
gene=c("PNKP","DDX41","IRF3"),cluster=c("1","1","1"))
colors <- list("1"="red")
TextMarkers(markers, col=colors, type="wc")
}
\seealso{
obtainMarkersWC
}
