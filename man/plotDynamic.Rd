% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compareWordNet.R
\name{plotDynamic}
\alias{plotDynamic}
\title{plotDynamic}
\usage{
plotDynamic(
  listOfNets,
  concat = "union",
  alpha = 0.5,
  titles = NULL,
  tag = FALSE,
  useDynamic = TRUE
)
}
\arguments{
\item{listOfNets}{list consisting results of wc* functions (plotType="network")}

\item{concat}{"union" or "intersection"}

\item{alpha}{pass to layout_as_dynamic}

\item{titles}{title to be shown on plot}

\item{tag}{show tag on plot}

\item{useDynamic}{use layout_as_dynamic}
}
\value{
plot comparing gene clusters
}
\description{
list network of words using graphlayouts::layout_as_dynamic
}
\examples{
library(igraph)
wc1 <- refseq(c("DDX41","IRF3"), plotType="network")
wc2 <- refseq(c("DDX41","PNKP"), plotType="network")
compare <- plotDynamic(list(wc1, wc2))
}
