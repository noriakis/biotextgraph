% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps.R
\name{obtain_enrich}
\alias{obtain_enrich}
\title{obtain_enrich}
\usage{
obtain_enrich(
  geneList,
  keyType = "SYMBOL",
  enrich = "kegg",
  org_db = org.Hs.eg.db,
  top_path = 30
)
}
\arguments{
\item{geneList}{gene list}

\item{keyType}{key type of gene list}

\item{enrich}{`kegg` or `reactome`}

\item{org_db}{used to convert ID}

\item{top_path}{the number of pathways to be obtained
sorted by p-values}
}
\value{
biotext class object
}
\description{
Obtain enrichment analysis description
}
\details{
The function performs enrichment analysis on
queried vector, and returns the biotext object using the
textual information from biological pathways.
}
\examples{
if (requireNamespace("clusterProfiler")) {
    testgenes <- c("IRF3","PNKP","DDX41","ERCC1","ERCC2","XRCC1")
    obtain_enrich(testgenes)
}
}
