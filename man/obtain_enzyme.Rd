% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steps.R
\name{obtain_enzyme}
\alias{obtain_enzyme}
\title{obtain_enzyme}
\usage{
obtain_enzyme(
  file,
  ec_num,
  only_term = FALSE,
  only_df = FALSE,
  target = "abstract",
  tax_ec = FALSE,
  tax_file = NULL,
  cand_tax = NULL,
  arg_list = list(),
  api_key = NULL
)
}
\arguments{
\item{file}{file downloaded from expasy}

\item{ec_num}{candidate ecnum, like those obtained from eggNOG-mapper}

\item{only_term}{only return quoted queries to pubmed}

\item{only_df}{only return ec description data.frame
if onlyTerm and onlyDf are both specified, onlyTerm have priority}

\item{target}{abstract or title}

\item{tax_ec}{link taxonomy to EC using UniProt Taxonomy ID file
If this is TRUE, data.frame is returned}

\item{tax_file}{UniProt organism ID file path}

\item{cand_tax}{when taxec=TRUE, search only for these species.}

\item{arg_list}{passed to obtain_pubmed()}

\item{api_key}{api key for PubMed}
}
\value{
biotext class object
}
\description{
Obtain EC-related text data from PubMed
}
\details{
Query the Enzyme Comission number and obtain description,
and search pubmed for these enzymes and make word cloud and
correlation network. Need to specify the path to "enzyme.dat"
downloaded from from ExPASy (https://enzyme.expasy.org/).
}
\seealso{
generalf
}
