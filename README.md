# biotextgraph

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/noriakis/biotextgraph/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/noriakis/biotextgraph/actions/workflows/R-CMD-check.yaml)
  <!-- badges: end -->

This package makes a summarized visualization of biological entities using textual information from various sources through network approach. The package makes word network using R libraries such as [`GeneSummary`](https://bioconductor.org/packages/release/data/annotation/html/GeneSummary.html), [`tm`](https://www.jstatsoft.org/article/view/v025i05), `ggraph` and `wordcloud`.

### Installation
```R
devtools::install_github("noriakis/biotextgraph")
library(biotextgraph)
```

### Documentation
The documentation using bookdown is available [here](https://noriakis.github.io/software/biotextgraph/), describing the detailed usage.

### Shiny web server
The application using Shiny is available [here](https://nsato.shinyapps.io/biotextgraphweb/), hosted by [shinyapps.io](https://shinyapps.io).


### Examples

<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/wcbn.png?raw=true" width="800px">
</p>
<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/biofabric.png?raw=true" width="800px">
</p>


### Bugs and errors
If you find bugs or errors, please kindly report them to Issues, or make a pull request, or report it directly to [e-mail](nori@hgc.jp).