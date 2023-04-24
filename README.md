# biotextgraph

  <!-- badges: start -->
  [![R-CMD-check](https://github.com/noriakis/biotextgraph/workflows/R-CMD-check/badge.svg)](https://github.com/noriakis/biotextgraph/actions)
  <!-- badges: end -->

This package makes a summarized visualization of biological entities using textual information from various sources through network approach. The package makes word network using R libraries such as [`GeneSummary`](https://bioconductor.org/packages/release/data/annotation/html/GeneSummary.html), [`tm`](https://www.jstatsoft.org/article/view/v025i05), `ggraph` and `wordcloud`.

### Installation
```R
devtools::install_github("noriakis/biotextgraph")
library(biotextgraph)
```

### [Documentation](https://noriakis.github.io/software/biotextgraph/)

### [Shiny web server](https://nsato.shinyapps.io/biotextgraphweb/)

<p align="center">
<img src="https://github.com/noriakis/software/blob/main/images/wcbn.png?raw=true" width="800px">
</p>