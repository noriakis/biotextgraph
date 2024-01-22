#' exportWCNetwork
#' 
#' export the wordcloud network provided the gene cluster network
#' Node size will be the number of genes in the cluster
#' By default use cola layout, and see below for the parameters
#' https://github.com/cytoscape/cytoscape.js-cola
#' 
#' If "strength" attribute is in the edge of igraph object,
#' the parameter is used to size the edge width.
#' 
#' @param g igraph (node name corresponding to gene cluster name)
#' @param geneList named list of gene IDs
#' @param dir directory to output the files, default to current directory
#' @param keyType type of gene IDs (default to SYMBOL)
#' @param colors color palettes to use for each cluster
#' named list, or random palette will be assigned if NULL
#' @param orgDb annotation database
#' @param wcArgs parameters for ggwordcloud
#' @param numWords number of words
#' @param sizeMin minimum scale size for nodes
#' @param sizeMax maximum scale size for nodes
#' @param wcScale scaling size for wordcloud
#' @return export the Cytoscape.js network
#' @export
exportWCNetwork <- function(g, geneList, dir="network", colors=NULL,
	sizeMin=50, sizeMax=200, wcScale=20,
	keyType="SYMBOL", orgDb=org.Hs.eg.db, wcArgs=list(), numWords=50) {

	## Check colors
	if (is.null(colors)) {
		colors <- list()
		for (gl in names(geneList)) {
			colors[[gl]] <- brewer.pal(10, sample(row.names(RColorBrewer::brewer.pal.info), 1))
		}
	}

	## Set size
	sizes <- lapply(geneList, function(x) length(x)) |> unlist()
	V(g)$size <- as.numeric(sizes[V(g)$name])

	## Export images
	dir.create(paste0(dir))
	dir.create(paste0(dir, "/images"))

	images <- NULL
	for (gl in names(geneList)) {
		wcArgs[["colors"]] <- colors[[gl]]
		wc <- refseq(geneList[[gl]],
			keyType=keyType,
			orgDb=orgDb,
			plotType="wc",
			argList=wcArgs,
			wcScale=wcScale,
			numWords=numWords)@wc
	    ggsave(paste0(dir, "/images/", gl ,".png"),
            wc, dpi=300, width=10, height=10)
	    images <- c(images, paste0("images/", gl ,".png"))
	}

	V(g)$image <- images
    V(g)$shape <- rep("circle", length(V(g)))

	## Scale the node size
	rawMin <- min(V(g)$size)
	rawMax <- max(V(g)$size)
	scf <- (sizeMax-sizeMin)/(rawMax-rawMin)

	V(g)$size <- scf * V(g)$size + sizeMin - scf * rawMin
	exportCyjs(g, ".", dir)

}