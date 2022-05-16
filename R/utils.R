#' makeBar
#' 
#' Makeing a barplot of word frequency from queried genes
#' 
#' @examples
#' geneList <- c("6346")
#' makeBar(geneList)
#' @param queries entrez IDs
#' @param top how many numbers of words to be shown
#' @param pal palette used in barplot
#' @param textSize text size in barplot
#' @param reorder order by frequency or not
#' @import org.Hs.eg.db
#' @return barplot of word frequency
#' @export
#' 
makeBar <- function(queries, top=10, pal=NULL, textSize=20, reorder=TRUE) {
    if (is.null(pal)) {
        # palNum <- sample(1:151,1)
        # pal <- pokepal(palNum)
        pal <- palette()
        if (length(pal)<top){
            pal <- rep(pal, ceiling(top/length(pal)))
        }
    }
    wc <- wcGeneSummary(queries,
                        madeUpper=c("dna","rna",
                        tolower(AnnotationDbi::keys(org.Hs.eg.db,
                            keytype="SYMBOL"))))
    barp <- utils::head(wc$df, n=top)
    if (reorder){
        plt <- ggplot(barp, aes(x=reorder(word, freq),
            y=freq, fill=word)) + coord_flip()
    } else {
        plt <- ggplot(barp, aes(x=word, y=freq, fill=word)) +
                coord_flip()
    }     
    plt <- plt +    
        geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
        theme_minimal() + scale_fill_manual(values=pal, guide="none") + 
        theme(axis.text = element_text(size = textSize))
    plt
}

#' exportCyjs
#' 
#' Export Cytoscape.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @import jsonlite
#' @importFrom cyjShiny dataFramesToJSON
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \dontrun{exportCyjs(g, "./", "net")}
#' @export
#' 
exportCyjs <- function(g, rootDir, netDir) {
    
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    
    nodes <- data.frame(
        id=names(V(g)),
        label=names(V(g)),
        image=V(g)$image,
        size=V(g)$size,
        shape=V(g)$shape
    )
    
    edgeList <- as_edgelist(g)
    edges <- data.frame(source=edgeList[,1],
        target=edgeList[,2], interaction=NA)
    edges$strength <- E(g)$strength
    
    pret <- prettify(dataFramesToJSON(edges, nodes))
    pret <- substr(pret, 18, nchar(pret)-3)
    pret
    
    js <- paste0("
    var cy = window.cy = cytoscape({
        container: document.getElementById('cy'),
          style: cytoscape.stylesheet()
            .selector('node')
            .css({
                      'content': 'data(label)',
                      'shape' : 'data(shape)',
                      'background-image': 'data(image)',
                      'text-valign': 'bottom',
                      'background-color': '#FFF',
                      'background-fit': 'cover',
                      'width': 'data(size)',
                      'height': 'data(size)',
                      'font-size' : 'mapData(size, 0, 100, 1, 20)',
                      'text-outline-width': 1,
                      'text-outline-color': '#FFF',
                      'border-color' : '#555',
                      'border-width': 1
                  })
            .selector('edge')
            .css({
                    'width' : '4',
                    'target-arrow-shape': 'triangle',
                    'curve-style': 'bezier',
                    'width' : 'mapData(strength, 0.5, 1, 0, 5)'
                  }),
        'elements':
    ", pret, ",
        layout:{
              name: 'cola',
              padding: 0.5,
              avoidOverlap: true, 
              nodeSpacing: function( node ){ return 0.1; },
              nodeDimensionsIncludeLabels: true
          }
        });
    ")
    
    ## Using cola layout by default.
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdn.jsdelivr.net/npm/cytoscape@3.21.1/dist/cytoscape.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/webcola@3.4.0/WebCola/cola.min.js"></script>
        <script src="https://cdn.jsdelivr.net/npm/cytoscape-cola@2.4.0/cytoscape-cola.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="cy"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- "
    body {
        font-family: helvetica, sans-serif;
        font-size: 14px;
    }
    
    #cy {
        position: absolute;
        left: 0;
        top: 0;
        right: 0;
        bottom: 0;
        z-index: 999;
    }
    
    h1 {
        opacity: 0.5;
        font-size: 1em;
    }"
    
    
    write(js, file = paste0(rootDir, netDir, "/script.js"))
    write(style, file = paste0(rootDir, netDir, "/style.css"))
    write(html, file = paste0(rootDir, netDir, "/index.html"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' exportVisjs
#' 
#' Export vis.js script, HTML and stylesheet for the graph and image
#' 
#' @param g igraph object
#' @param rootDir root directory path
#' @param netDir directory to store scripts
#' @import jsonlite
#' @return return nothing, export to a specified directory
#' @examples
#' library(igraph)
#' g <- graph_from_literal( ME1-+ME2 )
#' V(g)$image <- c("path1","path2")
#' V(g)$shape <- c("image","image")
#' V(g)$size <- c(1,1)
#' \dontrun{exportVisjs(g, "./", "net")}
#' @export
#' 
exportVisjs <- function(g, rootDir, netDir){
    if (is.null(V(g)$shape)){stop("No node shape specified")}
    if (is.null(V(g)$size)){stop("No node size specified")}
    if (is.null(V(g)$image)){stop("No image path specified")}
    if (is.null(E(g)$strength)){E(g)$strength <- rep(1, length(E(g)))}
    if (unique(V(g)$shape)=="rectangle"){
        visjsShape <- "image"
    } else {
        visjsShape <- "circularImage"
    }
    
    nodejson <- toJSON(data.frame(
            id=names(V(g)),
            label=names(V(g)),
            image=V(g)$image,
            shape=visjsShape,
            size=V(g)$size
        ))
    
    edgeList <- as_edgelist(g)
    edgejson <- toJSON(data.frame(from=edgeList[,1], to=edgeList[,2],
        width=E(g)$strength))
    
    # Make JS
    js <- paste0("
    var nodes = null;
    var edges = null;
    var network = null;
    
    function draw() {
    nodes = ", nodejson, ";
    edges = ", edgejson, ";
      var container = document.getElementById('mynetwork');
      var data = {
        nodes: nodes,
        edges: edges,
      };
      var options = {
        nodes: {
          borderWidth: 2,
          size: 30,
          color: {
            border: '#222222',
    background: 'white',
    },
    font: { color: 'black' },
    },
    edges: {
        length: 200,
        color: 'lightgray',
        arrows: { to: {enabled: true} }
    },
    layout: {
        improvedLayout: true
    }
    };
    network = new vis.Network(container, data, options);
    }
    
    window.addEventListener('load', () => {
        draw();
    });
    ")
    
    html <- '
    <!DOCTYPE html>
    <html lang="en">
    
    <head>
        <meta charset="UTF-8">
        <script src="https://cdnjs.cloudflare.com/ajax/libs/vis/4.21.0/vis.min.js"></script>
        <link rel="stylesheet" type="text/css" href="style.css">
    </head>
    
    <body>
        <div id="mynetwork"></div>
        <script src="script.js"></script>
    </body>
    
    </html>
    '
    
    style <- '
    body {
      font: 10pt arial;
    }
    #mynetwork {
      width: 1000px;
      height: 1000px;
      border: 1px solid lightgray;
      background-color: white;
    }
    '
    
    write(js, file = paste0(rootDir, netDir, "/script.js"))
    write(html, file = paste0(rootDir, netDir, "/index.html"))
    write(style, file = paste0(rootDir, netDir, "/style.css"))
    # message(paste0("Exported to ",rootDir,netDir))
}

#' returnExample
#' 
#' return an example dataset used in the analysis
#' 
#' @import org.Hs.eg.db
#' @return return example MEs and colors
#' @examples returnExample()
#' @export
#' 
returnExample <- function() {
    ## Simulate WGCNA results (three modules)
    ccls <- c()
    for (i in c(1,2,3,4,5,6,7,8,9)){
        ccls <- c(ccls, paste0("CCL",i))
    }
    cxcls <- c()
    for (i in c(1,2,3,5,6,8,9,10,11,12,13,14,16)){
        cxcls <- c(cxcls, paste0("CXCL",i))
    }
    erccs <- c("ERCC1","ERCC2","ERCC3","ERCC4","ERCC5","ERCC6","ERCC8")

    CCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=ccls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    CXCLensg <- AnnotationDbi::select(org.Hs.eg.db, keys=cxcls,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL
    ERCCensg <- AnnotationDbi::select(org.Hs.eg.db, keys=erccs,
        columns=c("ENSEMBL"), keytype="SYMBOL")$ENSEMBL

    CCLensg <- CCLensg[!is.na(CCLensg)]
    CXCLensg <- CXCLensg[!is.na(CXCLensg)]
    ERCCensg <- ERCCensg[!is.na(ERCCensg)]

    CCLcol <- rep(1, length(CCLensg))
    names(CCLcol) <- CCLensg
    CXCLcol <- rep(2, length(CXCLensg))
    names(CXCLcol) <- CXCLensg
    ERCCcol <- rep(3, length(ERCCensg))
    names(ERCCcol) <- ERCCensg
    modColors <- c(CCLcol, CXCLcol, ERCCcol)
    ensg <- names(modColors)

    MEs <- data.frame(
        ME1 = 1:10,
        ME2 = 2:11,
        ME3 = 10:1
    )

    mod <- list()
    mod[["MEs"]] <- MEs
    mod[["colors"]] <- modColors
    mod
}