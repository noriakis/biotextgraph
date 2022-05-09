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
                        madeUpper=c("dna","rna",tolower(AnnotationDbi::keys(org.Hs.eg.db, keytype="SYMBOL"))))
    barp <- head(wc$df, n=top)
    if (reorder){
        plt <- ggplot(barp, aes(x=reorder(word, freq), y=freq, fill=word)) + coord_flip()
    } else {
        plt <- ggplot(barp, aes(x=word, y=freq, fill=word)) + coord_flip()
    }     
    plt <- plt +    
        geom_bar(stat = "identity") + xlab("Word") + ylab("Frequency") +
        theme_minimal() + scale_fill_manual(values=pal, guide="none") + 
        theme(axis.text = element_text(size = textSize))
    plt
}