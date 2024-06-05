#' @title Plot prediction sets
#' @description This function takes as input a prediction set and an
#' ontology and plots the ontology, highlighting the labels included in the set.
#'
#' @param pred.set character vector containing the labels in the prediction set
#' @param onto ontology as an igraph
#' @param probs estimated probabilities for the classes
#' @param col.grad color to use to highlight the classes
#' @param attrs other grtaphical attributes
#' @param k number of decimal digits to consider in \code{probs}
#' @param title title of the plot
#' @param add.scores boolean. If \code{TRUE}, estimated probabilities are
#' added to the name of the classes
#' @return a plot of the ontology with the considered classes colored
#' @importFrom grDevices colorRampPalette
#' @importFrom igraph as_graphnel
#' @import Rgraphviz
#' @export
#'

plotResult <- function(pred.set, onto, probs=NULL,
                       col.grad=c("lemonchiffon", "orange", "darkred"),
                       attrs=NULL, k=4, title=NULL,
                       add.scores=TRUE,...){
  ## Function works with graphnel
  graph <- as_graphnel(onto)
  vec.col <- NULL
  if(!is.null(probs)){
    p <- round(probs, k)*10^k
    colfunc <- colorRampPalette(col.grad)
    for(i in pred.set)
      vec.col[i] <- ifelse(p[i]==0, colfunc(10^k)[1], colfunc(10^k)[p[i]])
  }
  else {
    for(i in pred.set)
      vec.col[i] <- col.grad
  }
  nAttrs <- list()
  nAttrs$fillcolor <- vec.col

  # Add scores if requested
  if(add.scores){
    # Need igraph for .scores function
    scores <- round(sapply(V(onto)$name, function(x) .scores(probs, x, onto)),3)
    labels <- NULL
    for(i in 1:length(V(onto)$name))
      labels[i] <- paste(V(onto)$name[i], scores[i], sep = ", ")
    names(labels) <- names(scores)
    nAttrs$label <- labels
  }
  plot(graph, attrs=attrs, nodeAttrs=nAttrs, main=title)
}



