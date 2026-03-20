#' @title Plot prediction sets
#' @description This function takes as input a prediction set and an
#' ontology and plots the ontology, highlighting the labels included in the set.
#'
#' @param pred_set character vector containing the labels in the prediction set
#' @param onto an `igraph` object representing the ontology
#' @param probs numeric vector of estimated probabilities for the classes. The
#'   names of `probs` should correspond to node names in `onto`
#' @param col_grad character vector of colors used to highlight the classes. If
#'   `probs` is provided, this should contain a color gradient; otherwise,
#'   a single color can be supplied
#' @param attrs attrs list of additional graphical attributes passed to
#' `plot()`
#' @param k integer number of decimal digits to consider in `probs`
#' @param title title of the plot
#' @param add_scores Logical. If `TRUE`, estimated probabilities are
#' added to the name of the classes
#' @param ... additional graphical parameters passed to `plot()`
#' @return A plot of the ontology with the labels in the prediction set
#'   highlighted.
#' @examples
#' library(igraph)
#' # Let's build a random ontology
#' onto <- graph_from_literal(
#'     animal -+dog:cat, cat -+british:persian,
#'     dog -+cocker:retriever, retriever -+golden:labrador
#' )
#' # Let's consider this prediction set
#' pred_set <- c("golden", "labrador", "cocker")
#' plotResult(pred_set, onto,
#'     col_grad = "pink", add_scores = FALSE,
#'     title = "Prediction set"
#' )
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom igraph as_graphnel V
#' @importFrom Rgraphviz plot
#' @export
#'

plotResult <- function(
    pred_set,
    onto,
    probs = NULL,
    col_grad = c("lemonchiffon", "orange", "darkred"),
    attrs = NULL,
    k = 4,
    title = NULL,
    add_scores = TRUE,
    ...) {
  if (add_scores && is.null(probs)) {
    stop("'probs' must be provided when 'add_scores = TRUE'.")
  }

  graph <- as_graphnel(onto)

  if (!is.null(probs)) {
    p <- round(probs, k) * 10^k
    colfunc <- colorRampPalette(col_grad)
    palette_cols <- colfunc(10^k)

    vec_col <- vapply(pred_set, function(i) {
      if (p[i] == 0) {
        palette_cols[1]
      } else {
        palette_cols[p[i]]
      }
    }, character(1))
  } else {
      vec_col <- rep(col_grad[1], length(pred_set))
  }

  names(vec_col) <- pred_set
  nAttrs <- list(fillcolor = vec_col)

  if (add_scores) {
    node_names <- V(onto)$name
    scores <- round(vapply(
      node_names,
      function(x) .scores(probs, x, onto),
      numeric(1)
    ), 3)

    labels <- paste(node_names, scores, sep = ", ")
    names(labels) <- node_names
    nAttrs$label <- labels
  }

  plot(graph, attrs = attrs, nodeAttrs = nAttrs, main = title, ...)
}
