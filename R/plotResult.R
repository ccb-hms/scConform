#' @title Plot prediction sets
#' @description This function takes as input a prediction set and an
#' ontology and plots the ontology, highlighting the labels included in the set.
#'
#' @param pred_set character vector containing the labels in the prediction set
#' @param onto ontology as an igraph
#' @param probs estimated probabilities for the classes
#' @param col_grad color to use to highlight the classes
#' @param attrs other grtaphical attributes
#' @param k number of decimal digits to consider in \code{probs}
#' @param title title of the plot
#' @param add_scores boolean. If \code{TRUE}, estimated probabilities are
#' added to the name of the classes
#' @param ... general commands to be sent to plot
#' @return a plot of the ontology with the considered classes colored
#' @examples
#' library(igraph)
#' # Let's build a random ontology
#' onto <- graph_from_literal(
#'     animal-+dog:cat, cat-+british:persian,
#'     dog-+cocker:retriever, retriever-+golden:labrador
#' )
#' # Let's consider this prediction set
#' pred_set <- c("golden", "labrador", "cocker")
#' plotResult(pred_set, onto,
#'     col_grad = "pink", add_scores = FALSE,
#'     title = "Prediction set"
#' )
#'
#' @importFrom grDevices colorRampPalette
#' @importFrom igraph as_graphnel
#' @import Rgraphviz
#' @export
#'

plotResult <- function(
        pred_set, onto, probs = NULL,
        col_grad = c("lemonchiffon", "orange", "darkred"),
        attrs = NULL, k = 4, title = NULL,
        add_scores = TRUE, ...) {
    ## Function works with graphnel
    graph <- as_graphnel(onto)
    vec_col <- NULL
    if (!is.null(probs)) {
        p <- round(probs, k) * 10^k
        colfunc <- colorRampPalette(col_grad)
        for (i in pred_set) {
            vec_col[i] <- ifelse(
                p[i] == 0,
                colfunc(10^k)[1],
                colfunc(10^k)[p[i]]
            )
        }
    } else {
        for (i in pred_set) {
            vec_col[i] <- col_grad
        }
    }
    nAttrs <- list()
    nAttrs$fillcolor <- vec_col

    # Add scores if requested
    if (add_scores) {
        # Need igraph for .scores function
        scores <- round(vapply(
            V(onto)$name,
            function(x) .scores(probs, x, onto),
            numeric(1)
        ), 3)
        labels <- NULL
        for (i in seq_along(V(onto)$name)) {
            labels[i] <- paste(V(onto)$name[i], scores[i], sep = ", ")
        }
        names(labels) <- names(scores)
        nAttrs$label <- labels
    }
    plot(graph, attrs = attrs, nodeAttrs = nAttrs, main = title)
}
