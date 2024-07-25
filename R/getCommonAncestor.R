#' @title Return the common ancestor of the labels in the prediction set
#' @description This function takes as input a prediction set and an
#' ontology and returns the common ancestor of the labels in the prediction set.
#' It is useful when using hierarchical
#'
#' @param pred_set character vector containing the labels in the prediction set
#' @param onto ontology as an igraph object
#' @return the common ancestor of the labels in \code{pred_set}, according to
#' the ontology \code{onto}
#' @examples
#' library(igraph)
#' # Let's build a random ontology
#' onto <- graph_from_literal(
#'     animal-+dog:cat, cat-+british:persian,
#'     dog-+cocker:retriever, retriever-+golden:labrador
#' )
#' # Let's consider this prediction set
#' pred_set <- c("golden", "labrador", "cocker")
#' com_anc <- getCommonAncestor(pred_set, onto)
#' @importFrom igraph V distances degree
#' @export

# Function to return the common ancestor instead of the single leaf nodes
getCommonAncestor <- function(pred_set, onto) {
    com_anc <- Reduce(intersect, lapply(pred_set, function(node) {
        .ancestors(node, onto)
    }))
    root <- V(onto)$name[degree(onto, mode = "in") == 0]
    first_anc <- com_anc[which.max(distances(onto, v = com_anc, to = root))]
    return(first_anc)
}
