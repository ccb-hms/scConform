#' @title Return the common ancestor of the labels in the prediction set
#' @description Given a prediction set and an ontology represented as a directed
#' graph, this function returns the most specific common ancestor of the labels
#' in the prediction set. It is mainly intended for hierarchical conformal
#' prediction, where a set of predicted labels can be summarized by a single
#' ontology term representing their common ancestor.
#' @param pred_set character vector of labels included in the prediction set.
#' These labels should correspond to node names in `onto`.
#' @param onto an `igraph` object representing the ontology.
#' @return A character string corresponding to the most specific common ancestor
#' of the labels in `pred_set` according to the ontology `onto`.
#' @examples
#' library(igraph)
#' # Let's build a random ontology
#' onto <- graph_from_literal(
#'     animal -+dog:cat, cat -+british:persian,
#'     dog -+cocker:retriever, retriever -+golden:labrador
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
