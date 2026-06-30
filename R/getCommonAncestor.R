#' @title Return the common ancestor of the labels in the prediction set
#' @description Given a prediction set and an ontology represented as a directed
#' graph, this function returns the most specific common ancestor of the labels
#' in the prediction set. It is mainly intended for hierarchical conformal
#' prediction, where a set of predicted labels can be summarized by a single
#' ontology term representing their common ancestor.
#' @param pred_set character vector of labels included in the prediction set.
#' These labels should correspond to node names in `onto`.
#' @param onto an `igraph` object representing the ontology.
#' @param onto_cache optional precomputed cache from `precomputeOnto()`. If `NULL`, 
#' the cache is computed automatically (and discarded after the call).
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
getCommonAncestor <- function(pred_set, onto, onto_cache = NULL) {
    # Check if the prediction set is empty
    if (length(pred_set) == 0) {
        stop("Prediction set is empty. Cannot compute common ancestor.")
    }
    if(is.null(onto_cache)) {
        message("Building ontology cache...")
        onto_cache <- precomputeOnto(onto)
    }
    # Check if all labels in the prediction set are present in the ontology
    if (!all(pred_set %in% onto_cache$node_names)) {
        missing_labels <- pred_set[!pred_set %in% onto_cache$node_names]
        stop("The following labels in the prediction set are not present in the ontology: ", 
             paste(missing_labels, collapse = ", "))
    }
    com_anc <- Reduce(intersect, lapply(pred_set, function(node) {
        onto_cache$ancestors_of[[node]]
    }))
    root <- onto_cache$root
    # Most specific = farthest from root = largest distance root -> ancestor
    # dist_mat[root, com_anc]: distance from root to each common ancestor
    dists <- onto_cache$dist_mat[root, com_anc]
    first_anc <- com_anc[which.max(dists)]
    return(first_anc)
}
