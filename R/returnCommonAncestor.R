#' @title Return the common ancestor of the labels in the prediction set
#' @description This function takes as input a prediction set and an
#' ontology and returns the common ancestor of the labels in the prediction set.
#' It is useful when using hierarchical
#'
#' @param pred.set character vector containing the labels in the prediction set
#' @param onto ontology as an igraph object
#' @return the common ancestor of the labels in \code{pred.set}, according to
#' the ontology \code{onto}
#' @examples
#' # Let's build a random ontology
#' onto <- graph_from_literal(animal-+dog:cat, cat-+british:persian,
#' dog-+cocker:retriever, retriever-+golden:labrador)
#' # Let's consider this prediction set
#' pred.set <- c("golden", "labrador", "cocker")
#' com.anc <- returnCommonAncestor(pred.set, onto)
#' @importFrom igraph V distances degree
#' @export

# Function to return the common ancestor instead of the single leaf nodes
returnCommonAncestor <- function(pred.set, onto){
    com.anc <- Reduce(intersect, lapply(pred.set, function(node) {
        .ancestors(node, onto)
    }))
    root <- V(onto)$name[degree(onto, mode = "in") == 0]
    first.anc <- com.anc[which.max(distances(onto, v = com.anc, to = root))]
    return(first.anc)
}
