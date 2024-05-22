#' @title Get children nodes of a given node
#' @param node the name of a node on the ontology
#' @param onto the considered section of the cell ontology as an igraph object
#' @param leaf if \code{TRUE}, then the function returns only the leaf nodes that are children of the given one.
#' Otherwise it returns all the children.
#' @author Daniela Corbetta
#' @return returns a vector with the names of the children.
#' @importFrom igraph V
#' @importFrom igraph degree
#' @importFrom igraph distances


.children <- function(node, onto, leaf = TRUE) {
    if (leaf) {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out")) & degree(onto, mode = "out") == 0])
    } else {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out"))])
    }
}
