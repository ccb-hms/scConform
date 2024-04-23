#' @title Get parent nodes of a given node
#' @param node the name of a node on the ontology
#' @param onto the considered section of the cell ontology as an igraph object
#' @param include_self if \code{TRUE}, then the function returns also \code{node}
#' @author Daniela Corbetta
#' @return returns a vector with the names of the parents
#' @importFrom igraph V
#' @importFrom igraph degree
#' @importFrom igraph distances

.ancestors <- function(node, onto, include_self=TRUE){
  if(include_self)
    return(V(onto)$name[is.finite(distances(onto, node, mode="in"))])
  else
    return(V(onto)$name[is.finite(distances(onto, node, mode="in")) & V(onto)$name!=node])
}
