#' @title Get scores associated with a node
#' @description The score is given by the sum of the estimated probabilities of the leaf nodes that
#' are children of the given node
#' @param pred a named vectors with the estimated probabilities for each cell type. Names
#' have to correspond to names of the leaf nodes in the ontology
#' @param int_node the name of a node on the ontology
#' @param onto the considered section of the cell ontology as an igraph object
#' @author Daniela Corbetta
#' @return returns the score

.scores <- function(pred, int_node, onto){
    c <- .children(node = int_node, onto = onto, leaf = TRUE)
    return(sum(pred[c]))
}
