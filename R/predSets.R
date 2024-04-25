#' @title Get hierarchical prediction sets with given lambda for a single line
#' @param lambda score/probability treshold
#' @param pred a named vectors with the estimated probabilities for each cell type. Names
#' have to correspond to names of the leaf nodes in the ontology
#' @param onto the considered section of the cell ontology as an igraph object
#' @author Daniela Corbetta
#' @return vector with names of the selected leaf nodes

.predSets <- function(lambda, pred, onto){
    # Get the predicted class and its ancestors
    pred_class <- names(pred)[which.max(pred)]
    anc <- .ancestors(node = pred_class, onto = onto, include_self = TRUE)

    # Compute scores for all the ancestor of the predicted class
    scores <- sapply(as.character(anc), function(i) scores(pred=pred, int_node=i, onto=onto))
    names(scores) <- anc

    # Select nodes with scores higher than lambda
    sel_scores <- scores[round(scores, 15) >= lambda]
    # Select the most external node (i.e. the one with smallest score)
    sel_node <- names(sel_scores)[length(sel_scores)]
    # Add also the subgraphs we would have obtained with smaller lambda
    selected <- c(lapply(anc[round(scores, 15) <= lambda], function(x) .children(node = x, onto = onto)),
                  list(.children(sel_node, onto)))

    return(Reduce(union, selected))
}
