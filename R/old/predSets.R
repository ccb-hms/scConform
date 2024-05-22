#' @title Get hierarchical prediction sets with given lambda for a single line
#' @param lambda score/probability treshold
#' @param pred a named vectors with the estimated probabilities for each cell type. Names
#' have to correspond to names of the leaf nodes in the ontology
#' @param onto the considered section of the cell ontology as an igraph object
#' @author Daniela Corbetta
#' @return vector with names of the selected leaf nodes

.predSets <- function(lambda, pred, onto) {
    # Get the predicted class and its ancestors
    pred_class <- names(pred)[which.max(pred)]
    anc <- .ancestors(node = pred_class, onto = onto, include_self = TRUE)

    # Compute scores for all the ancestor of the predicted class
    s <- sapply(as.character(anc), function(i) {
        .scores(
            pred = pred, int_node = i,
            onto = onto
        )
    })
    names(s) <- anc

    # Sort them by score and if there are ties by distance to the predicted class
    ## compute distance from predicted class
    pos <- distances(onto, v = anc, to = pred_class, mode = "out")
    tie_breaker <- as.vector(t(pos))
    names(tie_breaker) <- colnames(t(pos))
    sorted_indices <- order(s, tie_breaker, decreasing = FALSE)
    sorted_scores <- s[sorted_indices]

    # Select the first score that is geq than lambda
    sel_node <- names(sorted_scores)[round(sorted_scores, 15) >= lambda][1]


    # Add also the subgraphs we would have obtained with smaller lambda
    selected <- c(
        lapply(anc[round(s, 15) < lambda], function(x) {
            .children(node = x, onto = onto)
        }),
        list(.children(sel_node, onto))
    )

    return(Reduce(union, selected))
}
