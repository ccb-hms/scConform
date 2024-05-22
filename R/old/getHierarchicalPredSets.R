#' @title Get hierarchical prediction sets exploiting the cell ontology
#' @description Let K be the total number of distinct cell type labels and n, m
#' the number of cells in the calibration and in the test data, respectively.
#' This function takes as input two matrices: a matrix \code{n x K} and
#' a matrix \code{m x K} with the estimated
#' probabilities for each cell in the calibration and in the test data, respectively.
#' It returns a list with the prediction sets for each cell in the test data.
#'
#' @param p.cal a nxK matrix of category scores for the cells in the calibration data.
#' Needs to have colnames for each category
#' @param p.test a mxK matrix of category scores for the cells in the test data
#' @param y.cal a vector of length n with the labels of the cells in the calibration data
#' @param onto the considered section of the cell ontology
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @param lambdas a vector of possible lambda values to be considered
#' @author Daniela Corbetta
#' @return The function \code{getHierarchicalPredSets} returns a list of length equal to the number of cells in the test data.
#' Each element of the list contains the prediction set for that cell.
#' @references For reference on conformal risk control, see
#' Angelopoulos, Anastasios N., et al. "Conformal risk control." arXiv preprint arXiv:2208.02814 (2022).
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach


# Function to get the hierarchical prediction sets for the observation in p.test
# Needs a vector of lambda values. For each of the lambdas computes the prediction sets
# for the data in the calibration set (p.cal n_cal x K matrix that contains
# estimated probabilities for each label). Based on these sets that compute the
# loss table and then gets lambda hat based on equation (4) in Bates and
# Angelopoulus (2023), Conformal Risk Control. Finally, builds prediction sets
# for p.test based on the selected lambda value.

.getHierarchicalPredSets <- function(p.cal, p.test, y.cal, onto, alpha, lambdas) {
    y.cal <- as.character(y.cal)
    # Get prediction sets for each value of lambda for all the calibration data
    j <- NULL
    exportedFn <- c(".predSets", ".scores", ".children", ".ancestors")
    sets <- foreach(j = lambdas, .export = exportedFn) %dopar% {
        lapply(
            1:nrow(p.cal),
            function(i) .predSets(lambda = j, pred = p.cal[i, ], onto = onto)
        )
    }

    # Get the loss table (ncal x length(lambda) table with TRUE\FALSE)
    loss <- sapply(1:length(lambdas), function(lambda) {
        sapply(seq_along(y.cal), function(i) {
            !(y.cal[i] %in% sets[[lambda]][[i]])
        })
    })

    # Get lhat
    n <- nrow(loss)
    rhat <- colMeans(loss)
    lhat_idx <- min(which(((n / (n + 1)) * rhat + 1 / (n + 1)) <= alpha))
    lhat <- lambdas[lhat_idx]


    # Get prediction sets for test data
    sets.test <- apply(p.test, 1, function(x) .predSets(lambda = lhat, pred = x, onto = onto))

    return(list(sets.test = sets.test, lhat = lhat))
}


.ancestors <- function(node, onto, include_self = TRUE) {
    if (include_self) {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "in"))])
    } else {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "in")) & V(onto)$name != node])
    }
}

.children <- function(node, onto, leaf = TRUE) {
    if (leaf) {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out")) & degree(onto, mode = "out") == 0])
    } else {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out"))])
    }
}

.scores <- function(pred, int_node, onto) {
    c <- .children(node = int_node, onto = onto, leaf = TRUE)
    return(sum(pred[c]))
}

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
