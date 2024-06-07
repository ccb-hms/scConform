# Function to get the hierarchical prediction sets for the observation in p.test
# Needs a vector of lambda values. For each of the lambdas computes the
# prediction sets for the data in the calibration set
# (p.cal n_cal x K matrix that contains
# estimated probabilities for each label). Based on these sets that compute the
# loss table and then gets lambda hat based on equation (4) in Bates and
# Angelopoulus (2023), Conformal Risk Control. Finally, builds prediction sets
# for p.test based on the selected lambda value.

.getHierarchicalPredSets <- function(
        p.cal, p.test, y.cal, onto, alpha,
        lambdas, BPPARAM) {
    y.cal <- as.character(y.cal)
    # Get prediction sets for each value of lambda for all the calibration data
    j <- NULL

    sets <- bplapply(lambdas, function(j) {
      lapply(seq_len(nrow(p.cal)), function(i) {
        .predSets(lambda = j, pred = p.cal[i, ], onto = onto)
      })
    }, BPPARAM = BPPARAM)

    # Get the loss table (ncal x length(lambda) table with TRUE\FALSE)
    loss <- vapply(seq_along(lambdas), function(lambda) {
        vapply(seq_along(y.cal), function(i) {
            !(y.cal[i] %in% sets[[lambda]][[i]])
        }, logical(1))
    }, FUN.VALUE = logical(length(y.cal)))

    # Get lhat
    n <- nrow(loss)
    rhat <- colMeans(loss)
    lhat_idx <- min(which(((n / (n + 1)) * rhat + 1 / (n + 1)) <= alpha))
    lhat <- lambdas[lhat_idx]


    # Get prediction sets for test data
    sets.test <- apply(p.test, 1, function(x) {
        .predSets(
            lambda = lhat, pred = x,
            onto = onto
        )
    })

    return(sets.test)
}

# Function to retrieve all ancestors of a given node, starting from an igraph
# object. Default includes also the given node
.ancestors <- function(node, onto, include_self = TRUE) {
    if (include_self) {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "in"))])
    } else {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "in")) &
            V(onto)$name != node])
    }
}

# Function to retrieve children of a given node, starting from an igraph object.
# By default, gives only the children that are leaves in the ontology.
.children <- function(node, onto, leaf = TRUE) {
    if (leaf) {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out")) &
            degree(onto, mode = "out") == 0])
    } else {
        return(V(onto)$name[is.finite(distances(onto, node, mode = "out"))])
    }
}

# Function to compute scores of a given node in an ontology (sum of estimated
# probabilities for the leaf nodes that are children of the given one)
.scores <- function(pred, int_node, onto) {
    c <- .children(node = int_node, onto = onto, leaf = TRUE)
    return(sum(pred[c]))
}

# Function to get prediction sets following the ontology
.predSets <- function(lambda, pred, onto) {
    # Get the predicted class and its ancestors
    pred_class <- names(pred)[which.max(pred)]
    anc <- .ancestors(node = pred_class, onto = onto, include_self = TRUE)

    # Compute scores for all the ancestor of the predicted class
    s <- vapply(as.character(anc), function(i) {
        .scores(
            pred = pred, int_node = i,
            onto = onto
        )
    }, numeric(1))
    names(s) <- anc

    ## Sort them by score and if there are ties by distance to the predicted
    ## class
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

