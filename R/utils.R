###############################################################################
### 1. Utils for standard conformal inference
###############################################################################

# Function to get conformal prediction sets with split conformal inference

.getConformalPredSets <- function(p_cal, p_test, y_cal, alpha) {
    # Get calibration scores (1-predicted probability for the true class)
    true <- rep(NA, dim(p_cal)[1])
    for (i in seq_len(dim(p_cal)[1])) {
        true[i] <- p_cal[i, y_cal[i]]
    }
    s <- 1 - true

    # Get adjusted quantile
    n <- nrow(p_cal)
    q_level <- ceiling((n + 1) * (1 - alpha)) / n
    qhat <- quantile(s, q_level)

    # Get prediction sets
    prediction_sets <- p_test >= 1 - qhat
    pr_list <- lapply(seq_len(nrow(prediction_sets)), function(i) {
        colnames(prediction_sets)[prediction_sets[i, ]]
    })
    return(pr_list)
}


###############################################################################
### 2. Utils for hierarchical prediction set
###############################################################################

# Function to get the hierarchical prediction sets for the observations in
# p_test.
# Needs a vector of lambda values. For each of the lambdas computes the
# prediction sets for the data in the calibration set
# (p_cal n_cal x K matrix that contains
# estimated probabilities for each label). Based on these sets that compute the
# loss table and then gets lambda hat based on equation (4) in Bates and
# Angelopoulus (2023), Conformal Risk Control. Finally, builds prediction sets
# for p_test based on the selected lambda value.

.getHierarchicalPredSets <- function(
        p_cal, p_test, y_cal, onto, alpha,
        lambdas, BPPARAM) {
    y_cal <- as.character(y_cal)
    # Get prediction sets for each value of lambda for all the calibration data
    j <- NULL

    sets <- bplapply(lambdas, function(j) {
        lapply(seq_len(nrow(p_cal)), function(i) {
            .predSets(lambda = j, pred = p_cal[i, ], onto = onto)
        })
    }, BPPARAM = BPPARAM)

    # Get the loss table (ncal x length(lambda) table with TRUE\FALSE)
    loss <- vapply(seq_along(lambdas), function(lambda) {
        vapply(seq_along(y_cal), function(i) {
            !(y_cal[i] %in% sets[[lambda]][[i]])
        }, logical(1))
    }, FUN.VALUE = logical(length(y_cal)))

    # Get lhat
    n <- nrow(loss)
    rhat <- colMeans(loss)
    lhat_idx <- min(which(((n / (n + 1)) * rhat + 1 / (n + 1)) <= alpha))
    lhat <- lambdas[lhat_idx]


    # Get prediction sets for test data
    sets_test <- apply(p_test, 1, function(x) {
        .predSets(
            lambda = lhat, pred = x,
            onto = onto
        )
    })

    return(sets_test)
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

# Function to retrieve children of a given node, starting from an igraph object
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

###############################################################################
### 3. Utils for resampling
###############################################################################

# Function to implement the resampling strategy when calibration and test
# set supposedly have a different distribution of the cell types.
# Right now it implements a two-fold strategy, dividing randomly
# test data in two.

.resampleTwo <- function(p_cal, p_test, y_cal, labels) {
    s <- sample(seq_len(nrow(p_test)), round(nrow(p_test) / 2))
    test1 <- p_test[s, ]
    test2 <- p_test[-s, ]

    # Compute predicted class
    pr_class1 <- apply(test1, 1, function(row) colnames(test1)[which.max(row)])
    pr_class2 <- apply(test2, 1, function(row) colnames(test2)[which.max(row)])
    test_freq1 <- prop.table(table(pr_class1))
    test_freq2 <- prop.table(table(pr_class2))
    # Transform to absolute frequencies
    des_freq1 <- round(test_freq1 * length(y_cal))
    des_freq2 <- round(test_freq2 * length(y_cal))

    idx1 <- idx2 <- NULL
    for (i in labels) {
        category <- which(y_cal == i)
        if (!is.na(des_freq1[i])) {
            idx_category1 <- sample(category,
                size = des_freq1[i],
                replace = TRUE
            )
            idx1 <- c(idx1, idx_category1)
        }
        if (!is.na(des_freq2[i])) {
            idx_category2 <- sample(category,
                size = des_freq2[i],
                replace = TRUE
            )
            idx2 <- c(idx2, idx_category2)
        }
    }

    return(
        list(
            p_cal1 = p_cal[idx1, ],
            p_cal2 = p_cal[idx2, ],
            p_test1 = test1,
            p_test2 = test2,
            y_cal1 = y_cal[idx1],
            y_cal2 = y_cal[idx2],
            idx = c(s, setdiff(seq_len(nrow(p_test)), s))
        ) # index in the original data
    )
}

###############################################################################
### 4. General
###############################################################################

## function to retrieve prediction matrix from the colData of a
## SingleCellExperiment object

.retrievePredMatrix <- function(sc, K, labels) {
    n_sc <- ncol(sc)
    p_sc <- matrix(NA, nrow = n_sc, ncol = K)
    colnames(p_sc) <- labels
    for (i in labels) {
        p_sc[, i] <- colData(sc)[[i]]
    }
    return(p_sc)
}
