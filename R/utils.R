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

.getHierarchicalPredSets <- function(p_cal, p_test, y_cal, onto, onto_cache, alpha,
                                     lambdas, BPPARAM,
                                     method = "full") {
    # method <- match.arg(method)
    y_cal <- as.character(y_cal)

    # Select sets construction algorithm
    if (is.character(method)) {
        pred_fun <- switch(method,
            full   = .predSets,
            step   = .predSetsStep,
            rank   = .predSetsRank
        )
    } else if (is.function(method)) {
        pred_fun <- method
    } else {
        stop("Invalid 'method' argument")
    }

    # Precompute per-cell scores for all calibration cells:
    # For each calibration cell i and each of its ancestors a,
    # score g(a, x_i) = sum of pred probs for leaf children of a.
 
    # For calibration: for each cell, we only need ancestors of the pred class.
    # Shape: list of length n_cal, each element a named numeric vector of scores
    # indexed by ancestor node name.
    # Validate that at least some ontology leaves are present in the prediction
    # matrix columns. For method = "rank" this is required; for others it means
    # all conformity scores will be NA and calibration will fail anyway.
    common_leaves <- intersect(onto_cache$leaves, colnames(p_cal))
    if (length(common_leaves) == 0L) {
        stop("No ontology leaves are present in the column names of the ",
             "prediction matrix. Check that the ontology and the prediction ",
             "matrix share the same labels.")
    }

    n_cal  <- nrow(p_cal)
    n_test <- nrow(p_test)

    cal_scores <- lapply(seq_len(n_cal), function(i) {
        .precomputeCellScores(p_cal[i, ], onto_cache)
    })

    test_scores <- lapply(seq_len(n_test), function(i) {
        .precomputeCellScores(p_test[i, ], onto_cache)
    })

    j <- NULL
    sets <- bplapply(lambdas, function(j) {
        lapply(seq_len(n_cal), function(i) {
            args <- list(lambda = j, pred = p_cal[i, ], onto = onto)
            if ("onto_cache" %in% names(formals(pred_fun))) args$onto_cache <- onto_cache
            if ("cell_cache" %in% names(formals(pred_fun))) args$cell_cache <- cal_scores[[i]]
            do.call(pred_fun, args)
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
    message("Calibration complete. Selected lambda_hat = ", lhat)

    sets_test <- lapply(seq_len(n_test), function(i) {
        args <- list(lambda = lhat, pred = p_test[i, ], onto = onto)
        if ("onto_cache" %in% names(formals(pred_fun))) args$onto_cache <- onto_cache
        if ("cell_cache" %in% names(formals(pred_fun))) args$cell_cache <- test_scores[[i]]
        do.call(pred_fun, args)
    })

    return(sets_test)
}

# Function to compute scores of a given node in an ontology (sum of estimated
# probabilities for the leaf nodes that are children of the given one)
.scores <- function(pred, int_node, onto_cache) {
    lc <- onto_cache$leaf_children[[int_node]]
    # Only leaves present in pred vector
    lc_present <- intersect(lc, names(pred))
    if (length(lc_present) == 0L) {
        return(NA_real_)
    }
    return(sum(pred[lc_present]))
}

# Precompute per-cell ancestor scores.
# Returns a named list with:
#   $pred_class : name of the predicted class
#   $anc        : character vector of ancestor names (incl. self)
#   $s          : named numeric vector of scores for each ancestor
#   $tie_breaker: named numeric vector of distances to pred_class
.precomputeCellScores <- function(pred, onto_cache) {
    pred_class <- names(pred)[which.max(pred)]

    # pred_class may not be an ontology node (labels mismatch).
    # Return an empty cache entry; pred functions will produce empty sets.
    anc <- onto_cache$ancestors_of[[pred_class]]
    if (is.null(anc)) {
        return(list(
            pred_class  = pred_class,
            anc         = character(0),
            s           = numeric(0),
            tie_breaker = numeric(0)
        ))
    }

    # Score for each ancestor: sum of pred probs over its leaf children
    s <- vapply(anc, function(a) {
        .scores(pred, a, onto_cache)
    }, numeric(1))
    names(s) <- anc

    # Distance from each ancestor to pred_class (for tie-breaking)
    # dist_mat[ancestor, pred_class]: finite if ancestor can reach pred_class
    tie_breaker <- onto_cache$dist_mat[anc, pred_class]
    names(tie_breaker) <- anc

    list(
        pred_class   = pred_class,
        anc          = anc,
        s            = s,
        tie_breaker  = tie_breaker
    )
}

################################
### Prediction sets construction
################################

# Function to get prediction sets following the ontology (full version)
.predSets <- function(lambda, pred, onto, onto_cache = NULL, cell_cache = NULL) {
    pred_class  <- cell_cache$pred_class
    anc         <- cell_cache$anc
    s           <- cell_cache$s
    tie_breaker <- cell_cache$tie_breaker
    sorted_indices <- order(s, tie_breaker, decreasing = FALSE)
    sorted_scores  <- s[sorted_indices]

    sel_node <- names(sorted_scores)[round(sorted_scores, 15) >= lambda][1]

    selected <- c(
        lapply(anc[round(s, 15) < lambda], function(x) {
            onto_cache$leaf_children[[x]]
        }),
        list(onto_cache$leaf_children[[sel_node]])
    )
    return(Reduce(union, selected))
}

# Reduced (simple) hierarchical prediction sets (reduced version, no union
# with L(v))
.predSetsSimple <- function(lambda, pred, onto, onto_cache = NULL, cell_cache = NULL) {
    anc <- cell_cache$anc
    s   <- cell_cache$s
    
    # Thresholded ancestors (exclude NA scores)
    selected <- lapply(
        anc[!is.na(s) & round(s, 10) <= lambda],
        function(x) onto_cache$leaf_children[[x]]
    )
    
    return(Reduce(union, selected))
}

# Nested sets by number of steps up the DAG (k = lambda)
.predSetsStep <- function(lambda, pred, onto, onto_cache = NULL, cell_cache = NULL) {
    k <- as.integer(lambda)
    if (is.na(k) || k < 0){
        stop("For method = 'step', lambda must be a non-negative integer.")
    }

    pred_class  <- cell_cache$pred_class
    anc         <- cell_cache$anc
    tie_breaker <- cell_cache$tie_breaker   # dist from anc to pred_class
    anc_keep    <- anc[is.finite(tie_breaker) & tie_breaker <= k]
    # Union of leaves of kept ancestors
    selected    <- lapply(anc_keep, function(a) onto_cache$leaf_children[[a]])
    if (length(selected) == 0) {
        return(character(0))
    }
    return(Reduce(union, selected))
}

# Prediction sets by probability ranking + LCA
.predSetsRank <- function(lambda, pred, onto, onto_cache = NULL, cell_cache = NULL) {
    # lambda is a cumulative probability threshold in [0, 1]
    if (is.na(lambda) || lambda < 0 || lambda > 1) {
        stop("For method = 'rank', lambda must be in [0, 1].")
    }

    # Ontology leaves
    leaves <- onto_cache$leaves

    # Keep only leaves present in prediction vector
    leaves <- intersect(leaves, names(pred))
    if (length(leaves) == 0) {
        stop("No ontology leaves are present among prediction vector names.")
    }

    # Probabilities for leaves
    p <- pred[leaves]
    p[is.na(p)] <- 0

    # Rank leaves by decreasing probability
    ord <- order(p, decreasing = TRUE)
    leaves_ord <- leaves[ord]
    p_ord <- p[ord]

    # Select leaves until cumulative probability exceeds lambda
    cumprob <- cumsum(p_ord)
    m <- which(cumprob >= lambda)[1]

    if (is.na(m)) {
        m <- length(leaves_ord)
    }
    if (m < 1) {
        m <- 1
    }

    selected_leaves <- leaves_ord[seq_len(m)]

    # Lowest common ancestor of selected leaves
    lca <- getCommonAncestor(selected_leaves, onto, onto_cache = onto_cache)

    # Return all leaves under the LCA
    return(onto_cache$leaf_children[[lca]])
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
    test1 <- p_test[s, , drop = FALSE]
    test2 <- p_test[-s, , drop = FALSE]

    # Compute predicted class
    pr_class1 <- apply(test1, 1, function(row) colnames(test1)[which.max(row)])
    pr_class2 <- apply(test2, 1, function(row) colnames(test2)[which.max(row)])
    test_freq1 <- prop.table(table(pr_class1))
    test_freq2 <- prop.table(table(pr_class2))

    # Transform to absolute frequencies
    des_freq1 <- round(test_freq1 * length(y_cal))
    des_freq2 <- round(test_freq2 * length(y_cal))

    idx1 <- unlist(lapply(labels, function(i) {
        category <- which(y_cal == i)
        n_i <- des_freq1[i]

        if (is.na(n_i) || length(category) == 0) {
            return(integer(0))
        }

        sample(category, size = n_i, replace = TRUE)
    }), use.names = FALSE)

    idx2 <- unlist(lapply(labels, function(i) {
        category <- which(y_cal == i)
        n_i <- des_freq2[i]

        if (is.na(n_i) || length(category) == 0) {
            return(integer(0))
        }

        sample(category, size = n_i, replace = TRUE)
    }), use.names = FALSE)

    list(
        p_cal1 = p_cal[idx1, , drop = FALSE],
        p_cal2 = p_cal[idx2, , drop = FALSE],
        p_test1 = test1,
        p_test2 = test2,
        y_cal1 = y_cal[idx1],
        y_cal2 = y_cal[idx2],
        idx = c(s, setdiff(seq_len(nrow(p_test)), s))
    ) # index in the original data
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
