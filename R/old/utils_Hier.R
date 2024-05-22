################################################################################
###### File with utils functions for conformal risk control on graph structure
################################################################################


# Conformal risk control
# Function to retrieve children nodes of a given interior node.
# If leaf=T gives leaf nodes, else gives all descendants
children <- function(node, graph, leaf = T) {
    require(igraph)
    if (leaf) {
        return(V(graph)$name[is.finite(distances(graph, node, mode = "out")) & igraph::degree(graph, mode = "out") == 0])
    } else {
        return(V(graph)$name[is.finite(distances(graph, node, mode = "out"))])
    }
}

# Function to retrieve ancestors of a given node
ancestors <- function(node, graph, include_self = T) {
    require(igraph)
    if (include_self) {
        return(V(graph)$name[is.finite(distances(graph, node, mode = "in"))])
    } else {
        return(V(graph)$name[is.finite(distances(graph, node, mode = "in")) & V(graph)$name != node])
    }
}

# Define function to compute g(a,x), where a is a node
# (g(a,x)=sum_{y in P(a)}(pi(y)))
# @param pred named vector with predicted probabilities for one observation to be in one of the leaves class

g <- function(pred, int_node, graph) {
    c <- children(node = int_node, graph = graph, leaf = T)
    return(sum(pred[c]))
}

# Function to construct confidence sets (for one obs for now)
# @param lambda probability bound
# @param p predicted probabilities for leaf nodes
# @param graph graph

pred_sets1 <- function(lambda, p, graph) {
    pred_class <- names(p)[which.max(p)]
    anc <- ancestors(node = pred_class, graph = graph, include_self = TRUE)

    # Obtain scores for ancestors of the predicted class
    scores <- sapply(as.character(anc), function(i) g(p, i, graph))
    names(scores) <- anc

    # Sort them by score and if there are ties by distance to the predicted class
    ## compute distance from predicted class
    pos <- distances(graph, v = anc, to = pred_class, mode = "out")
    tie_breaker <- as.vector(t(pos))
    names(tie_breaker) <- colnames(t(pos))
    sorted_indices <- order(scores, tie_breaker, decreasing = F)
    sorted_scores <- scores[sorted_indices]

    # Select the first score that is geq than lambda
    sel_node <- names(sorted_scores)[round(sorted_scores, 15) >= lambda][1]

    # c(lapply(anc[round(scores, 15) <= lambda]
    selected <- c(
        lapply(anc[round(scores, 15) < lambda], function(x) children(node = x, graph = graph)),
        list(children(sel_node, graph))
    )

    return(Reduce(union, selected))
}

# Do not consider this
pred_sets2 <- function(lambda, p, graph) {
    pred_class <- names(p)[which.max(p)]
    anc <- ancestors(node = pred_class, graph = graph, include_self = TRUE)

    scores <- sapply(anc, function(i) g(p, i, graph))
    names(scores) <- anc

    selected <- lapply(anc[round(scores, 15) <= lambda], function(x) children(node = x, graph = graph))
    if (length(selected) == 0) selected <- pred_class

    return(Reduce(union, selected))
}

# Loss table with miscoverage
get_loss_table_mis <- function(lambdas, prediction, ycal, graph) {
    loss <- sapply(1:length(lambdas), function(lambda) {
        sapply(seq_along(ycal), function(i) {
            !(ycal[i] %in% prediction[[lambda]][[i]])
        })
    })

    return(loss)
}

# Hierarchical loss proposed by Angelopoulus and Bates
hier_loss <- function(set, true_class, graph) {
    anc <- ancestors(true_class, graph)
    set_distances <- distances(graph, v = anc, to = set)

    vroot <- V(graph)[degree(graph, mode = "in") == 0]
    depth <- max(distances(graph, vroot, mode = "out"))

    return(min(set_distances) / depth)
}

# Loss table with hierarchical loss
get_loss_table <- function(lambdas, prediction, ycal, graph) {
    loss <- sapply(1:length(lambdas), function(lambda) {
        sapply(seq_along(ycal), function(i) {
            hier_loss(prediction[[lambda]][[i]], true_class = ycal[i], graph = graph)
        })
    })

    return(loss)
}

# Find lambda hat on a grid of lambda values
get_lhat <- function(calib_loss_table, lambdas, alpha, B = 1) {
    n <- nrow(calib_loss_table)
    rhat <- colMeans(calib_loss_table)
    lhat_idx <- min(which(((n / (n + 1)) * rhat + B / (n + 1)) <= alpha))
    # Return the corresponding lambda value at lhat_idx
    return(lambdas[lhat_idx])
}
