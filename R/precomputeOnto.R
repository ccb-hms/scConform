#' Precompute ontology lookup structures for fast hierarchical inference
#'
#' Call this once after building your igraph ontology object. The result
#' can be passed to getPredictionSets() via the `onto_cache` argument.
#'
#' @param onto An igraph object representing the cell ontology.
#' @return A named list with precomputed structures:
#'   \item{dist_mat}{Full node x node distance matrix (rows = from, cols = to)}
#'   \item{leaf_children}{Named list: for each node, character vector of
#'     descendant leaf nodes}
#'   \item{ancestors_of}{Named list: for each node, character vector of
#'     ancestor nodes (including self)}
#'   \item{leaves}{Character vector of leaf node names}
#'   \item{root}{Character string: name of the root node}
#'   \item{node_names}{Character vector of all node names}
#' @examples
#' library(igraph)
#' # Let's build a random ontology
#' onto <- graph_from_literal(
#'     animal -+dog:cat, cat -+british:persian,
#'     dog -+cocker:retriever, retriever -+golden:labrador
#' )
#' # Precompute the ontology cache
#' onto_cache <- precomputeOnto(onto)
#' @importFrom igraph V distances degree
#' @export
precomputeOnto <- function(onto) {
    node_names <- V(onto)$name
    if (is.null(node_names) || length(node_names) == 0L || anyNA(node_names) || anyDuplicated(node_names)) {
        stop("Ontology vertices must have unique, non-missing names in V(onto)$name.")
    }

    # Full pairwise distance matrix (nNodes x nNodes)
    # dist_mat[i, j] = distance from node i to node j
    # Inf where not reachable
    dist_mat <- distances(onto, mode = "out")
    rownames(dist_mat) <- node_names
    colnames(dist_mat) <- node_names

    out_deg <- degree(onto, mode = "out")
    in_deg  <- degree(onto, mode = "in")

    leaves <- node_names[out_deg == 0]
    root_nodes <- node_names[in_deg == 0]
     if (length(root_nodes) != 1L) {
         stop("Ontology must have exactly one root node (in-degree 0). Found: ", length(root_nodes))
     }
     root <- root_nodes[[1]]

    # For each node: which leaves are reachable (descendants that are leaves)?
    # dist_mat[node, leaf] is finite  <=>  leaf is reachable from node
    leaf_children <- lapply(node_names, function(nd) {
        leaves[is.finite(dist_mat[nd, leaves])]
    })
    names(leaf_children) <- node_names

    # For each node: which nodes can reach it (ancestors including self)?
    # node j is an ancestor of nd  <=>  dist_mat[j, nd] is finite
    ancestors_of <- lapply(node_names, function(nd) {
        node_names[is.finite(dist_mat[, nd])]
    })
    names(ancestors_of) <- node_names

    list(
        dist_mat      = dist_mat,
        leaf_children = leaf_children,
        ancestors_of  = ancestors_of,
        leaves        = leaves,
        root          = root,
        node_names    = node_names
    )
}

