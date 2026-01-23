make_toy_onto <- function() {
  igraph::graph_from_literal(
    animal -+ dog:cat,
    cat    -+ british:persian,
    dog    -+ cocker:retriever,
    retriever -+ golden:labrador
  )
}

make_toy_probs <- function(n = 6, labels = c("golden","labrador","cocker","persian")) {
  set.seed(1)
  p <- matrix(stats::runif(n * length(labels)), ncol = length(labels))
  p <- p / rowSums(p)
  colnames(p) <- labels
  p
}

