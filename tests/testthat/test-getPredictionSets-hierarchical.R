test_that("getPredictionSets errors when follow_ontology=TRUE but onto is NULL", {
  p_cal  <- make_toy_probs(n = 10)
  p_test <- make_toy_probs(n = 3, labels = colnames(p_cal))
  y_cal  <- sample(colnames(p_cal), 10, replace = TRUE)

  expect_error(
    getPredictionSets(
      x_query = p_test, x_cal = p_cal, y_cal = y_cal,
      onto = NULL, follow_ontology = TRUE,
      return_sc = FALSE,
      labels = colnames(p_cal)
    ),
    "An ontology is required"
  )
})

test_that("method argument is validated", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()

  p_cal  <- make_toy_probs(n = 20, labels = c("golden","labrador","cocker","persian"))
  p_test <- make_toy_probs(n = 4, labels = colnames(p_cal))
  y_cal  <- sample(colnames(p_cal), 20, replace = TRUE)

  expect_error(
    getPredictionSets(
      x_query = p_test, x_cal = p_cal, y_cal = y_cal,
      onto = onto, follow_ontology = TRUE,
      method = "nope",
      return_sc = FALSE
    ),
    "must be one of"
  )

  expect_error(
    getPredictionSets(
      x_query = p_test, x_cal = p_cal, y_cal = y_cal,
      onto = onto, follow_ontology = TRUE,
      method = 123,
      return_sc = FALSE
    ),
    "either a character string or a function"
  )
})

test_that("hierarchical prediction sets are subsets of ontology leaves", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()

  # Use only leaf labels present in ontology for matrix columns
  labels <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  p_cal  <- make_toy_probs(n = 40, labels = labels)
  p_test <- make_toy_probs(n = 6, labels = labels)
  y_cal  <- sample(labels, 40, replace = TRUE)

  res <- getPredictionSets(
    x_query = p_test, x_cal = p_cal, y_cal = y_cal,
    onto = onto, alpha = 0.2,
    follow_ontology = TRUE,
    method = "full",
    return_sc = FALSE,
    lambdas = seq(0.2, 0.9, length.out = 5)  # keep it small/fast
  )

  expect_type(res, "list")
  expect_length(res, nrow(p_test))
  expect_true(all(vapply(res, function(s) all(s %in% labels), logical(1))))
  # non-empty sets are expected for "full"
  expect_true(all(vapply(res, function(s) length(s) >= 1, logical(1))))
})

test_that("method='rank' rejects lambda outside [0,1] (via calibration), so use valid lambdas", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()
  labels <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  p_cal  <- make_toy_probs(n = 30, labels = labels)
  p_test <- make_toy_probs(n = 5, labels = labels)
  y_cal  <- sample(labels, 30, replace = TRUE)

  res <- getPredictionSets(
    x_query = p_test, x_cal = p_cal, y_cal = y_cal,
    onto = onto, alpha = 0.2,
    follow_ontology = TRUE,
    method = "rank",
    lambdas = seq(0.2, 0.8, length.out = 4),
    return_sc = FALSE
  )
  expect_length(res, nrow(p_test))
})

test_that("simplify returns common ancestor only when it is not a ramification", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()
  labels <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  p_cal  <- make_toy_probs(n = 30, labels = labels)
  p_test <- make_toy_probs(n = 3, labels = labels)
  y_cal  <- sample(labels, 30, replace = TRUE)

  res <- getPredictionSets(
    x_query = p_test, x_cal = p_cal, y_cal = y_cal,
    onto = onto, alpha = 0.2,
    follow_ontology = TRUE,
    method = "full",
    lambdas = seq(0.2, 0.9, length.out = 4),
    simplify = TRUE,
    return_sc = FALSE
  )

  # After simplify, elements can be a single node name (ancestor) or a leaf set.
  # Check that any singletons are valid ontology nodes.
  all_nodes <- igraph::V(onto)$name
  ok <- vapply(res, function(x) {
    if (length(x) == 1) x %in% all_nodes else all(x %in% labels)
  }, logical(1))
  expect_true(all(ok))
})

test_that("method = 'step' rejects negative lambdas", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()
  labels <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  p_cal  <- make_toy_probs(n = 10, labels = labels)
  p_test <- make_toy_probs(n = 2, labels = labels)
  y_cal  <- sample(labels, 10, replace = TRUE)

  expect_error(
    getPredictionSets(
      x_query = p_test,
      x_cal   = p_cal,
      y_cal   = y_cal,
      onto    = onto,
      follow_ontology = TRUE,
      method  = "step",
      lambdas = c(-2, -1),
      return_sc = FALSE
    ),
    "non-negative integer"
  )
})

test_that("method = 'rank' rejects lambdas outside [0,1]", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()
  labels <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  p_cal  <- make_toy_probs(n = 10, labels = labels)
  p_test <- make_toy_probs(n = 2, labels = labels)
  y_cal  <- sample(labels, 10, replace = TRUE)

  expect_error(
    getPredictionSets(
      x_query = p_test,
      x_cal   = p_cal,
      y_cal   = y_cal,
      onto    = onto,
      follow_ontology = TRUE,
      method  = "rank",
      lambdas = c(-0.1, 0.2),
      return_sc = FALSE
    ),
    "\\[0, 1\\]"
  )
})

test_that("method = 'rank' errors when no ontology leaves match prediction names", {
  skip_if_not_installed("igraph")
  onto <- make_toy_onto()

  # prediction matrix with labels NOT in ontology
  p_cal  <- make_toy_probs(n = 10, labels = c("X", "Y"))
  p_test <- make_toy_probs(n = 2,  labels = c("X", "Y"))
  y_cal  <- sample(c("X", "Y"), 10, replace = TRUE)

  expect_error(
    getPredictionSets(
      x_query = p_test,
      x_cal   = p_cal,
      y_cal   = y_cal,
      onto    = onto,
      follow_ontology = TRUE,
      method  = "rank",
      lambdas = c(0.5),
      return_sc = FALSE
    ),
    "No ontology leaves are present"
  )
})

test_that("hierarchical (full) calibration always selects a lambda when calibration is perfect", {
  onto <- make_toy_onto()
  leaves <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

  n_cal <- 50
  set.seed(1)
  y_cal <- sample(leaves, n_cal, replace = TRUE)

  # Perfect calibration: true label prob = 1
  p_cal <- matrix(0, nrow = n_cal, ncol = length(leaves))
  colnames(p_cal) <- leaves
  for (i in seq_len(n_cal)) p_cal[i, y_cal[i]] <- 1

  p_test <- make_toy_probs(n = 10, labels = leaves)

  expect_message(
    res <- getPredictionSets(
      x_query = p_test,
      x_cal = p_cal,
      y_cal = y_cal,
      onto = onto,
      follow_ontology = TRUE,
      method = "full",
      lambdas = seq(0.2, 0.9, length.out = 5),
      alpha = 0.2,
      return_sc = FALSE
    ),
    "Selected lambda_hat"
  )

  expect_length(res, nrow(p_test))
  expect_true(all(vapply(res, function(s) all(s %in% leaves), logical(1))))
})

