test_that("resample=TRUE returns ordered prediction sets of correct length (conformal)", {
  labels <- c("A", "B", "C")

  set.seed(123)
  # p_test needs to be big enough so split isn't degenerate
  n_cal <- 60
  n_test <- 20

  # Make valid probability matrices
  p_cal <- matrix(runif(n_cal * length(labels)), nrow = n_cal)
  p_cal <- p_cal / rowSums(p_cal)
  colnames(p_cal) <- labels

  p_test <- matrix(runif(n_test * length(labels)), nrow = n_test)
  p_test <- p_test / rowSums(p_test)
  colnames(p_test) <- labels

  y_cal <- sample(labels, n_cal, replace = TRUE)

  res <- getPredictionSets(
    x_query = p_test,
    x_cal = p_cal,
    y_cal = y_cal,
    onto = NULL,
    follow_ontology = FALSE,
    resample = TRUE,
    labels = labels,
    return_sc = FALSE,
    alpha = 0.2
  )

  expect_type(res, "list")
  expect_length(res, n_test)
  expect_true(all(vapply(res, function(s) all(s %in% labels), logical(1))))
})
