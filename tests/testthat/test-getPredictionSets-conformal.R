test_that("getPredictionSets errors when follow_ontology=FALSE and simplify=TRUE", {
    p_cal <- make_toy_probs(n = 10)
    p_test <- make_toy_probs(n = 3, labels = colnames(p_cal))
    y_cal <- sample(colnames(p_cal), 10, replace = TRUE)

    expect_error(
        getPredictionSets(
            x_query = p_test, x_cal = p_cal, y_cal = y_cal,
            onto = NULL, alpha = 0.2,
            follow_ontology = FALSE,
            simplify = TRUE,
            labels = colnames(p_cal),
            return_sc = FALSE
        ),
        "simplify=FALSE"
    )
})

test_that("getPredictionSets returns list of leaf labels in conformal mode", {
    p_cal <- make_toy_probs(n = 50)
    p_test <- make_toy_probs(n = 5, labels = colnames(p_cal))
    y_cal <- sample(colnames(p_cal), 50, replace = TRUE)

    res <- getPredictionSets(
        x_query = p_test, x_cal = p_cal, y_cal = y_cal,
        onto = NULL, alpha = 0.1,
        follow_ontology = FALSE,
        resample = FALSE,
        labels = colnames(p_cal),
        return_sc = FALSE
    )

    expect_type(res, "list")
    expect_length(res, nrow(p_test))
    expect_true(all(vapply(res, is.character, logical(1))))
    # Each set should be subset of labels
    expect_true(all(vapply(res, function(s) all(s %in% colnames(p_test)), logical(1))))
})

test_that("getPredictionSets requires labels when onto is NULL", {
    p_cal <- make_toy_probs(n = 10)
    p_test <- make_toy_probs(n = 3, labels = colnames(p_cal))
    y_cal <- sample(colnames(p_cal), 10, replace = TRUE)

    expect_error(
        getPredictionSets(
            x_query = p_test, x_cal = p_cal, y_cal = y_cal,
            onto = NULL, follow_ontology = FALSE,
            return_sc = FALSE
        ),
        "Please provide cell labels"
    )
})

test_that("matrix input cannot return SingleCellExperiment", {
    p_cal <- make_toy_probs(n = 10)
    p_test <- make_toy_probs(n = 3, labels = colnames(p_cal))
    y_cal <- sample(colnames(p_cal), 10, replace = TRUE)

    expect_error(
        getPredictionSets(
            x_query = p_test, x_cal = p_cal, y_cal = y_cal,
            onto = NULL, alpha = 0.2,
            follow_ontology = FALSE,
            labels = colnames(p_cal),
            return_sc = TRUE
        ),
        "output has to be a list"
    )
})

test_that("invalid x_query type is rejected", {
    expect_error(
        getPredictionSets(
            x_query = list(a = 1),
            x_cal = matrix(0.5, nrow = 1, ncol = 1),
            y_cal = "A",
            labels = "A",
            follow_ontology = FALSE,
            return_sc = FALSE
        ),
        "Please provide as input in x_query"
    )
})

test_that("return_sc = TRUE returns a SingleCellExperiment with prediction sets in colData", {
    skip_if_not_installed("SingleCellExperiment")

    library(SingleCellExperiment)

    labels <- c("A", "B", "C")

    ## Construct a tiny SingleCellExperiment
    set.seed(1)
    n_cells <- 5
    p <- matrix(runif(n_cells * length(labels)), nrow = n_cells)
    p <- p / rowSums(p)
    colnames(p) <- labels

    sce <- SingleCellExperiment(
        assays = list(counts = matrix(0, nrow = 1, ncol = n_cells))
    )

    ## Put predicted probabilities in colData (as expected by .retrievePredMatrix)
    for (lab in labels) {
        colData(sce)[[lab]] <- p[, lab]
    }

    ## Calibration data (matrix is fine)
    p_cal <- matrix(runif(30 * length(labels)), nrow = 30)
    p_cal <- p_cal / rowSums(p_cal)
    colnames(p_cal) <- labels
    y_cal <- sample(labels, 30, replace = TRUE)

    res <- getPredictionSets(
        x_query = sce,
        x_cal = p_cal,
        y_cal = y_cal,
        onto = NULL,
        follow_ontology = FALSE,
        labels = labels,
        return_sc = TRUE
    )

    ## Assertions
    expect_s4_class(res, "SingleCellExperiment")
    expect_true("pred.set" %in% colnames(colData(res)))

    ps <- colData(res)[["pred.set"]]
    expect_type(ps, "list")
    expect_length(ps, n_cells)
    expect_true(all(vapply(ps, is.character, logical(1))))
})
