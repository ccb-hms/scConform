test_that("plotResult runs (suppressing possible graphviz warnings)", {
    skip_if_not_installed("Rgraphviz")

    onto <- make_toy_onto()

    pdf(NULL)

    expect_silent(
        suppressWarnings(
            plotResult(
                pred_set = c("golden", "labrador"),
                onto = onto,
                onto_cache = precomputeOnto(onto),
                probs = NULL,
                add_scores = FALSE,
                title = "test"
            )
        )
    )
})

test_that("plotResult with probs and add_scores runs (suppressing possible graphviz warnings)", {
    skip_if_not_installed("Rgraphviz")

    onto <- make_toy_onto()
    leaves <- igraph::V(onto)$name[igraph::degree(onto, mode = "out") == 0]

    probs <- setNames(rep(0, length(leaves)), leaves)
    probs["golden"] <- 0.7
    probs["labrador"] <- 0.2
    probs["cocker"] <- 0.1

    pdf(NULL)

    expect_silent(
        suppressWarnings(
            plotResult(
                pred_set = c("golden", "labrador"),
                onto = onto,
                onto_cache = precomputeOnto(onto),
                probs = probs,
                add_scores = TRUE,
                title = "test"
            )
        )
    )
})
