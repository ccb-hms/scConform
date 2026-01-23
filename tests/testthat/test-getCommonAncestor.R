test_that("getCommonAncestor returns correct ancestor on a toy ontology", {
  skip_if_not_installed("igraph")

  onto <- make_toy_onto()

  # All within "dog" subtree -> LCA should be "dog"
  pred_set <- c("golden", "labrador", "cocker")
  expect_equal(getCommonAncestor(pred_set, onto), "dog")

  # Two cats -> LCA is "cat"
  pred_set2 <- c("british", "persian")
  expect_equal(getCommonAncestor(pred_set2, onto), "cat")

  # Across dog & cat -> LCA is "animal"
  pred_set3 <- c("golden", "persian")
  expect_equal(getCommonAncestor(pred_set3, onto), "animal")
})

test_that("getCommonAncestor works when pred_set is length 1", {
  skip_if_not_installed("igraph")

  onto <- make_toy_onto()
  expect_equal(getCommonAncestor("golden", onto), "golden")
})
