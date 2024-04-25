#' @title Get hierarchical prediction sets exploiting the cell ontology
#' @description Let K be the total number of distinct cell type labels and n, m
#' the number of cells in the calibration and in the test data, respectively.
#' This function takes as input two matrices: a matrix \code{n x K} and
#' a matrix \code{m x K} with the estimated
#' probabilities for each cell in the calibration and in the test data, respectively.
#' It returns a list with the prediction sets for each cell in the test data.
#'
#' @param p.cal a nxK matrix of category scores for the cells in the calibration data.
#' Needs to have colnames for each category
#' @param p.test a mxK matrix of category scores for the cells in the test data
#' @param y.cal a vector of length n with the labels of the cells in the calibration data
#' @param onto the considered section of the cell ontology
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @param lambda a vector of possible lambda values to be considered
#' @author Daniela Corbetta
#' @return The function \code{getHierarchicalPredSets} returns a list of length equal to the number of cells in the test data.
#' Each element of the list contains the prediction set for that cell.
#' @references For reference on conformal risk control, see
#' Angelopoulos, Anastasios N., et al. "Conformal risk control." arXiv preprint arXiv:2208.02814 (2022).
#' @importFrom foreach %dopar%
#' @export


getHierarchicalPredSets <- function(p.cal, p.test, y.cal, onto, alpha, lambdas){
  # Get prediction sets for each value of lambda for all the calibration data
  sets <- foreach(lambda = lambdas) %dopar% {
            lapply(1:nrow(p.cal),
              function(i) .predSets(lambda=lambda, pred=p.cal[i, ], onto=onto))}

  # Get the loss table (ncal x length(lambda) table with TRUE\FALSE)
  loss <- sapply(1:length(lambdas), function(lambda) {
    sapply(seq_along(y.cal), function(i) {
      !(y.cal[i] %in% prediction[[lambda]][[i]])
    })
  })

  # Get lhat
  n <- nrow(loss)
  rhat <- colMeans(loss)
  lhat_idx <- min(which(((n/(n+1)) * rhat + 1/(n+1) ) <= alpha))
  lhat <- lambdas[lhat_idx]

  # Get prediction sets for test data
  sets.test <- apply(p.test, 1, function(x) .predSets(lambda=lhat, pred=x, onto=onto))

  return(sets.test)
}








