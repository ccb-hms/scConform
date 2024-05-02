#' @title Get prediction sets with split conformal inference
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
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @author Daniela Corbetta
#' @return The function \code{getConformalPredSets} returns a list of length equal to the number of cells in the test data.
#' Each element of the list contains the prediction set for that cell.
#' @importFrom stats quantile
#' @references For reference on split conformal prediction, refer to section 1 of
#' Angelopoulos, Anastasios N., and Stephen Bates. "A gentle introduction to conformal prediction and distribution-free uncertainty quantification." arXiv preprint arXiv:2107.07511 (2021).

.getConformalPredSets <- function(p.cal, p.test, y.cal, alpha){
    # Get calibration scores (1-predicted probability for the true class)
    true <- rep(NA, dim(p.cal)[1])
    for (i in 1:dim(p.cal)[1])
      true[i] <- p.cal[i, y.cal[i]]
    s <- 1-true

    # Get adjusted quantile
    n <- nrow(p.cal)
    q_level <- ceiling((n+1)*(1-alpha))/n
    qhat <- quantile(s, q_level)

    # Get prediction sets
    prediction_sets <- p.test >= 1-qhat
    pr.list <- lapply(1:nrow(prediction_sets), function(i) {
        colnames(prediction_sets)[prediction_sets[i, ]]
    })
    return(pr.list)
}





