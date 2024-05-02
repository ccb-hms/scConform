#' @title Get prediction sets
#' @description Let K be the total number of distinct cell type labels and n, m
#' the number of cells in the calibration and in the test data, respectively.
#' This function takes as input two matrices: a matrix \code{n x K} and
#' a matrix \code{m x K} with the estimated
#' probabilities for each cell in the calibration and in the test data, respectively.
#' It returns a list with the prediction sets for each cell in the test data.
#'
#' @param x.query query data for which we want to build prediction sets. Could be either a
#' SingleCellExperiment object with the estimated probabilities for each cell type
#' in the colData, or a matrix of dimension \code{n x K}, where \code{n} is the number
#' of cells and \code{K} is the number of different labels. The colnames of the matrix
#' have to correspond to the cell labels.
#' @param x.cal calibration data. Could be either a
#' SingleCellExperiment object with the estimated probabilities for each cell type
#' in the colData, or a named matrix of dimension \code{m x K}, where \code{m} is the number
#' of cells and \code{K} is the number of different labels. The colnames of the matrix
#' have to correspond to the cell labels.
#' @param y.cal a vector of length n with the labels of the cells in the calibration data
#' @param onto the considered section of the cell ontology as an igraph object.
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @param lambdas a vector of possible lambda values to be considered
#' @param follow_ontology If \code{TRUE}, then the function returns hierarchical
#' prediction sets that follow the cell ontology structure. If \code{FALSE}, it
#' returns classical conformal prediction sets.
#' @author Daniela Corbetta
#' @return The function \code{getPredictionSets} returns a list of length equal
#' to the number of cells in the test data.
#' Each element of the list contains the prediction set for that cell.
#' @references For an introduction to conformal prediction, see of
#' Angelopoulos, Anastasios N., and Stephen Bates. "A gentle introduction to
#' conformal prediction and distribution-free uncertainty quantification." arXiv preprint arXiv:2107.07511 (2021).
#' For reference on conformal risk control, see
#' Angelopoulos, Anastasios N., et al. "Conformal risk control." arXiv preprint arXiv:2208.02814 (2022).
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom SummarizedExperiment colData
#' @export

getPredictionSets <- function(x.query, x.cal, y.cal, onto=NULL, alpha = 0.1,
                              lambdas = lambdas <- seq(0.001,0.999,length.out=100),
                              follow_ontology=TRUE){

    # Add check to see if x.cal, x.query are SingleCell/SpatialExperiment or matrices
    # Retrieve labels from the ontology (need to add retrieval from y.cal/data
    # when follow_ontology=FALSE)
    labels <- V(onto)$name[degree(onto, mode="out")==0]
    K <- length(labels)
    if(!is.matrix(x.query)){
        n.query <- ncol(x.query)
        p.query <- matrix(NA, nrow=n.query, ncol=K)
        colnames(p.query) <- labels
        for(i in labels){
          p.query[,i] <- colData(x.query)[[i]]
        }
    }
    else p.query <- x.query

    if(!is.matrix(x.cal)){
      n.cal <- ncol(x.cal)
      p.cal <- matrix(NA, nrow=n.cal, ncol=K)
      colnames(p.cal) <- labels
      for(i in labels){
        p.cal[,i] <- colData(x.cal)[[i]]
      }
    }
    else p.cal <- x.cal

    if (follow_ontology){
        pred.sets <- .getHierarchicalPredSets(p.cal=p.cal, p.test=p.query,
                                              y.cal=y.cal, onto=onto, alpha=alpha,
                                              lambdas=lambdas)$sets.test
    }
    else
        pred.sets <- .getConformalPredSets(p.cal=p.cal, p.test=p.query,
                                           y.cal=y.cal, alpha=alpha)
    return(pred.sets)
}
