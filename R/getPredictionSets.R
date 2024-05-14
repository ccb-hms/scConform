#' @title Get conformal prediction sets
#' @description Let K be the total number of distinct cell type labels and n, m
#' the number of cells in the calibration and in the test data, respectively.
#' This function takes as input two matrices: a matrix \code{n x K} and
#' a matrix \code{m x K} with the estimated
#' probabilities for each cell in the calibration and in the test data,
#' respectively. It returns a list with the prediction sets for each cell in
#' the test data.
#'
#' @param x.query query data for which we want to build prediction sets. Could
#' be either a SingleCellExperiment object with the estimated probabilities for
#' each cell type in the colData, or a matrix of dimension \code{n x K},
#' where \code{n} is the number of cells and \code{K} is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param x.cal calibration data. Could be either a
#' SingleCellExperiment object with the estimated probabilities for each cell
#' type in the colData, or a named matrix of dimension \code{m x K}, where
#' \code{m} is the number of cells and \code{K} is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param y.cal a vector of length n with the labels of the cells in the
#' calibration data
#' @param onto the considered section of the cell ontology as an igraph object.
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @param lambdas a vector of possible lambda values to be considered. Necessary
#' only when \code{follow.ontology=TRUE}.
#' @param follow.ontology If \code{TRUE}, then the function returns hierarchical
#' prediction sets that follow the cell ontology structure. If \code{FALSE}, it
#' returns classical conformal prediction sets.
#' @param resample.cal Should the calibration dataset be resampled according to
#' the estimated relative frequencies of cell types in the query data?
#' @param labels labels of different considered cell types. Necessary if
#' \code{onto=NULL}, otherwise they are set to the leaf nodes of the provided
#' graph
#' @param return.sc if \code{TRUE}, the prediction sets are returned as a
#' SingleCellExperiment. If \code{FALSE}, they are returned as list. By default,
#' it is set to \code{TRUE} when \code{x.query} is a SingleCellExperiment or
#' SpatialExperiment object, to \code{FALSE} when \code{x.query} is a matrix.
#' @param pr.name name of the colData variable that will contain the prediction
#' sets. Default \code{pred.set}.
#' @author Daniela Corbetta
#' @return The function \code{getPredictionSets} returns
#' \item{\code{return.sc=TRUE}}{a SingleCellExperiment or SpatialExperiment
#' object with the prediction sets in the colData. The name of the variable
#' containing the prediction sets is given by the parameter \code{pr.name}}
#' \item{\code{return.sc=FALSE}}{a list of length equal
#' to the number of cells in the test data. Each element of the list contains
#' the prediction set for that cell.}
#' @references For an introduction to conformal prediction, see of
#' Angelopoulos, Anastasios N., and Stephen Bates. "A gentle introduction to
#' conformal prediction and distribution-free uncertainty quantification."
#' arXiv preprint arXiv:2107.07511 (2021).
#' For reference on conformal risk control, see
#' Angelopoulos, Anastasios N., et al. "Conformal risk control."
#' arXiv preprint arXiv:2208.02814 (2022).
#' @importFrom foreach %dopar%
#' @importFrom foreach foreach
#' @importFrom SummarizedExperiment colData
#' @importFrom igraph V
#' @importFrom igraph degree
#' @importFrom igraph distances
#' @importFrom stats quantile
#' @export

getPredictionSets <- function(x.query, x.cal, y.cal, onto=NULL, alpha = 0.1,
                              lambdas = seq(0.001,0.999,length.out=100),
                              follow.ontology=TRUE,
                              resample.cal=FALSE,
                              labels=NULL,
                              return.sc=NULL,
                              pr.name="pred.set"){
    if(follow.ontology & is.null(onto)){
        stop("An ontology is required for hierarchical prediction set.
             Please provide one or ask for conformal prediction set
             (follow.ontology=FALSE)")
    }
    if(is.null(onto) & is.null(labels)){
        stop("Please provide cell labels (labels parameter)")
    }
    if(!is.null(return.sc)){
        if(return.sc==TRUE & is.matrix(x.query)){
            stop("If x.query is a matrix output has to be a list
                 (return.sc=FALSE)")
        }
    }
    if(isa(x.query, "SpatialExperiment") |
       isa(x.query, "SingleCellExperiment") |
       isa(x.query, "SummarizedExperiment"))
        sc <- TRUE
    else if(is.matrix(x.query))
        sc <- FALSE
    else stop("Please provide as input in x.query a SpatialExperiment,
              SingleCellExperiment or a matrix")

    # Retrieve labels from the ontology (need to add retrieval from y.cal/data
    # when follow.ontology=FALSE)
    if(is.null(labels))
        labels <- V(onto)$name[degree(onto, mode="out")==0]
    K <- length(labels)

    # If input is not a matrix, retrieve prediction matrix from colData
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

    if(!resample.cal){
        if (follow.ontology){
            pred.sets <- .getHierarchicalPredSets(p.cal=p.cal, p.test=p.query,
                                                  y.cal=y.cal, onto=onto,
                                                  alpha=alpha,
                                                  lambdas=lambdas)
        }
        else
            pred.sets <- .getConformalPredSets(p.cal=p.cal, p.test=p.query,
                                               y.cal=y.cal, alpha=alpha)
    }

    if(resample.cal){
        data <- resample.two(p.cal=p.cal, p.test=p.query, y.cal=y.cal,
                             labels=labels)
        if (follow.ontology){
            pred.sets1 <- .getHierarchicalPredSets(p.cal=data$p.cal2,
                                                   p.test=data$p.test1,
                                                   y.cal=data$y.cal2,
                                                   onto=onto,
                                                   alpha=alpha,
                                                   lambdas=lambdas)
            pred.sets2 <- .getHierarchicalPredSets(p.cal=data$p.cal1,
                                                   p.test=data$p.test2,
                                                   y.cal=data$y.cal1,
                                                   onto=onto,
                                                   alpha=alpha,
                                                   lambdas=lambdas)
            pred.sets <- c(pred.sets1, pred.sets2)
        }
        else {
            pred.sets1 <- .getConformalPredSets(p.cal=data$p.cal2,
                                                p.test=data$p.test1,
                                                y.cal=data$y.cal2,
                                                alpha=alpha)
            pred.sets2 <- .getConformalPredSets(p.cal=data$p.cal1,
                                                p.test=data$p.test2,
                                                y.cal=data$y.cal1,
                                                alpha=alpha)
            pred.sets <- c(pred.sets1, pred.sets2)
        }
        # Order the prediction set
        pred.sets <- pred.sets[order(data$idx)]
    }

    # if not specified, return a sc object if the input was a sc object,
    # a matrix if the input was a matrix
    if(is.null(return.sc) & sc)
        return.sc <- TRUE
    if(is.null(return.sc) & !sc)
        return.sc <- FALSE
    if(return.sc){
        colData(x.query)[[pr.name]] <- pred.sets
        return(x.query)
    }
    if(!return.sc){
        return(pred.sets)
    }

}
