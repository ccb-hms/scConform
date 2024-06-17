#' @title Get conformal prediction sets
#' @description This function returns prediction sets for the cell
#' type of cells in a SingleCellExperiment objects.
#' It implements two methods: the first one uses standard conformal inference,
#' while the second one conformal risk control (see details). The output is
#' either a SingleCellExperiment object with the prediction sets in the colData
#' or a list.
#'
#' @param x_query query data for which we want to build prediction sets. Could
#' be either a SingleCellExperiment object with the estimated probabilities for
#' each cell type in the colData, or a named matrix of dimension \code{n x K},
#' where \code{n} is the number of cells and \code{K} is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param x_cal calibration data. Could be either a
#' SingleCellExperiment object with the estimated probabilities for each cell
#' type in the colData, or a named matrix of dimension \code{m x K}, where
#' \code{m} is the number of cells and \code{K} is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param y_cal a vector of length \code{m} with the true labels of the cells in
#' the calibration data.
#' @param onto the considered section of the cell ontology as an igraph object.
#' @param alpha a number between 0 and 1 that indicates the allowed miscoverage
#' @param lambdas a vector of possible lambda values to be considered. Necessary
#' only when \code{follow_ontology=TRUE}.
#' @param follow_ontology If \code{TRUE}, then the function returns hierarchical
#' prediction sets that follow the cell ontology structure. If \code{FALSE}, it
#' returns classical conformal prediction sets. See details.
#' @param resample_cal Should the calibration dataset be resampled according to
#' the estimated relative frequencies of cell types in the query data?
#' @param labels labels of different considered cell types. Necessary if
#' \code{onto=NULL}, otherwise they are set equal to the leaf nodes of the
#' provided graph.
#' @param return_sc parameter the controls the output. If \code{TRUE}, the
#' function returns a SingleCellExperiment.
#' If \code{FALSE}, the function returns a list. By default,
#' it is set to \code{TRUE} when \code{x_query} is a SingleCellExperiment or
#' SpatialExperiment object and to \code{FALSE} when \code{x_query} is a matrix.
#' @param pr_name name of the colData variable in the returned
#' SingleCellExperiment object that will contain the prediction
#' sets. The default name is \code{pred.set}.
#' @param simplify_pred if \code{TRUE}, the output will be the common ancestor
#' of the labels inserted into the prediction set. If \code{FALSE} (default),
#' the output will be the set of the leaf labels.
#' @param BPPARAM BiocParallel instance for parallel computing. Default is
#' \code{SerialParam()}.
#' @return
#' \item{\code{return_sc = TRUE}}{the function \code{getPredictionSets} returns
#' a SingleCellExperiment or SpatialExperiment
#' object with the prediction sets in the colData. The name of the variable
#' containing the prediction sets is given by the parameter \code{pr_name}}
#' \item{\code{return_sc = FALSE}}{the function \code{getPredictionSets} returns
#' a list of length equal
#' to the number of cells in the test data. Each element of the list contains
#' the prediction set for that cell.}
#' @details
#' \subsection{Split conformal sets}{Conformal inference is a statistical
#' framework that allows to build
#' prediction sets for any probabilistic or machine learning model. Suppose we
#' have a classification task with \eqn{K} classes. We fit a classification
#' model \eqn{\hat{f}} that outputs estimated probabilities for each class:
#' \eqn{\hat{f}(x) \in [0,1]^K}. Split conformal inference requires to reserve a
#' portion of the labelled training data, \eqn{(X_1,Y_1),\dots, (X_n,Y_n)}, to
#' be used as calibration data. Given \eqn{\hat{f}} and the calibration data,
#' the objective of conformal inference is to build, for a new observation
#' \eqn{X_{n+1},} a prediction set \eqn{C(X_{n+1}) \subseteq\{1,\dots,K\}} that
#' satisfies
#' \deqn{P\left(Y_{n+1}\in C(X_{n+1})\right) \geq 1-\alpha}
#' for a user-chosen error rate \eqn{\alpha}. Note that conformal inference is
#' distribution-free and the sets provided have finite-samples validity.
#' The only assumption is that the test data and the calibration data are
#' exchangeable. The algorithm of split conformal inference is the following:
#' \enumerate{
#'   \item For the data in the calibration set, \eqn{(X_1,Y_1),\dots, (X_n,Y_n)}
#'   , obtain the \emph{conformal scores}, \eqn{s_i=1-\hat{f}(X_i)_{Y_i},
#'   \;i=1,\dots,n}. These scores will be high when the model is assigning a
#'   small probability to the true class, and low otherwise.
#'   \item Obtain \eqn{\hat{q}}, the
#'   \eqn{\lceil(1-\alpha)(n+1)\rceil/n} empirical quantile of the conformal
#'   scores.
#'   \item Finally, for a new observation \eqn{X_{n+1}}, construct a prediction
#'   set by including all the classes for which the estimated probability is
#'   higher than \eqn{1-\hat{q}}:
#'   \deqn{C(X_{n+1})=\{y: \hat{f}(X_{n+1})_y\geq 1-\hat{q}\}.}
#' }}
#' \subsection{Hierarchical conformal sets}{
#' Let \eqn{\hat{y}(x)} be the class with maximum estimated probability.
#' Moreover, given a directed graph let \eqn{\mathcal{P}(v)} and
#' \eqn{\mathcal{A}(v)} be the set on children nodes and ancestor nodes of
#' \eqn{v}, respectively. Finally, for each node \eqn{v} define a score
#' \eqn{g(v,x)} as the sum of the predicted probabilities of the leaf nodes that
#' are children of \eqn{v}.
#' To build the sets we propose the following algorithm: \deqn{\mathcal{P}(v)
#' \cup \{\mathcal{P}(a): a\in\mathcal{A}(\hat{y}(x)): g(a,x)\leq\lambda \},}
#' where \eqn{v:v\in \mathcal{A}(\hat{y}(x)), \;g(v,x)\geq\lambda,\;
#' v=\arg\min_{u:g(u,x)\geq\lambda}g(u,x)}.
#' In words, we start from the predicted class and we go up in the graph until
#' we find an ancestor of \eqn{\hat{y}(x)} that has a score that is at least
#' \eqn{\lambda} and include in the prediction sets all its children.
#' For theoretical reasons, to this subgraph we add all the other
#' ones that contain \eqn{\hat{y}(x)} for which the score is less than
#' \eqn{\lambda}. To choose \eqn{\lambda}, we follow eq. (4) in Angelopoulus et
#' al. (2023), considering the miscoverage as loss function. In this way, it is
#' still guaranteed that
#' \deqn{P(Y_{n+1}\notin C_\lambda (X_{n+1})) \leq \alpha.}}
#' @references For an introduction to conformal prediction, see
#' Angelopoulos, Anastasios N., and Stephen Bates. "A gentle introduction to
#' conformal prediction and distribution-free uncertainty quantification."
#' arXiv preprint arXiv:2107.07511 (2021).
#' For reference on conformal risk control, see
#' Angelopoulos, Anastasios N., et al. "Conformal risk control."
#' arXiv preprint arXiv:2208.02814 (2023).
#' @examples
#' # random p matrix
#' set.seed(1040)
#' p <- matrix(rnorm(2000 * 4), ncol = 4)
#' # Normalize the matrix p to have all numbers between 0 and 1 that sum to 1
#' # by row
#' p <- exp(p - apply(p, 1, max))
#' p <- p / rowSums(p)
#' cell_types <- c("T (CD4+)", "T (CD8+)", "B", "NK")
#' colnames(p) <- cell_types
#'
#' # Take 1000 rows as calibration and 1000 as test
#' p_cal <- p[1:1000, ]
#' p_test <- p[1001:2000, ]
#'
#' # Randomly create the vector of real cell types for p_cal and p_test
#' y_cal <- sample(cell_types, 1000, replace = TRUE)
#' y_test <- sample(cell_types, 1000, replace = TRUE)
#'
#' # Obtain conformal prediction sets
#' conf_sets <- getPredictionSets(
#'     x_query = p_test,
#'     x_cal = p_cal,
#'     y_cal = y_cal,
#'     onto = NULL,
#'     alpha = 0.1,
#'     follow_ontology = FALSE,
#'     resample_cal = FALSE,
#'     labels = cell_types,
#'     return_sc = FALSE
#' )
#'
#' @importFrom foreach %dopar% foreach
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom igraph V distances degree
#' @importFrom stats quantile
#' @importFrom BiocParallel SerialParam bplapply
#' @export

getPredictionSets <- function(
        x_query, x_cal, y_cal, onto = NULL, alpha = 0.1,
        lambdas = seq(0.001, 0.999, length.out = 100),
        follow_ontology = TRUE,
        resample_cal = FALSE,
        labels = NULL,
        return_sc = NULL,
        pr_name = "pred.set",
        simplify_pred = FALSE,
        BPPARAM = SerialParam()) {
    ## Sanity checks

    if (follow_ontology & is.null(onto)) {
        stop("An ontology is required for hierarchical prediction set.
             Please provide one or ask for conformal prediction set
             (follow_ontology=FALSE)")
    }

    if (is.null(onto) & is.null(labels)) {
        stop("Please provide cell labels (labels parameter)")
    }

    if (isa(x_query, "SpatialExperiment") |
        isa(x_query, "SingleCellExperiment") |
        isa(x_query, "SummarizedExperiment")) {
        sc <- TRUE
    } else if (is.matrix(x_query)) {
        sc <- FALSE
    } else {
        stop("Please provide as input in x_query a SpatialExperiment,
              SingleCellExperiment or a matrix")
    }

    if (!is.null(return_sc)) {
        if (return_sc == TRUE & !sc) {
            stop("If x_query is a matrix output has to be a list
                 (return_sc=FALSE)")
        }
    }

    if (!follow_ontology & simplify_pred) {
        stop("If follow_ontology=FALSE, please set simplify_pred=FALSE")
    }

    ## If labels parameter is NULL, retrieve labels from the ontology
    if (is.null(labels)) {
        labels <- V(onto)$name[degree(onto, mode = "out") == 0]
    }
    K <- length(labels)

    ## If input is not a matrix, retrieve prediction matrix from colData
    if (!is.matrix(x_query)) {
        p_query <- .retrievePredMatrix(x_query, K = K, labels = labels)
    } else {
        p_query <- x_query
    }

    if (!is.matrix(x_cal)) {
        p_cal <- .retrievePredMatrix(x_cal, K = K, labels = labels)
    } else {
        p_cal <- x_cal
    }

    if (!resample_cal) {
        if (follow_ontology) {
            pred_sets <- .getHierarchicalPredSets(
                p_cal = p_cal, p_test = p_query,
                y_cal = y_cal, onto = onto,
                alpha = alpha,
                lambdas = lambdas,
                BPPARAM = BPPARAM
            )
        } else {
            pred_sets <- .getConformalPredSets(
                p_cal = p_cal, p_test = p_query,
                y_cal = y_cal, alpha = alpha
            )
        }
    }

    if (resample_cal) {
        data <- .resampleTwo(
            p_cal = p_cal, p_test = p_query, y_cal = y_cal,
            labels = labels
        )
        if (follow_ontology) {
            pred_sets1 <- .getHierarchicalPredSets(
                p_cal = data$p_cal2,
                p_test = data$p_test1,
                y_cal = data$y_cal2,
                onto = onto,
                alpha = alpha,
                lambdas = lambdas,
                BPPARAM = BPPARAM
            )
            pred_sets2 <- .getHierarchicalPredSets(
                p_cal = data$p_cal1,
                p_test = data$p_test2,
                y_cal = data$y_cal1,
                onto = onto,
                alpha = alpha,
                lambdas = lambdas,
                BPPARAM = BPPARAM
            )
            pred_sets <- c(pred_sets1, pred_sets2)
        } else {
            pred_sets1 <- .getConformalPredSets(
                p_cal = data$p_cal2,
                p_test = data$p_test1,
                y_cal = data$y_cal2,
                alpha = alpha
            )
            pred_sets2 <- .getConformalPredSets(
                p_cal = data$p_cal1,
                p_test = data$p_test2,
                y_cal = data$y_cal1,
                alpha = alpha
            )
            pred_sets <- c(pred_sets1, pred_sets2)
        }
        # Order the prediction set
        pred_sets <- pred_sets[order(data$idx)]
    }

    ## Transform prediction with leaf nodes to prediction sets with the
    ## common ancestor if simplify_pred=TRUE
    if (simplify_pred) {
        pred_sets1 <- vapply(
            pred_sets,
            function(x) returnCommonAncestor(x, onto),
            character(1)
        )
        ## Check for ramification. If there is a ramification in the ontology
        ## and the children of the common ancestor include also labels not
        ## in the prediction set, don't return the common ancestor
        for (i in seq_len(length(pred_sets1))) {
            if (length(.children(pred_sets1[[i]], onto)) ==
                length(pred_sets[[i]])) {
                pred_sets[[i]] <- pred_sets1[[i]]
            }
        }
    }

    ## if not specified, return a sc object if the input was a sc object,
    ## a matrix if the input was a matrix
    if (is.null(return_sc) & sc) {
        return_sc <- TRUE
    }
    if (is.null(return_sc) & !sc) {
        return_sc <- FALSE
    }
    if (return_sc) {
        colData(x_query)[[pr_name]] <- pred_sets
        return(x_query)
    }
    if (!return_sc) {
        return(pred_sets)
    }
}
