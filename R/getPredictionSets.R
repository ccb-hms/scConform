#' @title Get conformal prediction sets
#' @description This function returns conformal prediction sets for the cell
#' type of cells in a query dataset. It implements two methods: standard split
#' conformal inference and a hierarchical conformal risk-control approach that
#' incorporates the cell ontology structure. Depending on the input and on the
#' value of `return_sc`, the output is either a list of prediction sets or a
#' `SingleCellExperiment`/`SpatialExperiment` object with prediction sets stored
#' in the `colData`.
#'
#' @param x_query query data for which we want to build prediction sets. This can
#' be either a `SingleCellExperiment` (or `SpatialExperiment`) object with the
#' estimated probabilities for
#' each cell type in the `colData`, or a named numeric matrix
#' with `n` rows and `K` columns,
#' where `n` is the number of cells and `K` is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param x_cal calibration data. This can be either a
#' `SingleCellExperiment` object with the estimated probabilities for each cell
#' type in the `colData`, or a named matrix of dimension `m x K`, where
#' `m` is the number of cells and `K` is the number of different
#' labels. The colnames of the matrix have to correspond to the cell labels.
#' @param y_cal a vector of length `m` with the true labels of the cells in
#' the calibration data.
#' @param onto An `igraph` object representing the considered section of the
#' cell ontology.
#' @param alpha Numeric value between 0 and 1 that indicates the allowed
#' miscoverage
#' @param lambdas a numeric vector of possible lambda values to be considered.
#' Necessary only when `follow_ontology=TRUE`.
#' @param follow_ontology Logical. If `TRUE`, then the function returns
#' hierarchical
#' prediction sets that follow the cell ontology structure. If `FALSE`, it
#' returns classical conformal prediction sets. See Details.
#' @param resample Logical. If `TRUE`, the calibration dataset is resampled
#' according to the estimated relative frequencies of cell types in the query
#' data.
#' @param labels Character vector of labels of different considered cell types.
#' Necessary if
#' `onto=NULL`, otherwise they are set equal to the leaf nodes of the
#' provided graph.
#' @param return_sc Logical. Parameter the controls the output type. If
#' `TRUE`, the function returns a `SingleCellExperiment`.
#' If `FALSE`, the function returns a list. By default,
#' it is set to `TRUE` when `x_query` is a `SingleCellExperiment` (or
#' `SpatialExperiment`) object and to `FALSE` when `x_query` is a matrix.
#' @param pr_name Character string giving the name of the `colData` variable
#' in the returned
#' `SingleCellExperiment` object that will contain the prediction
#' sets. The default name is `pred.set`.
#' @param simplify Logical. If `TRUE`, the output will be the common
#' ancestor
#' of the labels inserted into the prediction set. If `FALSE` (default),
#' the output will be the set of the leaf labels.
#' @param method character string or function specifying how hierarchical
#' prediction sets are constructed when `follow_ontology=TRUE`.
#' If a character string, it must be one of:
#' \describe{
#'   \item{`"full"`}{the default hierarchical construction described in the
#'   Details section, which guarantees non-empty prediction sets;}
#'   \item{`"step"`}{a construction that includes all ancestors up to a
#'   fixed number of steps above the predicted class;}
#'   \item{`"rank"`}{a construction that ranks leaf nodes by predicted
#'   probability, accumulates probability until a threshold is reached, and
#'   then expands the lowest common ancestor of the selected leaves.}
#' }
#' Alternatively, `method` can be a user-defined function with signature
#' `function(lambda, pred, onto)`, returning a character vector of leaf
#' labels defining the prediction set.
#' @param BPPARAM BiocParallel instance for parallel computing. Default is
#' `SerialParam()`.
#' @return
#' \describe{
#'   \item{If `return_sc = TRUE`:}{A `SingleCellExperiment` or
#'   `SpatialExperiment` object with the prediction sets stored in the
#'   `colData`. The name of the corresponding variable is given by `pr_name`.}
#'   \item{If `return_sc = FALSE`:}{A list of length equal to the number of
#'   cells in the query data. Each element contains the prediction set for one
#'   cell.}
#' }
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
#'   , obtain the *conformal scores*, \eqn{s_i=1-\hat{f}(X_i)_{Y_i},
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
#' \eqn{\lambda}. To choose \eqn{\lambda}, we follow eq. (4) in Angelopoulos et
#' al. (2023), considering the miscoverage as loss function. In this way, it is
#' still guaranteed that
#' \deqn{P(Y_{n+1}\notin C_\lambda (X_{n+1})) \leq \alpha.}
#' The construction described above corresponds to the default choice
#' `method = "full"`. Other values of `method` implement alternative
#' nested prediction-set constructions that incorporate the ontology structure
#' in different ways. All methods are calibrated using the same conformal
#' risk-control procedure to select the threshold parameter \eqn{\lambda}.}
#' @references
#' Corbetta, D. et al. *Conformal inference for cell type annotation with
#' graph-structured constraints*. arXiv preprint arXiv:2410.23786.
#'
#' Angelopoulos, A. N. and Bates, S. *A gentle introduction to conformal
#' prediction and distribution-free uncertainty quantification*. arXiv preprint
#' arXiv:2107.07511.
#'
#' Angelopoulos, A. N. et al. *Conformal risk control*. arXiv preprint
#' arXiv:2208.02814.
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
#'     resample = FALSE,
#'     labels = cell_types,
#'     return_sc = FALSE
#' )
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom igraph V distances degree
#' @importFrom stats quantile
#' @importFrom BiocParallel SerialParam bplapply
#' @export

getPredictionSets <- function(
        x_query,
        x_cal,
        y_cal,
        onto = NULL,
        alpha = 0.1,
        lambdas = seq(0.001, 0.999, length.out = 100),
        follow_ontology = TRUE,
        resample = FALSE,
        labels = NULL,
        return_sc = NULL,
        pr_name = "pred.set",
        simplify = FALSE,
        method = "full",
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

    if (!follow_ontology & simplify) {
        stop("If follow_ontology=FALSE, please set simplify=FALSE")
    }

    # Validate method
    if (is.character(method)) {
        allowed_methods <- c("full", "step", "rank")
        if (length(method) != 1 || !method %in% allowed_methods) {
            stop(
                "If 'method' is a character, it must be one of: ",
                paste(allowed_methods, collapse = ", ")
            )
        }
    } else if (is.function(method)) {
        # ok: user-supplied prediction-set constructor
    } else {
        stop("'method' must be either a character string or a function")
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

    if (!resample) {
        if (follow_ontology) {
            pred_sets <- .getHierarchicalPredSets(
                p_cal = p_cal, p_test = p_query,
                y_cal = y_cal, onto = onto,
                alpha = alpha,
                lambdas = lambdas,
                BPPARAM = BPPARAM,
                method = method
            )
        } else {
            pred_sets <- .getConformalPredSets(
                p_cal = p_cal, p_test = p_query,
                y_cal = y_cal, alpha = alpha
            )
        }
    }

    if (resample) {
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
                BPPARAM = BPPARAM,
                method = method
            )
            pred_sets2 <- .getHierarchicalPredSets(
                p_cal = data$p_cal1,
                p_test = data$p_test2,
                y_cal = data$y_cal1,
                onto = onto,
                alpha = alpha,
                lambdas = lambdas,
                BPPARAM = BPPARAM,
                method = method
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
    ## common ancestor if simplify=TRUE
    if (simplify) {
        pred_sets1 <- vapply(
            pred_sets,
            function(x) getCommonAncestor(x, onto),
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
