#' scConform: Conformal prediction for cell-type annotation
#'
#' The \pkg{scConform} package provides conformal prediction methods for
#' uncertainty-aware cell-type annotation in single-cell RNA-seq data.
#' It constructs finite-sample valid prediction sets for cell types,
#' starting from model-based estimated class probabilities.
#'
#' The main entry point of the package is the function
#' [getPredictionSets()], which builds prediction sets for
#' query cells using a calibration dataset and a user-specified
#' miscoverage level \eqn{\alpha}.
#'
#' Two types of prediction sets are supported:
#' \itemize{
#'   \item *Classical conformal prediction sets*, which return subsets
#'   of cell-type labels without structural constraints;
#'   \item *Hierarchical conformal prediction sets*, which incorporate
#'   prior knowledge encoded in a cell ontology and return prediction sets
#'   that are consistent with the ontology structure.
#' }
#'
#' When hierarchical prediction is enabled, \pkg{scConform} uses a conformal
#' risk-control procedure to select a threshold parameter that guarantees
#' finite-sample coverage while producing interpretable, non-empty
#' ontology-aware prediction sets.
#'
#' The package is designed to work seamlessly with
#' `SingleCellExperiment` and `SpatialExperiment` objects, but
#' also supports matrix-based inputs of predicted class probabilities.
#'
#' The methodological framework implemented in \pkg{scConform} is described
#' in Corbetta et al. (2025).
#'
#' To get started, see the help page of [getPredictionSets()] and
#' the package vignette for a complete workflow.
#'
#' @references
#' Corbetta, D., Geistlinger L., Finos, L., and Risso, D. (2025).
#' *Conformal inference for cell type annotation with graph-structured
#' constraints*. arXiv preprint arXiv:2410.23786.
#'
#' @keywords internal
"_PACKAGE"
