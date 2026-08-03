#' The LongOmicsTraj Class
#'
#' S4 container for supervised and unsupervised longitudinal omics topology analysis.
#'
#' @slot mae MultiAssayExperiment object storing raw/normalised multi-omics data.
#' @slot vmm Visit Mean Matrix as S4Vectors::DataFrame.
#' @slot deltas Adjacent visit-to-visit changes.
#' @slot thresholds Feature/object-specific topology thresholds.
#' @slot ots Overall Trajectory State labels.
#' @slot relationships Child-parent relationships, e.g. gene-to-cluster.
#' @slot results List of query, enrichment, diagnostic, or model outputs.
#' @slot cache List for reusable intermediate objects.
#' @slot settings List of analysis settings.
#' @slot history List recording processing steps.
#'
#' @import methods
#' @importFrom MultiAssayExperiment MultiAssayExperiment
#' @importFrom S4Vectors DataFrame
#' @export
setClass(
  "LongOmicsTraj",
  slots = list(
    mae           = "MultiAssayExperiment",
    vmm           = "DataFrame",
    deltas        = "DataFrame",
    thresholds    = "DataFrame",
    ots           = "DataFrame",
    relationships = "DataFrame",
    results       = "list",
    cache         = "list",
    settings      = "list",
    history       = "list"
  ),
  prototype = list(
    mae           = MultiAssayExperiment::MultiAssayExperiment(),
    vmm           = S4Vectors::DataFrame(),
    deltas        = S4Vectors::DataFrame(),
    thresholds    = S4Vectors::DataFrame(),
    ots           = S4Vectors::DataFrame(),
    relationships = S4Vectors::DataFrame(),
    results       = list(),
    cache         = list(),
    settings      = list(),
    history       = list()
  )
)
