#' Construct a LongOmicsTraj object
#'
#' @export
LongOmicsTraj <- function(
    mae = MultiAssayExperiment::MultiAssayExperiment(),
    vmm = S4Vectors::DataFrame(),
    deltas = S4Vectors::DataFrame(),
    thresholds = S4Vectors::DataFrame(),
    ots = S4Vectors::DataFrame(),
    relationships = S4Vectors::DataFrame(),
    results = list(),
    cache = list(),
    settings = list(),
    history = list()
) {
  methods::new(
    "LongOmicsTraj",
    mae = mae,
    vmm = vmm,
    deltas = deltas,
    thresholds = thresholds,
    ots = ots,
    relationships = relationships,
    results = results,
    cache = cache,
    settings = settings,
    history = history
  )
}

setValidity("LongOmicsTraj", function(object) {
  errors <- character()

  req_vmm <- c(
    "object_id", "object_type", "assay", "group",
    "visit", "estimated_value", "estimator"
  )

  if (nrow(object@vmm) > 0) {
    missing_vmm <- setdiff(req_vmm, colnames(object@vmm))
    if (length(missing_vmm) > 0) {
      errors <- c(
        errors,
        paste("vmm is missing columns:", paste(missing_vmm, collapse = ", "))
      )
    }
  }

  req_rel <- c(
    "child_id", "parent_id", "relationship",
    "weight", "assay", "method"
  )

  if (nrow(object@relationships) > 0) {
    missing_rel <- setdiff(req_rel, colnames(object@relationships))
    if (length(missing_rel) > 0) {
      errors <- c(
        errors,
        paste("relationships is missing columns:", paste(missing_rel, collapse = ", "))
      )
    }
  }

  if (length(errors) == 0) TRUE else errors
})

#' @export
setMethod("show", "LongOmicsTraj", function(object) {
  cat("LongOmicsTraj object\n")
  cat("-------------------\n")
  cat("Assays:         ", length(MultiAssayExperiment::experiments(object@mae)), "\n")
  cat("VMM rows:       ", nrow(object@vmm), "\n")
  cat("Delta rows:     ", nrow(object@deltas), "\n")
  cat("Threshold rows: ", nrow(object@thresholds), "\n")
  cat("OTS rows:       ", nrow(object@ots), "\n")
  cat("Relationships:  ", nrow(object@relationships), "\n")
  cat("Results:        ", length(object@results), "\n")
})
