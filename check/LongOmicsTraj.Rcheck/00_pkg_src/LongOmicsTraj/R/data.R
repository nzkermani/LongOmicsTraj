#' Cleaned longitudinal transcriptomics dataset from the GLUCOLD study
#'
#' A pre-processed `LongOmicsTraj` object containing an integrated
#' `MultiAssayExperiment` of 71 candidate and hypothesis-driven genes across
#' longitudinal milestones from the GLUCOLD study.
#'
#' @format A `LongOmicsTraj` object containing a `MultiAssayExperiment`
#' with one experiment named `"transcriptomics"`.
#'
#' The internal `colData` contains:
#'
#' \describe{
#'   \item{sample_id}{
#'     Character. Unique sample identifier.
#'   }
#'   \item{subject_id}{
#'     Character. Unique participant identifier.
#'   }
#'   \item{visit}{
#'     Factor. Chronological study visit.
#'   }
#'   \item{group}{
#'     Factor. Treatment-allocation group.
#'   }
#'   \item{time}{
#'     Numeric or factor representation of study time.
#'   }
#'   \item{treatment}{
#'     Factor. Original clinical treatment allocation.
#'   }
#' }
#'
#' @source NCBI Gene Expression Omnibus accession GSE36221.
#'
#' @usage data(lot_glucold)
"lot_glucold"
#' GLUCOLD clinical metadata
#'
#' Clinical and longitudinal sample metadata used to construct the example
#' `LongOmicsTraj` objects.
#'
#' @format A data frame containing sample identifiers, subject identifiers,
#' study visits, treatment groups, time and treatment allocation.
#'
#' @source NCBI Gene Expression Omnibus accession GSE36221.
#'
#' @usage data(lot_glucold_clinical)
"lot_glucold_clinical"


#' GLUCOLD expression matrix
#'
#' Pre-processed longitudinal transcriptomic measurements for 71 selected
#' candidate and hypothesis-driven genes from the GLUCOLD study.
#'
#' @format A numeric matrix with molecular features in rows and samples
#' in columns.
#'
#' @source NCBI Gene Expression Omnibus accession GSE36221.
#'
#' @usage data(lot_glucold_expression)
"lot_glucold_expression"


#' GLUCOLD LongOmicsTraj object with OTS results
#'
#' A processed `LongOmicsTraj` object containing model-estimated molecular
#' trajectories, between-visit changes, thresholds and Omics Trajectory
#' Signature assignments.
#'
#' @format A `LongOmicsTraj` object.
#'
#' @source Derived from NCBI Gene Expression Omnibus accession GSE36221.
#'
#' @usage data(lot_glucold_ots)
"lot_glucold_ots"