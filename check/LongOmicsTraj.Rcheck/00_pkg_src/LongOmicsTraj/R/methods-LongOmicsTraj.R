#' Display a LongOmicsTraj object
#'
#' Prints a concise summary of the assays, study design and analysis
#' results stored in a `LongOmicsTraj` object.
#'
#' @param object A `LongOmicsTraj` object.
#'
#' @return The object is returned invisibly.
#'
#' @examples
#' data("lot_glucold")
#'
#' ## Printing the object invokes the show() method.
#' lot_glucold
#'
#' ## The object is returned invisibly.
#' invisible(lot_glucold)
#'
#' @export
methods::setMethod(
  "show",
  signature = "LongOmicsTraj",
  function(object) {
    
    # -------------------------------------------------------------------------
    # Extract experiments
    # -------------------------------------------------------------------------
    
    mae <- object@mae
    
    experiments <- MultiAssayExperiment::experiments(
      mae
    )
    
    assay_names <- names(experiments)
    
    if (is.null(assay_names)) {
      assay_names <- character(0)
    }
    
    # -------------------------------------------------------------------------
    # Summarise dimensions across assays
    # -------------------------------------------------------------------------
    
    feature_counts <- vapply(
      experiments,
      nrow,
      integer(1)
    )
    
    sample_counts <- vapply(
      experiments,
      ncol,
      integer(1)
    )
    
    total_features <- sum(
      feature_counts,
      na.rm = TRUE
    )
    
    unique_samples <- unique(
      unlist(
        lapply(
          experiments,
          colnames
        ),
        use.names = FALSE
      )
    )
    
    total_samples <- length(
      unique_samples
    )
    
    # -------------------------------------------------------------------------
    # Extract metadata from the first assay
    # -------------------------------------------------------------------------
    
    metadata <- NULL
    
    if (length(experiments) > 0L) {
      metadata <- as.data.frame(
        SummarizedExperiment::colData(
          experiments[[1L]]
        )
      )
    }
    
    subject_values <- character(0)
    visit_values <- character(0)
    group_values <- character(0)
    
    if (!is.null(metadata)) {
      
      if ("subject_id" %in% colnames(metadata)) {
        subject_values <- unique(
          as.character(
            metadata$subject_id[
              !is.na(metadata$subject_id)
            ]
          )
        )
      }
      
      if ("visit" %in% colnames(metadata)) {
        
        if (is.factor(metadata$visit)) {
          
          observed_visits <- unique(
            as.character(metadata$visit)
          )
          
          visit_values <- levels(
            droplevels(metadata$visit)
          )
          
          visit_values <- visit_values[
            visit_values %in% observed_visits
          ]
          
        } else {
          
          visit_values <- unique(
            as.character(
              metadata$visit[
                !is.na(metadata$visit)
              ]
            )
          )
        }
      }
      
      if ("group" %in% colnames(metadata)) {
        group_values <- unique(
          as.character(
            metadata$group[
              !is.na(metadata$group)
            ]
          )
        )
      }
    }
    
    # -------------------------------------------------------------------------
    # Identify available estimators
    # -------------------------------------------------------------------------
    
    estimator_values <- character(0)
    
    if (
      "ots" %in% methods::slotNames(object) &&
      nrow(as.data.frame(object@ots)) > 0L
    ) {
      
      ots_data <- as.data.frame(
        object@ots
      )
      
      if ("estimator" %in% colnames(ots_data)) {
        estimator_values <- sort(
          unique(
            as.character(
              ots_data$estimator[
                !is.na(ots_data$estimator)
              ]
            )
          )
        )
      }
    }
    
    # -------------------------------------------------------------------------
    # Check which result slots contain data
    # -------------------------------------------------------------------------
    
    result_slots <- c(
      trajectories = "vmm",
      deltas = "deltas",
      thresholds = "thresholds",
      OTS = "ots"
    )
    
    has_content <- function(slot_name) {
      
      if (!slot_name %in% methods::slotNames(object)) {
        return(FALSE)
      }
      
      value <- methods::slot(
        object,
        slot_name
      )
      
      if (is.null(value)) {
        return(FALSE)
      }
      
      if (
        methods::is(value, "DataFrame") ||
        is.data.frame(value) ||
        is.matrix(value)
      ) {
        return(nrow(value) > 0L)
      }
      
      if (is.list(value)) {
        return(length(value) > 0L)
      }
      
      length(value) > 0L
    }
    
    result_status <- vapply(
      result_slots,
      has_content,
      logical(1)
    )
    
    status_symbol <- ifelse(
      result_status,
      "\u2713",
      "\u2013"
    )
    
    # -------------------------------------------------------------------------
    # Print summary
    # -------------------------------------------------------------------------
    
    cat(
      "\nLongOmicsTraj object\n"
    )
    
    cat(
      paste0(
        strrep("\u2500", 44),
        "\n"
      )
    )
    
    cat("\nAssays\n")
    
    if (length(assay_names) == 0L) {
      
      cat("  None\n")
      
    } else {
      
      for (i in seq_along(assay_names)) {
        cat(
          "  ",
          assay_names[[i]],
          ": ",
          feature_counts[[i]],
          " features \u00d7 ",
          sample_counts[[i]],
          " samples\n",
          sep = ""
        )
      }
    }
    
    cat("\nDimensions\n")
    cat(
      "  Features : ",
      total_features,
      "\n",
      sep = ""
    )
    cat(
      "  Samples  : ",
      total_samples,
      "\n",
      sep = ""
    )
    
    if (length(subject_values) > 0L) {
      cat(
        "  Subjects : ",
        length(subject_values),
        "\n",
        sep = ""
      )
    }
    
    cat("\nStudy design\n")
    
    if (length(visit_values) > 0L) {
      cat(
        "  Visits : ",
        paste(
          visit_values,
          collapse = " \u2192 "
        ),
        "\n",
        sep = ""
      )
    } else {
      cat("  Visits : not available\n")
    }
    
    if (length(group_values) > 0L) {
      cat(
        "  Groups : ",
        length(group_values),
        "\n",
        sep = ""
      )
      
      cat(
        "           ",
        paste(
          group_values,
          collapse = ", "
        ),
        "\n",
        sep = ""
      )
    } else {
      cat("  Groups : not available\n")
    }
    
    cat("\nComputed results\n")
    
    for (i in seq_along(result_slots)) {
      cat(
        "  ",
        sprintf(
          "%-12s",
          names(result_slots)[[i]]
        ),
        ": ",
        status_symbol[[i]],
        "\n",
        sep = ""
      )
    }
    
    cat("\nEstimators\n")
    
    if (length(estimator_values) > 0L) {
      cat(
        "  ",
        paste(
          estimator_values,
          collapse = ", "
        ),
        "\n",
        sep = ""
      )
    } else {
      cat("  None recorded\n")
    }
    
    cat("\n")
    
    invisible(object)
  }
)