#' Summarise a LongOmicsTraj object
#'
#' Summarises assays, study design, available estimators and Omics
#' Trajectory Signature distributions stored in a `LongOmicsTraj` object.
#'
#' @param object A `LongOmicsTraj` object.
#' @param estimator Optional estimator used to filter OTS results.
#' @param top_n Number of most common OTS classes to retain per group.
#'
#' @return An object of class `lot_summary`. The object is a list containing
#'   dataset dimensions, study design information, result availability and
#'   OTS summaries.
#'
#' @export
#'
#' @examples
#' data("lot_glucold_ots")
#'
#' result <- lot_summary(
#'   lot_glucold_ots,
#'   estimator = "masigpro"
#' )
#'
#' result$assay_summary
#' result$top_ots
lot_summary <- function(
    object,
    estimator = NULL,
    top_n = 5L
) {
  
  # ---------------------------------------------------------------------------
  # Validate inputs
  # ---------------------------------------------------------------------------
  
  if (!methods::is(object, "LongOmicsTraj")) {
    stop(
      "`object` must be a LongOmicsTraj object.",
      call. = FALSE
    )
  }
  
  if (
    length(top_n) != 1L ||
    is.na(top_n) ||
    !is.numeric(top_n) ||
    top_n < 1L ||
    top_n != as.integer(top_n)
  ) {
    stop(
      "`top_n` must be one positive integer.",
      call. = FALSE
    )
  }
  
  top_n <- as.integer(top_n)
  
  # ---------------------------------------------------------------------------
  # Extract assays
  # ---------------------------------------------------------------------------
  
  experiments <- MultiAssayExperiment::experiments(
    object@mae
  )
  
  assay_names <- names(experiments)
  
  if (is.null(assay_names)) {
    assay_names <- character(0)
  }
  
  assay_summary <- data.frame(
    assay = assay_names,
    features = vapply(
      experiments,
      nrow,
      integer(1)
    ),
    samples = vapply(
      experiments,
      ncol,
      integer(1)
    ),
    stringsAsFactors = FALSE
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
  
  # ---------------------------------------------------------------------------
  # Extract study metadata
  # ---------------------------------------------------------------------------
  
  metadata <- NULL
  
  if (length(experiments) > 0L) {
    metadata <- as.data.frame(
      SummarizedExperiment::colData(
        experiments[[1L]]
      )
    )
  }
  
  subjects <- character(0)
  visits <- character(0)
  groups <- character(0)
  
  if (!is.null(metadata)) {
    
    if ("subject_id" %in% colnames(metadata)) {
      subjects <- unique(
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
        
        visits <- levels(
          droplevels(metadata$visit)
        )
        
        visits <- visits[
          visits %in% observed_visits
        ]
        
      } else {
        
        visits <- unique(
          as.character(
            metadata$visit[
              !is.na(metadata$visit)
            ]
          )
        )
      }
    }
    
    if ("group" %in% colnames(metadata)) {
      groups <- unique(
        as.character(
          metadata$group[
            !is.na(metadata$group)
          ]
        )
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Summarise OTS results
  # ---------------------------------------------------------------------------
  
  ots_data <- as.data.frame(
    object@ots
  )
  
  estimators <- character(0)
  ots_distribution <- data.frame()
  top_ots <- data.frame()
  
  if (nrow(ots_data) > 0L) {
    
    required_columns <- c(
      "group",
      "topology_label"
    )
    
    missing_columns <- setdiff(
      required_columns,
      colnames(ots_data)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "The OTS results are missing: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        call. = FALSE
      )
    }
    
    if ("estimator" %in% colnames(ots_data)) {
      
      estimators <- sort(
        unique(
          as.character(
            ots_data$estimator[
              !is.na(ots_data$estimator)
            ]
          )
        )
      )
      
      if (!is.null(estimator)) {
        
        if (
          !is.character(estimator) ||
          length(estimator) != 1L ||
          is.na(estimator) ||
          !nzchar(estimator)
        ) {
          stop(
            "`estimator` must be one non-empty character string.",
            call. = FALSE
          )
        }
        
        if (!estimator %in% estimators) {
          stop(
            "Estimator `",
            estimator,
            "` was not found. Available estimators: ",
            paste(
              estimators,
              collapse = ", "
            ),
            call. = FALSE
          )
        }
        
        ots_data <- ots_data[
          ots_data$estimator == estimator,
          ,
          drop = FALSE
        ]
      }
      
    } else if (!is.null(estimator)) {
      
      stop(
        "The OTS results do not contain an `estimator` column.",
        call. = FALSE
      )
    }
    
    ots_data <- ots_data[
      !is.na(ots_data$group) &
        !is.na(ots_data$topology_label),
      ,
      drop = FALSE
    ]
    
    if (nrow(ots_data) > 0L) {
      
      ots_distribution <- ots_data |>
        dplyr::count(
          .data$group,
          .data$topology_label,
          name = "count"
        ) |>
        dplyr::group_by(
          .data$group
        ) |>
        dplyr::mutate(
          proportion = .data$count /
            sum(.data$count)
        ) |>
        dplyr::ungroup() |>
        dplyr::arrange(
          .data$group,
          dplyr::desc(.data$count),
          .data$topology_label
        )
      

        top_ots <- ots_distribution |>
        dplyr::group_by(
          .data$group
        ) |>
        dplyr::slice_head(
          n = top_n
        ) |>
        dplyr::ungroup() 
    }
  }
  
  # ---------------------------------------------------------------------------
  # Check which analysis results are available
  # ---------------------------------------------------------------------------
  
  result_slots <- c(
    trajectories = "vmm",
    deltas = "deltas",
    thresholds = "thresholds",
    ots = "ots"
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
  
  # ---------------------------------------------------------------------------
  # Build the summary object
  # ---------------------------------------------------------------------------
  
  out <- list(
    assay_summary = assay_summary,
    n_assays = length(experiments),
    n_features = sum(
      assay_summary$features,
      na.rm = TRUE
    ),
    n_samples = length(unique_samples),
    n_subjects = length(subjects),
    subjects = subjects,
    visits = visits,
    groups = groups,
    estimators = estimators,
    selected_estimator = estimator,
    result_status = result_status,
    n_ots_classes = if (nrow(ots_distribution) > 0L) {
      dplyr::n_distinct(
        ots_distribution$topology_label
      )
    } else {
      0L
    },
    ots_distribution = ots_distribution,
    top_ots = top_ots
  )
  
  class(out) <- "lot_summary"
  
  out
}


#' Print a LongOmicsTraj summary
#'
#' @param x A `lot_summary` object.
#' @param ... Additional arguments, currently ignored.
#'
#' @return The summary object, returned invisibly.
#'
#' @export
print.lot_summary <- function(
    x,
    ...
) {
  
  cat("\n")
  cat("LongOmicsTraj summary\n")
  cat(
    strrep("\u2500", 60),
    "\n",
    sep = ""
  )
  
  # ---------------------------------------------------------------------------
  # Dataset overview
  # ---------------------------------------------------------------------------
  
  cat("\nDataset\n")
  
  cat(
    sprintf(
      "  Assays      : %d\n",
      x$n_assays
    )
  )
  
  cat(
    sprintf(
      "  Features    : %d\n",
      x$n_features
    )
  )
  
  cat(
    sprintf(
      "  Samples     : %d\n",
      x$n_samples
    )
  )
  
  cat(
    sprintf(
      "  Subjects    : %d\n",
      x$n_subjects
    )
  )
  
  # ---------------------------------------------------------------------------
  # Assays
  # ---------------------------------------------------------------------------
  
  cat("\nAssays\n")
  
  if (nrow(x$assay_summary) == 0L) {
    
    cat("  None\n")
    
  } else {
    
    for (i in seq_len(nrow(x$assay_summary))) {
      
      cat(
        sprintf(
          "  %-18s %5d features \u00d7 %4d samples\n",
          x$assay_summary$assay[[i]],
          x$assay_summary$features[[i]],
          x$assay_summary$samples[[i]]
        )
      )
    }
  }
  
  # ---------------------------------------------------------------------------
  # Study design
  # ---------------------------------------------------------------------------
  
  cat("\nStudy design\n")
  
  if (length(x$visits) > 0L) {
    
    cat(
      sprintf(
        "  Visits      : %s\n",
        paste(
          x$visits,
          collapse = " \u2192 "
        )
      )
    )
    
  } else {
    
    cat("  Visits      : not available\n")
  }
  
  if (length(x$groups) > 0L) {
    
    cat(
      sprintf(
        "  Groups      : %d\n",
        length(x$groups)
      )
    )
    
    cat(
      sprintf(
        "                 %s\n",
        paste(
          x$groups,
          collapse = ", "
        )
      )
    )
    
  } else {
    
    cat("  Groups      : not available\n")
  }
  
  # ---------------------------------------------------------------------------
  # Analysis status
  # ---------------------------------------------------------------------------
  
  cat("\nAnalysis status\n")
  
  status_symbol <- function(value) {
    if (isTRUE(value)) {
      "\u2713"
    } else {
      "\u2013"
    }
  }
  
  for (nm in names(x$result_status)) {
    
    cat(
      sprintf(
        "  %-12s : %s\n",
        nm,
        status_symbol(
          x$result_status[[nm]]
        )
      )
    )
  }
  
  # ---------------------------------------------------------------------------
  # Estimators
  # ---------------------------------------------------------------------------
  
  cat("\nEstimators\n")
  
  if (length(x$estimators) > 0L) {
    
    cat(
      sprintf(
        "  %s\n",
        paste(
          x$estimators,
          collapse = ", "
        )
      )
    )
    
  } else {
    
    cat("  None available\n")
  }
  
  # ---------------------------------------------------------------------------
  # OTS summary
  # ---------------------------------------------------------------------------
  
  cat("\nOmics Trajectory Signatures\n")
  
  cat(
    sprintf(
      "  Unique trajectory classes : %d\n",
      x$n_ots_classes
    )
  )
  
  if (nrow(x$top_ots) > 0L) {
    
    display_data <- x$top_ots |>
      dplyr::mutate(
        percentage = sprintf(
          "%.1f%%",
          100 * .data$proportion
        )
      )
    
    for (group_name in unique(display_data$group)) {
      
      cat("\n")
      
      cat(
        sprintf(
          "  %s\n",
          group_name
        )
      )
      
      group_data <- display_data[
        display_data$group == group_name,
        ,
        drop = FALSE
      ]
      
      for (i in seq_len(nrow(group_data))) {
        
        cat(
          sprintf(
            "    %-18s %5d (%6s)\n",
            group_data$topology_label[[i]],
            group_data$count[[i]],
            group_data$percentage[[i]]
          )
        )
      }
    }
    
  } else {
    
    cat("  No OTS results available.\n")
  }
  
  cat("\n")
  
  invisible(x)
}
