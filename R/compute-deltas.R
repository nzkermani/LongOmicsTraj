#' Compute deltas between visits
#' @export
setMethod(
  "lot_compute_deltas",
  "LongOmicsTraj",
  function(
    object,
    assay,
    visits,
    overwrite = FALSE
  ) {

    if (!overwrite && nrow(object@deltas) > 0) {
      stop("Deltas already exist. Use overwrite=TRUE.", call. = FALSE)
    }

    vmm <- as.data.frame(object@vmm)

    if (nrow(vmm) == 0) {
      stop("VMM is empty. Run a lot_from_* method first.", call. = FALSE)
    }

    # Ensure visit ordering
    vmm$visit <- factor(vmm$visit, levels = visits, ordered = TRUE)

    # Order properly
    vmm <- vmm[order(vmm$object_id, vmm$group, vmm$visit), ]

    # Compute deltas
    delta_list <- split(vmm, list(vmm$object_id, vmm$group), drop = TRUE)

    res <- lapply(delta_list, function(df) {

      if (nrow(df) < 2) return(NULL)

      df <- df[order(df$visit), ]

      data.frame(
        object_id = df$object_id[-1],
        object_type = df$object_type[-1],
        assay = df$assay[-1],
        group = df$group[-1],
        estimator = df$estimator[-1],
        transition_index = seq_len(nrow(df) - 1),
        from_visit = df$visit[-nrow(df)],
        to_visit = df$visit[-1],
        transition = paste(df$visit[-nrow(df)], df$visit[-1], sep = "_"),
        delta_hat = diff(df$estimated_value),
        stringsAsFactors = FALSE
      )
    })

    deltas <- do.call(rbind, res)

    object@deltas <- S4Vectors::DataFrame(deltas)

    object
  }
)
