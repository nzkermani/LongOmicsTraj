setMethod(
  "lot_compare_estimators",
  signature(... = "ANY"),
  function(..., .names = NULL) {

    objs <- list(...)

    if (length(objs) < 2) {
      stop("Provide at least two LongOmicsTraj objects.", call. = FALSE)
    }

    # Validate inputs
    if (!all(vapply(objs, function(x) is(x, "LongOmicsTraj"), logical(1)))) {
      stop("All inputs must be LongOmicsTraj objects.", call. = FALSE)
    }

    # Assign estimator names
    if (is.null(.names)) {
      .names <- paste0("estimator_", seq_along(objs))
    }

    if (length(.names) != length(objs)) {
      stop(".names must match number of objects.", call. = FALSE)
    }

    # Extract OTS and tag estimator
    ots_list <- mapply(function(obj, nm) {

      df <- as.data.frame(obj@ots)

      if (nrow(df) == 0) {
        stop("One object has empty OTS. Run lot_map_to_ots() first.", call. = FALSE)
      }

      df$estimator <- nm
      df

    }, objs, .names, SIMPLIFY = FALSE)

    # Combine
    combined <- do.call(rbind, ots_list)

    # Return as DataFrame
    S4Vectors::DataFrame(combined)
  }
)
