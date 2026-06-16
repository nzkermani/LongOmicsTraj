setMethod(
  "lot_consensus",
  signature(x = "DFrame"),
  function(
    x,
    min_agreement = 2,
    mode = c("strict", "majority")
  ) {

    mode <- match.arg(mode)

    df <- as.data.frame(x)

    required_cols <- c("object_id", "estimator", "topology_label")
    if (!all(required_cols %in% colnames(df))) {
      stop("Input must come from lot_compare_estimators()", call. = FALSE)
    }

    library(dplyr)
    library(tidyr)



    df_collapsed <- df %>%
      group_by(object_id, estimator, group) %>%
      summarise(
        topology_label = first(topology_label),
        .groups = "drop"
      )

    wide <- df_collapsed %>%
      pivot_wider(
        names_from = c(estimator, group),
        values_from = topology_label
      )

    topology_matrix <- wide[, -1, drop = FALSE]

    consensus_list <- apply(wide[, -1, drop = FALSE], 1, function(row) {

      vals <- as.character(row)
      vals <- vals[!is.na(vals)]

      if (length(vals) == 0) return(NA)

      tab <- table(vals)
      top <- names(tab)[which.max(tab)]
      count <- max(tab)

      if (mode == "majority" && count >= min_agreement) {
        return(top)
      } else {
        return(NA)
      }
    })

    result <- data.frame(
      object_id = wide$object_id,
      consensus_topology = consensus_list,
      stringsAsFactors = FALSE
    )

    S4Vectors::DataFrame(result)
  }
)
