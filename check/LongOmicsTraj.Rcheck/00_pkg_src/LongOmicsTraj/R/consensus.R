#' Calculate consensus trajectory assignments
#'
#' Combines OTS assignments stored in a `DFrame` and returns a consensus
#' topology label for each molecular feature.
#'
#' Multiple estimator and group-specific assignments are first collapsed to
#' one topology label per molecular feature, estimator and group. Consensus is
#' then calculated across the resulting assignments.
#'
#' @param x A `DFrame` containing OTS assignments from multiple estimators.
#' @param min_agreement Minimum number of assignments that must agree before a
#'   consensus topology is assigned.
#' @param mode Consensus rule. `"majority"` assigns the most frequent topology
#'   when its count is at least `min_agreement`. `"strict"` assigns a topology
#'   only when all available assignments agree and the number of assignments is
#'   at least `min_agreement`.
#'
#' @return A `DFrame` containing one row per molecular feature with the
#'   consensus topology, agreement count and number of available assignments.
#'
#'
#' @examples
#' assignments <- S4Vectors::DataFrame(
#'   object_id = rep(
#'     c("gene_1", "gene_2"),
#'     each = 3
#'   ),
#'   estimator = rep(
#'     c("empirical", "lmm", "gam"),
#'     times = 2
#'   ),
#'   group = "treated",
#'   topology_label = c(
#'     "up_down",
#'     "up_down",
#'     "flat_flat",
#'     "down_up",
#'     "down_up",
#'     "down_up"
#'   )
#' )
#'
#' consensus <- lot_consensus(
#'   assignments,
#'   min_agreement = 2L,
#'   mode = "majority"
#' )
#'
#' as.data.frame(consensus)
#'
#' @rdname lot_consensus
#' @export
methods::setMethod(
  "lot_consensus",
  signature(
    x = "DFrame"
  ),
  function(
    x,
    min_agreement = 2L,
    mode = c(
      "strict",
      "majority"
    )
  ) {
    
    mode <- match.arg(
      mode
    )
    
    if (
      length(min_agreement) != 1L ||
      !is.numeric(min_agreement) ||
      is.na(min_agreement) ||
      !is.finite(min_agreement) ||
      min_agreement < 1L ||
      min_agreement != as.integer(min_agreement)
    ) {
      stop(
        "`min_agreement` must be one positive integer.",
        call. = FALSE
      )
    }
    
    min_agreement <- as.integer(
      min_agreement
    )
    
    df <- as.data.frame(
      x
    )
    
    required_columns <- c(
      "object_id",
      "estimator",
      "group",
      "topology_label"
    )
    
    missing_columns <- setdiff(
      required_columns,
      colnames(df)
    )
    
    if (length(missing_columns) > 0L) {
      stop(
        "Input is missing required columns: ",
        paste(
          missing_columns,
          collapse = ", "
        ),
        ". The input should come from `lot_compare_estimators()`.",
        call. = FALSE
      )
    }
    
    df <- df[
      !is.na(df$object_id) &
        !is.na(df$estimator) &
        !is.na(df$group) &
        !is.na(df$topology_label),
      ,
      drop = FALSE
    ]
    
    if (nrow(df) == 0L) {
      stop(
        "No complete OTS assignments were available.",
        call. = FALSE
      )
    }
    
    df_collapsed <- df |>
      dplyr::group_by(
        .data$object_id,
        .data$estimator,
        .data$group
      ) |>
      dplyr::summarise(
        topology_label = dplyr::first(
          .data$topology_label
        ),
        .groups = "drop"
      )
    
    wide <- df_collapsed |>
      tidyr::pivot_wider(
        names_from = c(
          .data$estimator,
          .data$group
        ),
        values_from = .data$topology_label,
        names_sep = "__"
      )
    
    topology_matrix <- wide[
      ,
      setdiff(
        colnames(wide),
        "object_id"
      ),
      drop = FALSE
    ]
    
    consensus_result <- lapply(
      seq_len(
        nrow(topology_matrix)
      ),
      function(row_index) {
        
        values <- as.character(
          unlist(
            topology_matrix[
              row_index,
              ,
              drop = FALSE
            ],
            use.names = FALSE
          )
        )
        
        values <- values[
          !is.na(values) &
            nzchar(values)
        ]
        
        n_available <- length(
          values
        )
        
        if (n_available == 0L) {
          return(
            list(
              consensus_topology = NA_character_,
              agreement_count = 0L,
              n_available = 0L
            )
          )
        }
        
        frequency_table <- table(
          values
        )
        
        dominant_topology <- names(
          frequency_table
        )[
          which.max(
            frequency_table
          )
        ]
        
        agreement_count <- as.integer(
          max(
            frequency_table
          )
        )
        
        consensus_topology <- if (mode == "strict") {
          
          if (
            length(frequency_table) == 1L &&
            n_available >= min_agreement
          ) {
            dominant_topology
          } else {
            NA_character_
          }
          
        } else {
          
          if (agreement_count >= min_agreement) {
            dominant_topology
          } else {
            NA_character_
          }
        }
        
        list(
          consensus_topology = consensus_topology,
          agreement_count = agreement_count,
          n_available = n_available
        )
      }
    )
    
    result <- data.frame(
      object_id = as.character(
        wide$object_id
      ),
      consensus_topology = vapply(
        consensus_result,
        function(item) {
          item$consensus_topology
        },
        character(1)
      ),
      agreement_count = vapply(
        consensus_result,
        function(item) {
          item$agreement_count
        },
        integer(1)
      ),
      n_available = vapply(
        consensus_result,
        function(item) {
          item$n_available
        },
        integer(1)
      ),
      stringsAsFactors = FALSE
    )
    
    S4Vectors::DataFrame(
      result
    )
  }
)