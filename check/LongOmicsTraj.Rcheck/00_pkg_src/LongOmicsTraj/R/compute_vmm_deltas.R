#' Compute adjacent visit-to-visit deltas from a Visit Mean Matrix
#'
#' @param vmm A data.frame with columns object_id, object_type, group, visit,
#' estimated_value, estimator.
#' @param visits Character vector giving the ordered visits.
#'
#' @return A data.frame of adjacent deltas.
#'
#' @examples
#' visits <- c(
#'   "Baseline",
#'   "Week 4",
#'   "Week 12"
#' )
#'
#' vmm <- data.frame(
#'   object_id = rep(
#'     c("gene_1", "gene_2"),
#'     each = 3
#'   ),
#'   object_type = "gene",
#'   group = "treated",
#'   visit = rep(
#'     visits,
#'     times = 2
#'   ),
#'   estimated_value = c(
#'     1.0, 1.6, 1.2,
#'     2.0, 2.0, 2.5
#'   ),
#'   estimator = "empirical",
#'   stringsAsFactors = FALSE
#' )
#'
#' deltas <- compute_vmm_deltas(
#'   vmm = vmm,
#'   visits = visits
#' )
#'
#' deltas[
#'   ,
#'   c(
#'     "object_id",
#'     "transition",
#'     "delta_hat"
#'   )
#' ]
#'
#' @export
compute_vmm_deltas <- function(vmm, visits) {

  required <- c(
    "object_id",
    "object_type",
    "group",
    "visit",
    "estimated_value",
    "estimator"
  )

  missing_cols <- setdiff(required, colnames(vmm))

  if (length(missing_cols) > 0) {
    stop(
      "vmm is missing required columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  if (length(visits) < 2) {
    stop("At least two ordered visits are required.", call. = FALSE)
  }

  vmm <- vmm[vmm$visit %in% visits, , drop = FALSE]
  vmm$visit <- factor(vmm$visit, levels = visits, ordered = TRUE)

  vmm <- vmm[order(
    vmm$object_id,
    vmm$object_type,
    vmm$group,
    vmm$visit
  ), , drop = FALSE]

  split_key <- paste(
    vmm$object_id,
    vmm$object_type,
    vmm$group,
    vmm$estimator,
    sep = "___"
  )

  vmm_split <- split(vmm, split_key)

  delta_list <- lapply(vmm_split, function(df) {

    df <- df[order(df$visit), , drop = FALSE]

    if (nrow(df) != length(visits)) {
      return(NULL)
    }

    data.frame(
      object_id = df$object_id[1],
      object_type = df$object_type[1],
      group = df$group[1],
      estimator = df$estimator[1],
      transition_index = seq_len(length(visits) - 1),
      from_visit = head(visits, -1),
      to_visit = tail(visits, -1),
      transition = paste(head(visits, -1), tail(visits, -1), sep = "_to_"),
      delta_hat = diff(df$estimated_value),
      stringsAsFactors = FALSE
    )
  })

  deltas <- do.call(rbind, delta_list)

  if (is.null(deltas)) {
    deltas <- data.frame()
  }

  rownames(deltas) <- NULL
  deltas
}
