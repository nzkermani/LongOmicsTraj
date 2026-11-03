#' Choose a depth-model restart by agreement across restarts
#'
#' Depth models such as GASTON are trained from several random restarts, and
#' the restart with the lowest training loss is not always the one with the
#' most faithful depth map. This function picks, without labels, the restart
#' whose depth coordinate agrees best with all other restarts.
#'
#' @param depths A numeric matrix with one row per spot (or cell) and one
#'   column per restart, holding each restart's depth coordinate.
#' @param absolute Logical; use absolute Spearman correlations, so restarts
#'   whose depth runs in the opposite direction still count as agreeing.
#'
#' @return A list with `best` (column index of the chosen restart),
#'   `agreement` (each restart's mean agreement with the others) and
#'   `correlation` (the restart-by-restart Spearman correlation matrix).
#'
#' @details A restart whose depth map is folded or misoriented agrees poorly
#' with the majority of restarts and is therefore not chosen. The rule needs
#' at least three restarts; ten or more are recommended.
#'
#' @examples
#' set.seed(1)
#' truth <- runif(200)
#' depths <- cbind(sapply(1:8, function(i) truth + rnorm(200, 0, 0.05)),
#'                 folded = abs(truth - 0.5))
#' lot_consensus_restart(depths)$best
#'
#' @export
lot_consensus_restart <- function(depths, absolute = TRUE) {
  depths <- as.matrix(depths)
  if (!is.numeric(depths)) {
    stop("`depths` must be a numeric matrix (spots x restarts).", call. = FALSE)
  }
  if (ncol(depths) < 3L) {
    stop("At least three restarts are needed.", call. = FALSE)
  }
  r <- stats::cor(depths, method = "spearman", use = "pairwise.complete.obs")
  if (absolute) r <- abs(r)
  diag(r) <- NA
  agreement <- rowMeans(r, na.rm = TRUE)
  if (is.null(colnames(depths))) names(agreement) <- seq_len(ncol(depths))
  list(best = unname(which.max(agreement)), agreement = agreement, correlation = r)
}


#' Build depth-window pseudobulk for LongOmicsTraj from spot-level data
#'
#' Converts a spot- or cell-level count matrix into pseudobulk profiles per
#' tissue section and ordered depth window, with the metadata that
#' [lot_create()] expects. Windows are treated as visits, sections as
#' subjects and donors or disease groups as groups.
#'
#' @param counts A numeric matrix (or sparse matrix) of raw counts with one
#'   row per gene and one column per spot. Column names must match
#'   `rownames(spots)` or `spots[[spot_col]]`.
#' @param spots A data frame with one row per spot.
#' @param section_col Column of `spots` identifying the tissue section.
#' @param group_col Column of `spots` identifying the donor or disease group.
#' @param window_col Column of `spots` holding a window label for each spot,
#'   for example annotated layers grouped into windows. Use this or
#'   `depth_col`.
#' @param depth_col Column of `spots` holding a continuous depth coordinate,
#'   for example a GASTON isodepth. Spots are cut into windows at `breaks`.
#' @param breaks Numeric cut points for `depth_col`, including the outer
#'   limits.
#' @param window_levels Window labels in depth order (superficial first). With
#'   `depth_col`, one label per interval of `breaks`.
#' @param spot_col Optional column of `spots` with spot identifiers; by
#'   default `rownames(spots)`.
#' @param min_spots Minimum spots for a section-by-window profile; smaller
#'   profiles are dropped with a message.
#' @param scale_factor Library-size scale for normalisation (counts per
#'   `scale_factor`).
#'
#' @return A list with `expression` (genes x profiles, mean log1p-normalised
#'   expression of the spots in each profile) and `metadata` (columns
#'   `sample_id`, `subject_id`, `visit`, `group`, `n_spots`), ready for
#'   `lot_create(expression = , metadata = , visit_levels = window_levels)`.
#'
#' @details Each spot is normalised to `scale_factor` total counts and
#' log1p-transformed, and the profile is the mean over the spots of that
#' section and window. Sections missing any window are kept; LongOmicsTraj
#' estimates group trajectories from the sections available at each window.
#'
#' @examples
#' set.seed(1)
#' counts <- matrix(rpois(50 * 120, 5), 50, 120,
#'                  dimnames = list(paste0("g", 1:50), paste0("s", 1:120)))
#' spots <- data.frame(section = rep(c("A", "B", "C", "D"), each = 30),
#'                     donor = rep(c("D1", "D2"), each = 60),
#'                     depth = runif(120), row.names = colnames(counts))
#' pb <- lot_spatial_pseudobulk(counts, spots, section_col = "section",
#'                              group_col = "donor", depth_col = "depth",
#'                              breaks = c(0, 1/3, 2/3, 1),
#'                              window_levels = c("upper", "middle", "deep"),
#'                              min_spots = 5)
#' head(pb$metadata)
#'
#' @export
lot_spatial_pseudobulk <- function(counts, spots, section_col, group_col,
                                   window_col = NULL, depth_col = NULL,
                                   breaks = NULL, window_levels = NULL,
                                   spot_col = NULL, min_spots = 20L,
                                   scale_factor = 1e4) {
  spots <- as.data.frame(spots)
  ids <- if (is.null(spot_col)) rownames(spots) else as.character(spots[[spot_col]])
  if (is.null(colnames(counts)) || !all(ids %in% colnames(counts))) {
    stop("Spot identifiers in `spots` must match `colnames(counts)`.", call. = FALSE)
  }
  for (cc in c(section_col, group_col)) {
    if (!cc %in% names(spots)) stop("Column '", cc, "' not found in `spots`.", call. = FALSE)
  }
  if (is.null(window_col) == is.null(depth_col)) {
    stop("Give exactly one of `window_col` or `depth_col`.", call. = FALSE)
  }
  if (!is.null(depth_col)) {
    if (is.null(breaks) || is.null(window_levels) || length(window_levels) != length(breaks) - 1L) {
      stop("With `depth_col`, give `breaks` and one `window_levels` label per interval.",
           call. = FALSE)
    }
    window <- as.character(cut(spots[[depth_col]], breaks = breaks, labels = window_levels,
                               include.lowest = TRUE))
  } else {
    window <- as.character(spots[[window_col]])
    if (is.null(window_levels)) {
      stop("Give `window_levels` (window labels in depth order).", call. = FALSE)
    }
  }
  keep <- !is.na(window) & window %in% window_levels
  section <- as.character(spots[[section_col]])
  group <- as.character(spots[[group_col]])
  dup <- tapply(group, section, function(x) length(unique(x)))
  if (any(dup > 1)) {
    stop("Each section must belong to a single group.", call. = FALSE)
  }
  key <- paste(section, window, sep = "\r")
  prof <- list(); meta <- list()
  for (k in unique(key[keep])) {
    idx <- which(keep & key == k)
    if (length(idx) < min_spots) {
      message("Dropping ", sub("\r", " / ", k, fixed = TRUE), ": ", length(idx), " spots")
      next
    }
    m <- as.matrix(counts[, ids[idx], drop = FALSE])
    lib <- colSums(m); lib[lib == 0] <- 1
    prof[[k]] <- rowMeans(log1p(sweep(m, 2, lib, "/") * scale_factor))
    meta[[k]] <- data.frame(subject_id = section[idx[1]], visit = window[idx[1]],
                            group = group[idx[1]], n_spots = length(idx),
                            stringsAsFactors = FALSE)
  }
  if (!length(prof)) stop("No section-by-window profile has `min_spots` spots.", call. = FALSE)
  expression <- do.call(cbind, prof)
  metadata <- do.call(rbind, meta)
  metadata$sample_id <- paste(metadata$subject_id, metadata$visit, sep = "_")
  colnames(expression) <- metadata$sample_id
  rownames(metadata) <- NULL
  metadata <- metadata[order(metadata$group, metadata$subject_id,
                             match(metadata$visit, window_levels)),
                       c("sample_id", "subject_id", "visit", "group", "n_spots")]
  rownames(metadata) <- NULL
  list(expression = expression[, metadata$sample_id, drop = FALSE], metadata = metadata)
}


#' Cross-group reproducibility of trajectory topologies, with a permutation test
#'
#' Runs the LongOmicsTraj workflow (empirical trajectories, changes between
#' consecutive visits, thresholds, OTS labels) and measures how often each
#' feature receives the same topology in every group. Significance comes from
#' re-running the workflow after shuffling visit labels within each subject.
#'
#' @param expression,metadata As for [lot_create()], with columns
#'   `sample_id`, `subject_id`, `visit` and `group` in `metadata`.
#' @param visits Visits (or spatial windows) in order.
#' @param center Passed to [lot_compute_deltas()]; `"relative"` is
#'   recommended for spatial windows with three or more groups.
#' @param method,k Threshold rule and multiplier, passed to
#'   [lot_compute_thresholds()].
#' @param nperm Number of permutations (0 for none).
#' @param seed Random seed for the permutations.
#' @param assay Assay name used internally.
#'
#' @return A list with `agreement` (fraction of features with the same
#'   topology in all groups), `kappa` (Fleiss kappa across groups),
#'   `null` (permutation agreements), `z` (standardised excess over the null),
#'   `p` (permutation p-value, (b + 1) / (nperm + 1)) and `labels` (features x
#'   groups topology labels).
#'
#' @details Agreement depends on the number of groups and features, so judge
#' it against the permutation null rather than in absolute terms. The
#' smallest attainable p-value is 1 / (nperm + 1).
#'
#' @export
lot_topology_agreement <- function(expression, metadata, visits,
                                   center = c("relative", "median", "mean", "none"),
                                   method = "hybrid", k = 0.25, nperm = 0L,
                                   seed = 1L, assay = "assay") {
  center <- match.arg(center)
  metadata <- as.data.frame(metadata)
  run <- function(md) {
    l <- lot_create(expression = expression, metadata = md, assay = assay,
                    visit_levels = visits, feature_type = "gene")
    l <- lot_from_empirical(l, assay = assay, visits = visits, group_col = "group",
                            visit_col = "visit", estimator = "empirical",
                            overwrite = TRUE, verbose = FALSE)
    l <- suppressWarnings(lot_compute_deltas(l, assay = assay, visits = visits,
                                             overwrite = TRUE, center = center))
    l <- lot_compute_thresholds(l, assay = assay, visits = visits, method = method,
                                k = k, overwrite = TRUE, verbose = FALSE)
    l <- lot_map_to_ots(l, assay = assay, visits = visits, overwrite = TRUE,
                        verbose = FALSE)
    o <- as.data.frame(l@ots)[, c("object_id", "group", "topology_label")]
    w <- stats::reshape(o, idvar = "object_id", timevar = "group", direction = "wide")
    w <- w[stats::complete.cases(w), ]
    lab <- as.matrix(w[, -1, drop = FALSE])
    rownames(lab) <- w$object_id
    colnames(lab) <- sub("^topology_label\\.", "", colnames(lab))
    lab
  }
  agree <- function(lab) mean(apply(lab, 1, function(r) length(unique(r)) == 1L))
  fleiss <- function(lab) {
    cats <- sort(unique(as.vector(lab))); n <- ncol(lab)
    nij <- t(apply(lab, 1, function(r) table(factor(r, levels = cats))))
    pbar <- mean((rowSums(nij^2) - n) / (n * (n - 1)))
    pj <- colSums(nij) / (nrow(nij) * n)
    (pbar - sum(pj^2)) / (1 - sum(pj^2))
  }
  lab <- run(metadata)
  if (ncol(lab) < 2L) stop("At least two groups are needed.", call. = FALSE)
  a <- agree(lab)
  null <- numeric(0)
  if (nperm > 0L) {
    null <- withr::with_seed(seed, vapply(seq_len(nperm), function(b) {
      md <- metadata
      for (s in unique(md$subject_id)) {
        i <- which(md$subject_id == s)
        md$visit[i] <- md$visit[i][sample.int(length(i))]
      }
      agree(run(md))
    }, numeric(1)))
  }
  list(agreement = a, kappa = fleiss(lab), null = null,
       z = if (length(null) > 1) (a - mean(null)) / stats::sd(null) else NA_real_,
       p = if (length(null)) (sum(null >= a) + 1) / (length(null) + 1) else NA_real_,
       labels = lab)
}
