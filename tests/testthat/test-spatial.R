test_that("lot_consensus_restart rejects a folded restart", {
  set.seed(1)
  truth <- runif(300)
  depths <- cbind(sapply(1:6, function(i) truth + rnorm(300, 0, 0.05)),
                  folded = abs(truth - 0.5), flipped = -truth)
  res <- lot_consensus_restart(depths)
  expect_true(res$best %in% c(1:6, 8))           # a faithful (or flipped) restart, never the folded one
  expect_lt(res$agreement[["folded"]], min(res$agreement[1:6]))
  expect_error(lot_consensus_restart(depths[, 1:2]), "three restarts")
})

test_that("lot_spatial_pseudobulk builds lot_create-ready input", {
  set.seed(2)
  counts <- matrix(rpois(40 * 160, 4), 40, 160,
                   dimnames = list(paste0("g", 1:40), paste0("s", 1:160)))
  spots <- data.frame(section = rep(c("A", "B", "C", "D"), each = 40),
                      donor = rep(c("D1", "D2"), each = 80),
                      depth = rep(seq(0, 1, length.out = 40), 4),
                      row.names = colnames(counts))
  pb <- lot_spatial_pseudobulk(counts, spots, section_col = "section", group_col = "donor",
                               depth_col = "depth", breaks = c(0, 0.5, 1),
                               window_levels = c("upper", "deep"), min_spots = 5)
  expect_equal(nrow(pb$metadata), 8L)
  expect_equal(colnames(pb$expression), pb$metadata$sample_id)
  expect_true(all(pb$metadata$visit %in% c("upper", "deep")))
  expect_equal(sum(pb$metadata$n_spots), 160L)
  # normalisation: mean of log1p(CP10k) over the spots of one profile
  i <- rownames(spots)[spots$section == "A" & spots$depth <= 0.5]
  m <- counts[, i]; ref <- rowMeans(log1p(sweep(m, 2, colSums(m), "/") * 1e4))
  expect_equal(unname(pb$expression[, "A_upper"]), unname(ref))
  # window labels instead of a depth coordinate
  spots$layer <- ifelse(spots$depth <= 0.5, "upper", "deep")
  pb2 <- lot_spatial_pseudobulk(counts, spots, section_col = "section", group_col = "donor",
                                window_col = "layer", window_levels = c("upper", "deep"),
                                min_spots = 5)
  expect_equal(pb2$expression, pb$expression)
  expect_error(lot_spatial_pseudobulk(counts, spots, "section", "donor"), "exactly one")
  # small profiles are dropped; nothing left is an error
  spots$layer[rownames(spots) %in% paste0("s", 1:15)] <- "outside"
  expect_message(pb3 <- lot_spatial_pseudobulk(counts, spots, "section", "donor", window_col = "layer",
                                               window_levels = c("upper", "deep"), min_spots = 10),
                 "Dropping A / upper")
  expect_equal(nrow(pb3$metadata), 7L)
  expect_error(suppressMessages(lot_spatial_pseudobulk(counts, spots, "section", "donor", window_col = "layer",
                                      window_levels = c("upper", "deep"), min_spots = 30)),
               "No section-by-window")
})

test_that("lot_topology_agreement recovers shared topologies and a group-specific offset", {
  set.seed(3)
  G <- 120; visits <- c("upper", "deep", "WM")
  truth <- matrix(sample(c(-1, 0, 1), 2 * G, TRUE, c(0.2, 0.6, 0.2)), G)
  prof <- cbind(0, 0.5 * truth[, 1], 0.5 * (truth[, 1] + truth[, 2])) + 2
  meta <- expand.grid(window = visits, section = 1:3, donor = 1:3, stringsAsFactors = FALSE)
  E <- sapply(seq_len(nrow(meta)), function(i) {
    w <- match(meta$window[i], visits)
    prof[, w] + rnorm(G, 0, 0.1) + if (meta$donor[i] == 2 && w > 1) 0.4 else 0
  })
  rownames(E) <- paste0("g", 1:G)
  md <- data.frame(sample_id = paste0("d", meta$donor, "s", meta$section, meta$window),
                   subject_id = paste0("d", meta$donor, "s", meta$section),
                   visit = meta$window, group = paste0("D", meta$donor))
  colnames(E) <- md$sample_id
  none <- lot_topology_agreement(E, md, visits, center = "none")
  rel <- lot_topology_agreement(E, md, visits, center = "relative", nperm = 9)
  expect_gt(rel$agreement, none$agreement)
  expect_equal(ncol(rel$labels), 3L)
  expect_length(rel$null, 9L)
  expect_true(rel$p >= 0.1 && rel$p <= 1)
  expect_equal(rel$p, 0.1)                       # observed exceeds every permutation
})
