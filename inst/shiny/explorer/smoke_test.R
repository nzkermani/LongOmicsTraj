# Run from the repository root after installing LongOmicsTraj.
# This smoke test does not require a shinyapps.io account.
stopifnot(requireNamespace("LongOmicsTraj", quietly = TRUE))
stopifnot(requireNamespace("shiny", quietly = TRUE))
stopifnot(requireNamespace("DT", quietly = TRUE))
parse(file = "inst/shiny/explorer/app.R")
data("lot_glucold_ots", package = "LongOmicsTraj")
x <- as.data.frame(LongOmicsTraj::lot_ots(lot_glucold_ots))
stopifnot(nrow(x) > 0L)
stopifnot(all(c("object_id", "group", "assay", "estimator", "topology_label") %in% names(x)))
one <- x[1L, ]
ans <- LongOmicsTraj::lot_query(
  lot_glucold_ots,
  assay = as.character(one$assay),
  estimator = as.character(one$estimator),
  conditions = stats::setNames(list(as.character(one$topology_label)), as.character(one$group))
)
stopifnot(nrow(ans) > 0L, as.character(one$object_id) %in% as.character(ans$object_id))
cat("Shiny Explorer syntax and core OTS query smoke test passed.\n")
