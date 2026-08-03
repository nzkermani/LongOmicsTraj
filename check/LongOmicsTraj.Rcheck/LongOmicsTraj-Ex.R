pkgname <- "LongOmicsTraj"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
base::assign(".ExTimings", "LongOmicsTraj-Ex.timings", pos = 'CheckExEnv')
base::cat("name\tuser\tsystem\telapsed\n", file=base::get(".ExTimings", pos = 'CheckExEnv'))
base::assign(".format_ptime",
function(x) {
  if(!is.na(x[4L])) x[1L] <- x[1L] + x[4L]
  if(!is.na(x[5L])) x[2L] <- x[2L] + x[5L]
  options(OutDec = '.')
  format(x[1L:3L], digits = 7L)
},
pos = 'CheckExEnv')

### * </HEADER>
library('LongOmicsTraj')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
base::assign(".old_wd", base::getwd(), pos = 'CheckExEnv')
cleanEx()
nameEx("compute_vmm_deltas")
### * compute_vmm_deltas

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: compute_vmm_deltas
### Title: Compute adjacent visit-to-visit deltas from a Visit Mean Matrix
### Aliases: compute_vmm_deltas

### ** Examples

visits <- c(
  "Baseline",
  "Week 4",
  "Week 12"
)

vmm <- data.frame(
  object_id = rep(
    c("gene_1", "gene_2"),
    each = 3
  ),
  object_type = "gene",
  group = "treated",
  visit = rep(
    visits,
    times = 2
  ),
  estimated_value = c(
    1.0, 1.6, 1.2,
    2.0, 2.0, 2.5
  ),
  estimator = "empirical",
  stringsAsFactors = FALSE
)

deltas <- compute_vmm_deltas(
  vmm = vmm,
  visits = visits
)

deltas[
  ,
  c(
    "object_id",
    "transition",
    "delta_hat"
  )
]




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("compute_vmm_deltas", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_accessors")
### * lot_accessors

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_mae
### Title: Access components of a LongOmicsTraj object
### Aliases: lot_mae lot_vmm lot_deltas lot_thresholds lot_ots
###   lot_relationships lot_results lot_mae,LongOmicsTraj-method
###   lot_vmm,LongOmicsTraj-method lot_deltas,LongOmicsTraj-method
###   lot_thresholds,LongOmicsTraj-method lot_ots,LongOmicsTraj-method
###   lot_relationships,LongOmicsTraj-method
###   lot_results,LongOmicsTraj-method lot_accessors

### ** Examples

data("lot_glucold_ots")

## Extract the underlying experiment container.
mae <- lot_mae(lot_glucold_ots)
names(MultiAssayExperiment::experiments(mae))

## Inspect the visit-level molecular estimates.
vmm <- lot_vmm(lot_glucold_ots)
dim(vmm)
utils::head(as.data.frame(vmm))

## Examine the estimated changes between adjacent visits.
deltas <- lot_deltas(lot_glucold_ots)
dim(deltas)
utils::head(
  as.data.frame(deltas)[
    ,
    c(
      "object_id",
      "group",
      "transition",
      "delta_hat"
    )
  ]
)

## Inspect the thresholds used to classify molecular changes.
thresholds <- lot_thresholds(lot_glucold_ots)
utils::head(
  as.data.frame(thresholds)[
    ,
    c(
      "object_id",
      "threshold_method",
      "delta_g"
    )
  ]
)

## Inspect the resulting Omics Trajectory Signatures.
ots <- lot_ots(lot_glucold_ots)
table(
  as.data.frame(ots)$topology_label
)

## Relationships and additional results may be empty when no
## clustering analysis has been performed.
relationships <- lot_relationships(lot_glucold_ots)
nrow(relationships)

results <- lot_results(lot_glucold_ots)
names(results)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_accessors", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_cluster_purity")
### * lot_cluster_purity

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_cluster_purity
### Title: Calculate topology purity within clusters
### Aliases: lot_cluster_purity

### ** Examples

assignments <- data.frame(
  gene_id = paste0("G", 1:10),
  cluster = c(rep(1, 5), rep(2, 5)),
  group = "Withdrawal",
  topology_label = c(
    rep("up_down", 4),
    "flat_flat",
    rep("down_up", 5)
  )
)

result <- lot_cluster_purity(
  assignments,
  purity_threshold = 0.90
)

result$cluster_group_purity
result$overall



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_cluster_purity", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_cluster_topology_score")
### * lot_cluster_topology_score

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_cluster_topology_score
### Title: Calculate a topology-aware cluster score
### Aliases: lot_cluster_topology_score

### ** Examples

assignments <- data.frame(
  gene_id = paste0("G", 1:12),
  cluster = c(
    rep(1, 4),
    rep(2, 4),
    rep(3, 4)
  ),
  group = "Withdrawal",
  topology_label = c(
    rep("flat_flat", 4),
    rep("flat_flat", 3),
    "flat_up",
    rep("down_up", 4)
  )
)

purity <- lot_cluster_purity(
  assignments,
  purity_threshold = 0.90
)

score <- lot_cluster_topology_score(
  purity,
  k = 3,
  bic = 125.4,
  minimum_cluster_size = 4
)

score$group_summary
score$overall



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_cluster_topology_score", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_compare_estimators")
### * lot_compare_estimators

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_compare_estimators
### Title: Compare trajectory estimators
### Aliases: lot_compare_estimators
###   lot_compare_estimators,LongOmicsTraj-method

### ** Examples

data("lot_glucold_ots")

## For illustration, compare two stored trajectory analyses.
comparison <- lot_compare_estimators(
  lot_glucold_ots,
  lot_glucold_ots,
  .names = c(
    "analysis_1",
    "analysis_2"
  )
)

comparison <- as.data.frame(comparison)

table(
  comparison$estimator
)

utils::head(
  comparison[
    ,
    c(
      "object_id",
      "group",
      "topology_label",
      "estimator"
    )
  ]
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_compare_estimators", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_compute_deltas")
### * lot_compute_deltas

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_compute_deltas
### Title: Compute between-visit molecular changes
### Aliases: lot_compute_deltas lot_compute_deltas,LongOmicsTraj-method

### ** Examples

data("lot_glucold_ots")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

fitted <- lot_compute_deltas(
  object = lot_glucold_ots,
  assay = "transcriptomics",
  visits = visits,
  overwrite = TRUE
)

deltas <- as.data.frame(
  lot_deltas(fitted)
)

utils::head(
  deltas[
    ,
    c(
      "object_id",
      "group",
      "transition",
      "delta_hat"
    )
  ]
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_compute_deltas", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_compute_thresholds")
### * lot_compute_thresholds

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_compute_thresholds
### Title: Compute molecular-change thresholds
### Aliases: lot_compute_thresholds
###   lot_compute_thresholds,LongOmicsTraj-method

### ** Examples

data("lot_glucold_ots")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

fitted <- lot_compute_thresholds(
  object = lot_glucold_ots,
  assay = "transcriptomics",
  method = "hybrid",
  visits = visits,
  overwrite = TRUE
)

thresholds <- as.data.frame(
  lot_thresholds(fitted)
)

utils::head(
  thresholds[
    ,
    c(
      "object_id",
      "delta_g",
      "threshold_method"
    )
  ]
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_compute_thresholds", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_consensus")
### * lot_consensus

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_consensus
### Title: Calculate consensus trajectory assignments
### Aliases: lot_consensus lot_consensus,DFrame-method

### ** Examples

assignments <- S4Vectors::DataFrame(
  object_id = rep(
    c("gene_1", "gene_2"),
    each = 3
  ),
  estimator = rep(
    c("empirical", "lmm", "gam"),
    times = 2
  ),
  group = "treated",
  topology_label = c(
    "up_down",
    "up_down",
    "flat_flat",
    "down_up",
    "down_up",
    "down_up"
  )
)

consensus <- lot_consensus(
  assignments,
  min_agreement = 2L,
  mode = "majority"
)

as.data.frame(consensus)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_consensus", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_create")
### * lot_create

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_create
### Title: Create a LongOmicsTraj object
### Aliases: lot_create

### ** Examples

data("lot_glucold_expression")
data("lot_glucold_clinical")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

time_values <- c(
  "Baseline" = 0,
  "6 Months" = 6,
  "30 Months" = 30
)

lot <- lot_create(
  expression = lot_glucold_expression,
  metadata = lot_glucold_clinical,
  assay = "transcriptomics",
  visit_levels = visits,
  time_values = time_values
)

methods::is(
  lot,
  "LongOmicsTraj"
)

names(
  MultiAssayExperiment::experiments(
    lot_mae(lot)
  )
)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_create", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_from_flexmix_clusters")
### * lot_from_flexmix_clusters

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_from_flexmix_clusters
### Title: Cluster molecular features using FlexMix trajectory models
### Aliases: lot_from_flexmix_clusters

### ** Examples

## Not run: 
##D data("lot_glucold_ots")
##D 
##D lot_flexmix <- lot_from_flexmix_clusters(
##D   object = lot_glucold_ots,
##D   assay = "transcriptomics",
##D   group = "ICS 6 months then withdrawal",
##D   visits = c(
##D     "Baseline",
##D     "6 Months",
##D     "30 Months"
##D   ),
##D   trajectory_source = "estimated",
##D   trajectory_estimator = "masigpro",
##D   k = 2:5,
##D   purity_threshold = 0.90,
##D   required_passing_fraction = 0.80,
##D   min_cluster_size = 5
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_from_flexmix_clusters", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_map_to_ots")
### * lot_map_to_ots

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_map_to_ots
### Title: Map molecular changes to Omics Trajectory Signatures
### Aliases: lot_map_to_ots lot_map_to_ots,LongOmicsTraj-method

### ** Examples

data("lot_glucold_ots")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

fitted <- lot_map_to_ots(
  object = lot_glucold_ots,
  visits = visits,
  overwrite = TRUE
)

utils::head(
  as.data.frame(
    lot_ots(fitted)
  )
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_map_to_ots", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_cluster_selection")
### * lot_plot_cluster_selection

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_cluster_selection
### Title: Plot topology-aware cluster selection
### Aliases: lot_plot_cluster_selection

### ** Examples

data("lot_glucold_ots")

example_object <- lot_glucold_ots

example_results <- lot_results(
  example_object
)

example_results$flexmix_model_comparison <- S4Vectors::DataFrame(
  k = 2:4,
  BIC = c(540, 505, 512),
  weighted_purity = c(0.76, 0.91, 0.88),
  passing_fraction = c(0.50, 1.00, 0.75),
  minimum_cluster_size = c(18, 12, 6),
  minimum_purity = c(0.68, 0.87, 0.72),
  mean_purity = c(0.78, 0.92, 0.86)
)

example_results$flexmix_selected_k <- 3L

methods::slot(
  example_object,
  "results"
) <- example_results

plot <- lot_plot_cluster_selection(
  object = example_object
)

plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_cluster_selection", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_cluster_topology")
### * lot_plot_cluster_topology

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_cluster_topology
### Title: Plot trajectory topology profiles for fitted clusters
### Aliases: lot_plot_cluster_topology

### ** Examples

## Not run: 
##D lot_plot_cluster_topology(
##D   lot_flexmix,
##D   show_features = TRUE,
##D   sort_by = "purity"
##D )
## End(Not run)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_cluster_topology", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_cluster_trajectory_arms")
### * lot_plot_cluster_trajectory_arms

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_cluster_trajectory_arms
### Title: Plot trajectory arms within molecular-feature clusters
### Aliases: lot_plot_cluster_trajectory_arms

### ** Examples

data("lot_glucold_ots")

example_object <- lot_glucold_ots

example_results <- lot_results(
  example_object
)

example_results$cluster_ots_distribution <-
  S4Vectors::DataFrame(
    cluster_id = c(
      "cluster_1",
      "cluster_1",
      "cluster_2",
      "cluster_2"
    ),
    group = "treated",
    topology_label = c(
      "up_down",
      "flat_flat",
      "down_up",
      "flat_down"
    ),
    count = c(
      8L,
      2L,
      7L,
      3L
    ),
    proportion = c(
      0.8,
      0.2,
      0.7,
      0.3
    )
  )

methods::slot(
  example_object,
  "results"
) <- example_results

plot <- lot_plot_cluster_trajectory_arms(
  object = example_object,
  groups = "treated"
)

plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_cluster_trajectory_arms", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_ots_distribution")
### * lot_plot_ots_distribution

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_ots_distribution
### Title: Plot the distribution of Omics Trajectory Signatures
### Aliases: lot_plot_ots_distribution

### ** Examples

data("lot_glucold_ots")

plot <- lot_plot_ots_distribution(
  object = lot_glucold_ots,
  estimator = "masigpro",
  type = "proportion"
)

plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_ots_distribution", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_topology_deviation")
### * lot_plot_topology_deviation

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_topology_deviation
### Title: Plot topology deviations within clusters
### Aliases: lot_plot_topology_deviation

### ** Examples

data("lot_glucold_ots")

example_object <- lot_glucold_ots

example_results <- lot_results(
  example_object
)

example_results$flexmix_ots_distribution <-
  S4Vectors::DataFrame(
    cluster = c(
      "1",
      "1",
      "2",
      "2"
    ),
    group = "treated",
    topology_label = c(
      "up_down",
      "flat_flat",
      "up_down",
      "flat_flat"
    ),
    n_genes = c(
      8L,
      2L,
      3L,
      7L
    ),
    proportion = c(
      0.8,
      0.2,
      0.3,
      0.7
    )
  )

methods::slot(
  example_object,
  "results"
) <- example_results

plot <- lot_plot_topology_deviation(
  object = example_object,
  group = "treated"
)

plot



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_topology_deviation", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_plot_trajectory")
### * lot_plot_trajectory

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_plot_trajectory
### Title: Plot longitudinal trajectories for selected molecular features
### Aliases: lot_plot_trajectory

### ** Examples

data("lot_glucold_ots")

lot_plot_trajectory(
  lot_glucold_ots,
  features = c("ATP2A3", "CAT"),
  groups = c(
    "ICS 30 months",
    "ICS 6 months then withdrawal",
    "Placebo"
  ),
  estimator = "masigpro",
  show_subjects = TRUE
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_plot_trajectory", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_query")
### * lot_query

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_query
### Title: Query LongOmicsTraj results
### Aliases: lot_query NOT lot_query,LongOmicsTraj-method

### ** Examples

data("lot_glucold_ots")

flat <- lot_query(
  object = lot_glucold_ots,
  topology = "flat_flat"
)

nrow(flat)

utils::head(
  as.data.frame(flat)
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_query", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_summary")
### * lot_summary

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_summary
### Title: Summarise a LongOmicsTraj object
### Aliases: lot_summary

### ** Examples

data("lot_glucold_ots")

result <- lot_summary(
  lot_glucold_ots,
  estimator = "masigpro"
)

result$assay_summary
result$top_ots



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_summary", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("lot_trajectory_estimators")
### * lot_trajectory_estimators

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: lot_from_empirical
### Title: Estimate trajectories using a maSigPro-style polynomial model
### Aliases: lot_from_empirical lot_from_lmm lot_from_gam lot_from_masigpro
###   lot_from_empirical,LongOmicsTraj-method
###   lot_from_gam,LongOmicsTraj-method lot_from_lmm,LongOmicsTraj-method
###   lot_from_masigpro,LongOmicsTraj-method lot_trajectory_estimators

### ** Examples

data("lot_glucold")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

fitted <- lot_from_empirical(
  object = lot_glucold,
  assay = "transcriptomics",
  visits = visits,
  overwrite = TRUE
)

utils::head(
  as.data.frame(
    lot_vmm(fitted)
  )
)



base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("lot_trajectory_estimators", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
cleanEx()
nameEx("show-LongOmicsTraj-method")
### * show-LongOmicsTraj-method

flush(stderr()); flush(stdout())

base::assign(".ptime", proc.time(), pos = "CheckExEnv")
### Name: show,LongOmicsTraj-method
### Title: Display a LongOmicsTraj object
### Aliases: show,LongOmicsTraj-method

### ** Examples

data("lot_glucold")

## Printing the object invokes the show() method.
lot_glucold

## The object is returned invisibly.
invisible(lot_glucold)




base::assign(".dptime", (proc.time() - get(".ptime", pos = "CheckExEnv")), pos = "CheckExEnv")
base::cat("show-LongOmicsTraj-method", base::get(".format_ptime", pos = 'CheckExEnv')(get(".dptime", pos = "CheckExEnv")), "\n", file=base::get(".ExTimings", pos = 'CheckExEnv'), append=TRUE, sep="\t")
### * <FOOTER>
###
cleanEx()
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
