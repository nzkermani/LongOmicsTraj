## ----setup, message=FALSE, warning=FALSE--------------------------------------
library(LongOmicsTraj)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(S4Vectors)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(scales)

knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.align = "center",
  message = FALSE,
  warning = FALSE
)

## ----load-data----------------------------------------------------------------
# Load the example GLUCOLD expression matrix included with the package.
# The demo contains 71 selected genes measured across 221 longitudinal samples.
# The reduced feature set keeps the example fast to run.
data(
  "lot_glucold_expression",
  package = "LongOmicsTraj"
)

# Load the matching clinical metadata.
# There is one metadata row for each expression-matrix sample.
#
# Required columns:
#   sample_id  = unique sample identifier; must match colnames(expr_use)
#   subject_id = participant identifier; repeated across longitudinal visits
#   visit      = visit label, for example Baseline, 6 Months, 30 Months
#   group      = treatment or clinical group
#
# A numeric time column may also be added later for models such as GAM:
#   time       = numerical visit time, for example 0, 6, 30
data(
  "lot_glucold_clinical",
  package = "LongOmicsTraj"
)

# Replace these two objects with your own data.
#
# Expression matrix:
#   rows    = genes or molecular features
#   columns = samples
expr_use <- lot_glucold_expression

# Clinical metadata:
#   rows = samples
#   required columns = sample_id, subject_id, visit, group
#   sample_id values must match the expression-matrix column names
meta_use <- lot_glucold_clinical

## ----inspect-expression-------------------------------------------------------
dim(expr_use)

expr_use[1:5, 1:5]

## ----inspect-metadata---------------------------------------------------------
dim(meta_use)

head(meta_use)

## ----metadata-columns---------------------------------------------------------
colnames(meta_use)

## ----user-input-template, eval=FALSE------------------------------------------
# expr_use <- your_expression_matrix
# meta_use <- your_clinical_metadata

## ----align-samples------------------------------------------------------------
expr_use <- as.matrix(expr_use)
storage.mode(expr_use) <- "numeric"

meta_use <- as.data.frame(meta_use)

required_columns <- c(
  "sample_id",
  "subject_id",
  "visit",
  "group"
)

missing_columns <- setdiff(
  required_columns,
  colnames(meta_use)
)

if (length(missing_columns) > 0L) {
  stop(
    "Clinical metadata are missing: ",
    paste(missing_columns, collapse = ", ")
  )
}

if (is.null(colnames(expr_use))) {
  stop(
    "The expression matrix must have sample IDs as column names."
  )
}

if (!setequal(
  colnames(expr_use),
  meta_use$sample_id
)) {
  stop(
    "Expression sample names do not match metadata sample IDs."
  )
}

meta_use <- meta_use[
  match(
    colnames(expr_use),
    meta_use$sample_id
  ),
  ,
  drop = FALSE
]

rownames(meta_use) <- meta_use$sample_id

stopifnot(
  identical(
    colnames(expr_use),
    meta_use$sample_id
  )
)

## ----define-visits------------------------------------------------------------
clinical_visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

meta_use$visit <- factor(
  meta_use$visit,
  levels = clinical_visits,
  ordered = TRUE
)

time_lookup <- c(
  "Baseline" = 0,
  "6 Months" = 6,
  "30 Months" = 30
)

meta_use$time <- unname(
  time_lookup[
    as.character(meta_use$visit)
  ]
)

if (anyNA(meta_use$time)) {
  stop(
    "One or more visit labels could not be converted to time."
  )
}

table(meta_use$visit)

table(meta_use$group)

## ----create-se----------------------------------------------------------------
se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    transcriptomics = expr_use
  ),
  colData = S4Vectors::DataFrame(
    meta_use,
    row.names = meta_use$sample_id
  )
)

se

## ----create-mae---------------------------------------------------------------
mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments =
    MultiAssayExperiment::ExperimentList(
      transcriptomics = se
    )
)

mae

## ----create-lot---------------------------------------------------------------
lot_object <- methods::new(
  "LongOmicsTraj",
  mae = mae
)

methods::validObject(lot_object)

lot_object

## ----empirical-fit------------------------------------------------------------
lot_emp <- lot_from_empirical(
  object = lot_object,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_emp <- lot_compute_deltas(
  object = lot_emp,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_emp <- lot_compute_thresholds(
  object = lot_emp,
  assay = "transcriptomics",
  method = "hybrid",
  delta_min = 0.15,
  k = 0.5,
  visits = clinical_visits,
  overwrite = TRUE
)

lot_emp <- lot_map_to_ots(
  object = lot_emp,
  visits = clinical_visits,
  overwrite = TRUE
)

head(
  as.data.frame(lot_emp@ots)
)

## ----masigpro-fit-------------------------------------------------------------
lot_mas <- lot_from_masigpro(
  object = lot_object,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_mas <- lot_compute_deltas(
  object = lot_mas,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_mas <- lot_compute_thresholds(
  object = lot_mas,
  assay = "transcriptomics",
  method = "hybrid",
  delta_min = 0.15,
  k = 0.5,
  visits = clinical_visits,
  overwrite = TRUE
)

lot_mas <- lot_map_to_ots(
  object = lot_mas,
  visits = clinical_visits,
  overwrite = TRUE
)

head(
  as.data.frame(lot_mas@ots)
)

## ----ots-distribution-data----------------------------------------------------
ots_results <- as.data.frame(
  lot_mas@ots
) |>
  dplyr::filter(
    estimator == "masigpro",
    !is.na(group),
    !is.na(topology_label)
  )

ots_distribution <- ots_results |>
  dplyr::count(
    group,
    topology_label,
    name = "n"
  ) |>
  dplyr::group_by(group) |>
  dplyr::mutate(
    proportion = n / sum(n)
  ) |>
  dplyr::ungroup()

## ----ots-distribution-plot, fig.width=11, fig.height=6------------------------
p_ots_distribution <- ggplot2::ggplot(
  ots_distribution,
  ggplot2::aes(
    x = group,
    y = proportion,
    fill = topology_label
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(),
    expand = ggplot2::expansion(
      mult = c(0, 0.02)
    )
  ) +
  ggplot2::scale_fill_brewer(
    palette = "Set3",
    name = "OTS"
  ) +
  ggplot2::labs(
    title = "Distribution of Omics Trajectory Signatures",
    subtitle = "maSigPro estimates across GLUCOLD treatment groups",
    x = "Treatment group",
    y = "Proportion of genes"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

p_ots_distribution

## ----treatment-query----------------------------------------------------------
withdrawal_specific <- lot_query(
  lot_mas,
  estimator = "masigpro",
  conditions = list(
    "ICS 6 months then withdrawal" = c(
      "up_down",
      "down_up"
    ),
    "ICS 30 months" = c(
      "up_up",
      "up_flat",
      "down_down",
      "down_flat"
    ),
    "Placebo" = NOT(
      c(
        "up_down",
        "down_up"
      )
    )
  )
)

withdrawal_specific <- as.data.frame(
  withdrawal_specific
)

head(withdrawal_specific)

## ----query-summary-data-------------------------------------------------------
query_distribution <- withdrawal_specific |>
  dplyr::count(
    group,
    topology_label,
    name = "n"
  )

## ----query-summary-plot, fig.width=10, fig.height=6---------------------------
p_query <- ggplot2::ggplot(
  query_distribution,
  ggplot2::aes(
    x = topology_label,
    y = n,
    fill = group
  )
) +
  ggplot2::geom_col(
    position = "dodge",
    width = 0.72
  ) +
  ggplot2::labs(
    title = "Treatment-specific trajectory query",
    subtitle = "Genes satisfying the predefined GLUCOLD trajectory conditions",
    x = "Omics Trajectory Signature",
    y = "Number of genes",
    fill = "Treatment group"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 35,
      hjust = 1
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

p_query

## ----select-genes-------------------------------------------------------------
example_genes <- withdrawal_specific |>
  dplyr::distinct(
    object_id,
    topology_label
  ) |>
  dplyr::group_by(topology_label) |>
  dplyr::slice_head(n = 1) |>
  dplyr::ungroup() |>
  dplyr::slice_head(n = 6) |>
  dplyr::pull(object_id) |>
  unique()

example_genes <- intersect(
  example_genes,
  rownames(expr_use)
)

if (length(example_genes) == 0L) {
  stop(
    "No query-selected genes were found in the expression matrix."
  )
}

example_genes

## ----prepare-gene-data--------------------------------------------------------
expr <- expr_use
clinical_data <- meta_use

expression_long <- as.data.frame(
  expr[
    example_genes,
    ,
    drop = FALSE
  ]
) |>
  tibble::rownames_to_column(
    "gene"
  ) |>
  tidyr::pivot_longer(
    cols = -gene,
    names_to = "sample_id",
    values_to = "expression"
  ) |>
  dplyr::left_join(
    clinical_data |>
      dplyr::select(
        sample_id,
        subject_id,
        visit,
        group
      ),
    by = "sample_id"
  )

selected_groups <- c(
  "ICS 30 months",
  "ICS 6 months then withdrawal",
  "Placebo"
)

expression_long <- expression_long |>
  dplyr::filter(
    group %in% selected_groups
  ) |>
  dplyr::mutate(
    group = factor(
      group,
      levels = selected_groups
    ),
    visit = factor(
      visit,
      levels = clinical_visits,
      ordered = TRUE
    )
  )

expression_summary <- expression_long |>
  dplyr::group_by(
    gene,
    group,
    visit
  ) |>
  dplyr::summarise(
    n = sum(!is.na(expression)),
    mean_expression = mean(
      expression,
      na.rm = TRUE
    ),
    sd_expression = stats::sd(
      expression,
      na.rm = TRUE
    ),
    se_expression = dplyr::if_else(
      n > 1,
      sd_expression / sqrt(n),
      0
    ),
    .groups = "drop"
  )

## ----prepare-labels-----------------------------------------------------------
gene_ots_labels <- ots_results |>
  dplyr::filter(
    object_id %in% example_genes,
    group == "ICS 6 months then withdrawal"
  ) |>
  dplyr::distinct(
    object_id,
    topology_label
  ) |>
  dplyr::mutate(
    facet_label = paste0(
      object_id,
      "\nWithdrawal OTS: ",
      topology_label
    )
  )

facet_labels <- stats::setNames(
  gene_ots_labels$facet_label,
  gene_ots_labels$object_id
)

missing_labels <- setdiff(
  example_genes,
  names(facet_labels)
)

facet_labels[missing_labels] <- missing_labels

## ----representative-trajectories, fig.width=11, fig.height=8------------------
p_gene_trajectories <- ggplot2::ggplot() +

  ggplot2::geom_line(
    data = expression_long,
    ggplot2::aes(
      x = visit,
      y = expression,
      group = interaction(
        group,
        subject_id
      ),
      colour = group
    ),
    alpha = 0.10,
    linewidth = 0.25,
    show.legend = FALSE
  ) +

  ggplot2::geom_errorbar(
    data = expression_summary,
    ggplot2::aes(
      x = visit,
      ymin = mean_expression - se_expression,
      ymax = mean_expression + se_expression,
      colour = group
    ),
    width = 0.08,
    linewidth = 0.45
  ) +

  ggplot2::geom_line(
    data = expression_summary,
    ggplot2::aes(
      x = visit,
      y = mean_expression,
      group = group,
      colour = group
    ),
    linewidth = 1
  ) +

  ggplot2::geom_point(
    data = expression_summary,
    ggplot2::aes(
      x = visit,
      y = mean_expression,
      colour = group
    ),
    size = 2.2
  ) +

  ggplot2::facet_wrap(
    ~ gene,
    scales = "free_y",
    ncol = 3,
    labeller = ggplot2::as_labeller(
      facet_labels
    )
  ) +

  ggplot2::labs(
    title = "Representative longitudinal gene trajectories",
    subtitle = "Genes selected from the treatment-specific maSigPro query",
    x = "Visit",
    y = "Expression",
    colour = "Treatment group",
    caption = paste0(
      "Faint lines show individual participants. ",
      "Solid lines show group means; error bars show ±SE. ",
      "Panels use independent y-axis ranges."
    )
  ) +

  ggplot2::theme_bw(base_size = 11) +

  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold"
    ),
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 30,
      hjust = 1
    ),
    plot.caption = ggplot2::element_text(
      hjust = 0
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

p_gene_trajectories

## ----own-data-template, eval=FALSE--------------------------------------------
# expr_use <- your_expression_matrix
# meta_use <- your_clinical_metadata
# 
# meta_use <- meta_use[
#   match(
#     colnames(expr_use),
#     meta_use$sample_id
#   ),
#   ,
#   drop = FALSE
# ]
# 
# rownames(meta_use) <- meta_use$sample_id
# 
# se <- SummarizedExperiment(
#   assays = list(
#     transcriptomics = expr_use
#   ),
#   colData = DataFrame(meta_use)
# )
# 
# mae <- MultiAssayExperiment(
#   experiments = ExperimentList(
#     transcriptomics = se
#   )
# )
# 
# lot_object <- new(
#   "LongOmicsTraj",
#   mae = mae
# )

## ----alternative-assays, eval=FALSE-------------------------------------------
# assays = list(
#   proteomics = protein_matrix
# )

## ----multiple-assays, eval=FALSE----------------------------------------------
# mae <- MultiAssayExperiment(
#   experiments = ExperimentList(
#     transcriptomics = transcriptomics_se,
#     proteomics = proteomics_se,
#     metabolomics = metabolomics_se
#   )
# )

## ----quick-reference, eval=FALSE----------------------------------------------
# # Load data
# data("lot_glucold_expression")
# data("lot_glucold_clinical")
# 
# # Create object
# se <- SummarizedExperiment(
#   assays = list(
#     transcriptomics = lot_glucold_expression
#   ),
#   colData = DataFrame(
#     lot_glucold_clinical
#   )
# )
# 
# mae <- MultiAssayExperiment(
#   experiments = ExperimentList(
#     transcriptomics = se
#   )
# )
# 
# lot <- new(
#   "LongOmicsTraj",
#   mae = mae
# )
# 
# # Fit trajectories
# lot <- lot_from_masigpro(
#   lot,
#   assay = "transcriptomics",
#   visits = c(
#     "Baseline",
#     "6 Months",
#     "30 Months"
#   )
# )
# 
# # Compute OTS
# lot <- lot_compute_deltas(
#   lot,
#   assay = "transcriptomics",
#   visits = c(
#     "Baseline",
#     "6 Months",
#     "30 Months"
#   )
# )
# 
# lot <- lot_compute_thresholds(
#   lot,
#   assay = "transcriptomics",
#   method = "hybrid",
#   delta_min = 0.15,
#   k = 0.5,
#   visits = c(
#     "Baseline",
#     "6 Months",
#     "30 Months"
#   )
# )
# 
# lot <- lot_map_to_ots(
#   lot,
#   visits = c(
#     "Baseline",
#     "6 Months",
#     "30 Months"
#   )
# )
# 
# # Inspect results
# head(
#   as.data.frame(lot@ots)
# )

