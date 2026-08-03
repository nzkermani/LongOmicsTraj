###############################################################
## LongOmicsTraj: GLUCOLD example analysis
###############################################################

library(LongOmicsTraj)
library(SummarizedExperiment)
library(MultiAssayExperiment)
library(S4Vectors)
library(mgcv)

data("lot_glucold_expression")
data("lot_glucold_clinical")

expr_use <- lot_glucold_expression
meta_use <- lot_glucold_clinical

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
    paste(missing_columns, collapse = ", "),
    call. = FALSE
  )
}

if (is.null(colnames(expr_use))) {
  stop(
    "The expression matrix must have sample IDs as column names.",
    call. = FALSE
  )
}

if (!setequal(colnames(expr_use), meta_use$sample_id)) {
  stop(
    "The expression-matrix column names do not match ",
    "the clinical-data sample IDs.",
    call. = FALSE
  )
}

meta_use <- meta_use[
  match(colnames(expr_use), meta_use$sample_id),
  ,
  drop = FALSE
]

rownames(meta_use) <- meta_use$sample_id

stopifnot(
  identical(
    colnames(expr_use),
    rownames(meta_use)
  )
)

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
  time_lookup[as.character(meta_use$visit)]
)

if (anyNA(meta_use$time)) {
  stop(
    "Some visit labels could not be converted to numerical time.",
    call. = FALSE
  )
}

se <- SummarizedExperiment::SummarizedExperiment(
  assays = list(
    expression = expr_use
  ),
  colData = S4Vectors::DataFrame(
    meta_use,
    row.names = meta_use$sample_id
  )
)

mae <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = MultiAssayExperiment::ExperimentList(
    transcriptomics = se
  )
)

lot_object <- methods::new(
  "LongOmicsTraj",
  mae = mae
)

methods::validObject(lot_object)

message("Processing empirical estimator...")

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

message("Processing linear mixed model estimator...")

lot_lmm <- lot_from_lmm(
  object = lot_object,
  assay = "transcriptomics",
  formula_fixed = y ~ visit + group + (1 | subject_id),
  visits = clinical_visits,
  overwrite = TRUE
)

lot_lmm <- lot_compute_deltas(
  object = lot_lmm,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_lmm <- lot_compute_thresholds(
  object = lot_lmm,
  assay = "transcriptomics",
  method = "hybrid",
  delta_min = 0.15,
  k = 0.5,
  visits = clinical_visits,
  overwrite = TRUE
)

lot_lmm <- lot_map_to_ots(
  object = lot_lmm,
  visits = clinical_visits,
  overwrite = TRUE
)

message("Processing GAM estimator...")

clinical_formula_gam <-
  y ~ group +
  mgcv::s(time, by = group, k = 3, bs = "cr") +
  mgcv::s(subject_id, bs = "re")

lot_gam <- lot_from_gam(
  object = lot_object,
  assay = "transcriptomics",
  formula_fixed = clinical_formula_gam,
  visits = clinical_visits,
  group_col = "group",
  subject_col = "subject_id",
  overwrite = TRUE
)

lot_gam <- lot_compute_deltas(
  object = lot_gam,
  assay = "transcriptomics",
  visits = clinical_visits,
  overwrite = TRUE
)

lot_gam <- lot_compute_thresholds(
  object = lot_gam,
  assay = "transcriptomics",
  method = "hybrid",
  delta_min = 0.15,
  k = 0.5,
  visits = clinical_visits,
  overwrite = TRUE
)

lot_gam <- lot_map_to_ots(
  object = lot_gam,
  visits = clinical_visits,
  overwrite = TRUE
)

message("Processing maSigPro estimator...")

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

message("Creating consensus trajectories...")

compared_df <- lot_compare_estimators(
  lot_gam,
  lot_mas,
  lot_lmm,
  lot_emp
)

consensus_output <- lot_consensus(
  x = compared_df,
  min_agreement = 2L,
  mode = "majority"
)

extract_ots <- function(object, method_name) {
  ots_data <- as.data.frame(object@ots)

  required_ots_columns <- c(
    "object_id",
    "group",
    "topology_label"
  )

  missing_ots_columns <- setdiff(
    required_ots_columns,
    colnames(ots_data)
  )

  if (length(missing_ots_columns) > 0L) {
    stop(
      "OTS results are missing: ",
      paste(missing_ots_columns, collapse = ", "),
      call. = FALSE
    )
  }

  ots_data |>
    dplyr::select(
      dplyr::all_of(required_ots_columns)
    ) |>
    dplyr::mutate(
      Method = method_name
    )
}

df_emp <- extract_ots(lot_emp, "Empirical")
df_lmm <- extract_ots(lot_lmm, "LMM")
df_gam <- extract_ots(lot_gam, "GAM")
df_mas <- extract_ots(lot_mas, "maSigPro")

consensus_data <- as.data.frame(consensus_output)

if ("group" %in% colnames(consensus_data)) {
  df_cons <- consensus_data |>
    dplyr::transmute(
      object_id = .data$object_id,
      group = .data$group,
      topology_label = .data$consensus_topology,
      Method = "Consensus"
    )
} else {
  df_cons <- consensus_data |>
    dplyr::transmute(
      object_id = .data$object_id,
      group = NA_character_,
      topology_label = .data$consensus_topology,
      Method = "Consensus"
    )
}

clinical_ots_compiled <- dplyr::bind_rows(
  df_emp,
  df_lmm,
  df_gam,
  df_mas,
  df_cons
) |>
  dplyr::rename(
    Feature = .data$object_id,
    Treatment_Group = .data$group,
    Predicted_OTS = .data$topology_label
  )

output_dir <- file.path(
  getwd(),
  "glucold_results"
)

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

utils::write.csv(
  clinical_ots_compiled,
  file = file.path(
    output_dir,
    "glucold_predicted_ots_distribution.csv"
  ),
  row.names = FALSE
)

save(
  lot_object,
  lot_emp,
  lot_lmm,
  lot_gam,
  lot_mas,
  compared_df,
  consensus_output,
  clinical_ots_compiled,
  file = file.path(
    output_dir,
    "glucold_longomicstraj_analysis.RData"
  )
)

plot_data <- clinical_ots_compiled |>
  dplyr::filter(
    !is.na(.data$Treatment_Group),
    !is.na(.data$Predicted_OTS)
  ) |>
  dplyr::count(
    .data$Method,
    .data$Treatment_Group,
    .data$Predicted_OTS,
    name = "Count"
  ) |>
  dplyr::group_by(
    .data$Method,
    .data$Treatment_Group
  ) |>
  dplyr::mutate(
    Proportion = .data$Count / sum(.data$Count)
  ) |>
  dplyr::ungroup()

p_distribution <- ggplot2::ggplot(
  plot_data,
  ggplot2::aes(
    x = .data$Treatment_Group,
    y = .data$Proportion,
    fill = .data$Predicted_OTS
  )
) +
  ggplot2::geom_col(
    width = 0.72,
    colour = "white",
    linewidth = 0.25
  ) +
  ggplot2::facet_wrap(
    facets = ggplot2::vars(.data$Method),
    ncol = 5
  ) +
  ggplot2::scale_y_continuous(
    labels = scales::percent_format(),
    expand = ggplot2::expansion(
      mult = c(0, 0.02)
    )
  ) +
  ggplot2::scale_fill_brewer(
    palette = "Set3",
    name = "Assigned OTS"
  ) +
  ggplot2::labs(
    title = "Distribution of Omics Trajectory Signatures in GLUCOLD",
    subtitle = paste0(
      "Comparison of trajectory estimators across three visits; ",
      length(unique(clinical_ots_compiled$Feature)),
      " features analysed"
    ),
    x = "Treatment group",
    y = "Proportion of classified features"
  ) +
  ggplot2::theme_bw(base_size = 11) +
  ggplot2::theme(
    plot.title = ggplot2::element_text(
      face = "bold",
      hjust = 0.5
    ),
    plot.subtitle = ggplot2::element_text(
      hjust = 0.5,
      margin = ggplot2::margin(b = 12)
    ),
    strip.background = ggplot2::element_rect(
      fill = "grey95"
    ),
    strip.text = ggplot2::element_text(
      face = "bold"
    ),
    axis.text.x = ggplot2::element_text(
      angle = 45,
      hjust = 1
    ),
    panel.grid.minor = ggplot2::element_blank(),
    legend.position = "bottom"
  )

print(p_distribution)

ggplot2::ggsave(
  filename = file.path(
    output_dir,
    "figure_glucold_ots_distribution.png"
  ),
  plot = p_distribution,
  width = 13,
  height = 6.5,
  units = "in",
  dpi = 300,
  bg = "white"
)
