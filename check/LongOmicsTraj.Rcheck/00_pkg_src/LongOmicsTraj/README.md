# LongOmicsTraj

**Topology-based analysis of longitudinal omics trajectories**

<p align="center">
  <img src="man/figures/LongOmicsTraj_overview.png"
       alt="LongOmicsTraj workflow overview"
       width="1000">
</p>

## Why LongOmicsTraj?

Longitudinal omics studies measure molecular features across multiple visits,
treatment periods or disease stages.

Traditional methods can identify statistical changes or smooth temporal
patterns, but their outputs are often difficult to translate into concise,
biologically interpretable trajectory classes.

LongOmicsTraj converts consecutive molecular changes into **Omics Trajectory
Signatures (OTS)**, providing an interpretable description of longitudinal
behaviour.

## Core idea

Many commonly used approaches answer related but different questions:

- Mixed-effects models ask whether expression changes over time.
- Generalised additive models estimate smooth temporal patterns.
- Clustering methods identify features with similar profiles.
- LongOmicsTraj identifies features that share the same trajectory topology.

LongOmicsTraj therefore asks:

> **Which molecular features share the same longitudinal trajectory topology?**

## Omics Trajectory Signatures

For three ordered visits, each feature receives a two-transition topology
label.

Examples include:

| OTS label | Interpretation |
|---|---|
| `down_up` | Rebound |
| `up_down` | Transient response |
| `flat_flat` | Stable trajectory |
| `down_down` | Sustained suppression |
| `up_up` | Sustained activation |
| `flat_up` | Delayed activation |
| `flat_down` | Delayed suppression |

The number of transition states increases naturally when more visits are
included.

## Basic workflow

```r
library(LongOmicsTraj)

data("lot_glucold")

visits <- c(
  "Baseline",
  "6 Months",
  "30 Months"
)

lot <- lot_from_empirical(
  object = lot_glucold,
  assay = "transcriptomics",
  visits = visits,
  overwrite = TRUE
)

lot <- lot_compute_deltas(
  object = lot,
  assay = "transcriptomics",
  visits = visits,
  overwrite = TRUE
)

lot <- lot_compute_thresholds(
  object = lot,
  assay = "transcriptomics",
  method = "hybrid",
  visits = visits,
  overwrite = TRUE
)

lot <- lot_map_to_ots(
  object = lot,
  visits = visits,
  overwrite = TRUE
)
```

## Inspect the results

```r
lot_summary(
  lot
)
```

```r
head(
  as.data.frame(
    lot_ots(lot)
  )
)
```

## Query trajectory classes

```r
stable_features <- lot_query(
  object = lot,
  topology = "flat_flat"
)
```

Cross-group conditions can also be specified:

```r
selected_features <- lot_query(
  object = lot,
  conditions = list(
    Placebo = "flat_flat",
    `ICS 30 months` = NOT("flat_flat")
  )
)
```

## Visualise trajectories

```r
lot_plot_ots_distribution(
  object = lot,
  type = "proportion"
)
```

```r
lot_plot_trajectory(
  object = lot,
  features = c(
    "ATP2A3",
    "CAT"
  ),
  assay = "transcriptomics"
)
```

## Trajectory estimators

LongOmicsTraj supports several approaches for estimating visit-level molecular
trajectories:

- Empirical visit means
- Linear mixed-effects models
- Generalised additive models
- maSigPro-style polynomial regression
- Molecular-feature trajectory clustering
- FlexMix trajectory clustering

Estimator outputs are converted into a common visit mean matrix before
computing deltas, thresholds and OTS assignments.

## Cluster-level topology analysis

LongOmicsTraj can cluster molecular features according to their longitudinal
profiles and quantify the topology composition within each cluster.

```r
lot_clusters <- lot_from_clusters(
  object = lot,
  assay = "transcriptomics",
  visits = visits,
  k = 4,
  estimator = "kmeans",
  overwrite = TRUE
)
```

Available cluster visualisations include:

```r
lot_plot_cluster_trajectory_arms(
  object = lot_clusters
)
```

```r
lot_plot_topology_deviation(
  object = lot_clusters
)
```

FlexMix model-selection and topology plots are also available through:

```r
lot_plot_cluster_selection()
lot_plot_cluster_topology()
```

## Comparison with related approaches

| Method | Primary output | Limitation addressed by LongOmicsTraj |
|---|---|---|
| Mixed-effects models | Statistical longitudinal effects | Trajectory topology may remain implicit |
| GAMs | Smooth temporal curves | Shared discrete trajectory structures are not directly assigned |
| K-means or PAM | Similarity-based clusters | Cluster membership does not itself define topology |
| LongOmicsTraj | Interpretable trajectory topology | Converts longitudinal changes into explicit OTS labels |

## Installation

Install the development version from GitHub:

```r
remotes::install_github(
  "nzkermani/LongOmicsTraj"
)
```

## Example data

The package includes processed GLUCOLD longitudinal transcriptomic data:

```r
data("lot_glucold")
data("lot_glucold_expression")
data("lot_glucold_clinical")
data("lot_glucold_ots")
```

The example data are derived from NCBI Gene Expression Omnibus accession
GSE36221.

## Package status

- Tested using `testthat`
- Includes reproducible example data
- Supports multiple trajectory estimators
- Supports gene-level and cluster-level topology analysis
- Under active development

## Author

**Nazanin Zounemat-Kermani**

Imperial College London