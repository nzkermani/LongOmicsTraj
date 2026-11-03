# LongOmicsTraj NEWS

## Development (spatial-visit-centring)

### New features

* `lot_compute_deltas()` gains a `center` argument (`"none"`, the default, `"median"` or `"mean"`). Centring removes visit-wide offsets that shift every feature in the same direction before changes are computed, which is needed when ordered "visits" are spatial windows (for example tissue depth) that can differ globally in capture efficiency or composition. The removed offsets are stored in `results$visit_centring`. Default behaviour is unchanged.
* `center = "relative"` removes only each group's deviation from the change shared by all groups, keeping shared (likely biological) shifts and removing group-specific (likely technical) ones. Recommended over `"median"` when there are three or more groups and many features may shift in one direction.
* New spatial helpers: `lot_consensus_restart()` picks a depth-model restart by agreement across restarts (label-free); `lot_spatial_pseudobulk()` turns spot- or cell-level counts plus a depth coordinate or layer labels into depth-window pseudobulk ready for `lot_create()`; `lot_topology_agreement()` measures how often features share the same topology across groups, with a permutation test.
* New vignette, "Spatial omics with LongOmicsTraj: a practical guide", taking spot-level counts through depth selection, pseudobulk, reproducibility, centring and cross-group queries.

## LongOmicsTraj 0.99.0

### New features

* Initial Bioconductor submission.
* Added methods for estimating longitudinal molecular trajectories using empirical, linear mixed-effects, generalized additive model (GAM), and maSigPro estimators.
* Introduced Omics Trajectory Signatures (OTS) for classifying temporal molecular behaviour.
* Implemented topology-aware clustering and cluster evaluation utilities.
* Added a flexible topology-based query framework for retrieving molecular features matching simple or complex trajectory patterns across experimental groups and treatment arms.
* Included publication-quality visualisation functions for trajectories, topology composition, and cluster summaries.
* Added comprehensive documentation, runnable examples, unit tests, and vignettes.