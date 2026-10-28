# LongOmicsTraj NEWS

## Development (spatial-visit-centring)

### New features

* `lot_compute_deltas()` gains a `center` argument (`"none"`, the default, `"median"` or `"mean"`). Centring removes visit-wide offsets that shift every feature in the same direction before changes are computed, which is needed when ordered "visits" are spatial windows (for example tissue depth) that can differ globally in capture efficiency or composition. The removed offsets are stored in `results$visit_centring`. Default behaviour is unchanged.

## LongOmicsTraj 0.99.0

### New features

* Initial Bioconductor submission.
* Added methods for estimating longitudinal molecular trajectories using empirical, linear mixed-effects, generalized additive model (GAM), and maSigPro estimators.
* Introduced Omics Trajectory Signatures (OTS) for classifying temporal molecular behaviour.
* Implemented topology-aware clustering and cluster evaluation utilities.
* Added a flexible topology-based query framework for retrieving molecular features matching simple or complex trajectory patterns across experimental groups and treatment arms.
* Included publication-quality visualisation functions for trajectories, topology composition, and cluster summaries.
* Added comprehensive documentation, runnable examples, unit tests, and vignettes.