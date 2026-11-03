# WSSmod 0.1.0

First release.

## Module scores

* `calculate_WSS()` computes Weighted Sum Score module-level summaries
  from proteomic data, using STRING network connectivity and
  biomarker-association weighting to align and aggregate protein-level
  signal into module scores. Supports `dedup_by` and `min_module_size`
  options for controlling how genes assigned to multiple modules are
  deduplicated and how small modules are filtered.
* `load_prebuilt_wss()`, `list_prebuilt_wss()`, and `load_prebuilt_terms()`
  provide access to a bundled reference module/weighting table
  (`"core_AD_plasma_biomarkers"`) and its module term labels.
* `calculate_WSS()` accepts a data frame as `result` (converted to a
  data.table) and rejects other classes with an informative error.
* `calculate_WSS()` stops with a clear error that names the affected
  proteins when a protein used by a module has a missing value. This also
  covers the case where every module has one, which previously returned
  `NA` scores without any message.
* The documentation now states that a module with none of its proteins in the
  matrix is reported with all-`NA` scores and a `module_sd` of 0.

## Prediction from core biomarkers

* `predict_WSS()` predicts all 75 module scores for the
  `"core_AD_plasma_biomarkers"` reference set from Age, Gender, and 8 core
  plasma biomarkers, using a bundled `joinet` two-layer elastic net model,
  without requiring the full proteomic panel.
* `normalize_wss_biomarkers()` builds the rank-normalized biomarker inputs
  `predict_WSS()` requires from raw biomarker values, via either
  `method = "project"` (projecting onto a bundled, anonymized baseline
  reference distribution, using `project_rank_norm()` -- supports a single
  new patient) or `method = "self"` (rank-normalizing within your own
  cohort).
* `project_rank_norm()` and `load_prebuilt_biomarker_reference()` expose
  the underlying projection utility and reference distribution directly.

## Documentation

* Added a citation (`citation("WSSmod")`) and links to the accompanying preprint.
* New vignette, "Calculating Weighted Sum Scores", covering input normalization,
  the bundled and custom module tables, and scoring new samples against a
  reference cohort.
* New vignette, "Predicting Weighted Sum Scores from plasma biomarkers", covering
  the biomarker inputs, the two normalization methods, and reading the
  predicted module scores.
* New vignette, "Exploring the bundled datasets", describing each bundled file, its
  loader and what it can be used for.
* A package-level help page (`?WSSmod`) gives an overview of the functions
  and bundled data.
* `calculate_WSS()` now documents what input scale it expects.

## Example data

* `wss_example_expression` and `wss_example_biomarkers` are small, fully
  simulated datasets (no participant data) for trying out `calculate_WSS()`,
  `normalize_wss_biomarkers()` and `predict_WSS()`.

The bundled model and reference data were verified end-to-end against the
original analysis's held-out follow-up cohort, reproducing published
predictions to floating-point precision.
