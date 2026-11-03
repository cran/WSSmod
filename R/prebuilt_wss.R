# Registry of prebuilt WSS reference sets shipped in inst/extdata.
# Each entry names the scores file (module/protein weighting table, required),
# an optional terms file (human-readable module labels), an optional model
# file (a fitted joinet prediction model, see predict_WSS()), and an optional
# biomarker_reference file (see load_prebuilt_biomarker_reference() and
# normalize_wss_biomarkers()).
.wss_prebuilt_registry <- list(
  core_AD_plasma_biomarkers = list(
    scores = "core_AD_plasma_biomarkers.STRING_merged_scaled_scores.csv",
    terms = "core_AD_plasma_biomarkers.STRING_INT_WSS_terms_Zsum.csv",
    model = "core_AD_plasma_biomarkers.prediction_model.rds",
    biomarker_reference = "core_AD_plasma_biomarkers.baseline_biomarker_reference.csv"
  )
)

.wss_prebuilt_entry <- function(prebuilt) {
  entry <- .wss_prebuilt_registry[[prebuilt]]
  if (is.null(entry)) {
    stop(sprintf(
      "Unknown prebuilt WSS reference set '%s'. Available sets: %s",
      prebuilt,
      paste(names(.wss_prebuilt_registry), collapse = ", ")
    ), call. = FALSE)
  }
  entry
}

.wss_prebuilt_path <- function(file) {
  path <- system.file("extdata", file, package = "WSSmod")
  if (!nzchar(path)) {
    stop(sprintf("Could not find packaged reference file '%s'", file), call. = FALSE)
  }
  path
}

#' Load a Prebuilt WSS Module/Weighting Table
#'
#' Reads one of the module weighting tables shipped with WSSmod. The raw
#' table may assign the same protein to more than one module (e.g. once per
#' biomarker it was associated with); [calculate_WSS()] deduplicates and
#' filters this table before use via `dedup_by` and `min_module_size`.
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A data.table with (at least) columns `module`, `symbol`,
#'   `mean_beta`, and `mean_alpha_scaled`, suitable for use as the `result`
#'   argument of [calculate_WSS()].
#'
#' @seealso [calculate_WSS()], [list_prebuilt_wss()], [load_prebuilt_terms()]
#'
#' @examples
#' ref <- load_prebuilt_wss()
#' head(ref[, c("module", "symbol", "mean_beta", "mean_alpha_scaled")])
#'
#' @importFrom data.table fread
#' @export
load_prebuilt_wss <-function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  result <- fread(.wss_prebuilt_path(entry$scores))
  result[, module := paste(Biomarker, cluster, sep = ".")]

  result[]
}

#' List Available Prebuilt WSS Reference Sets
#'
#' @return A character vector of prebuilt reference set names usable as the
#'   `prebuilt` argument of [calculate_WSS()] and [load_prebuilt_wss()].
#'
#' @examples
#' list_prebuilt_wss()
#'
#' @export
list_prebuilt_wss <-function() {
  names(.wss_prebuilt_registry)
}

#' Get Module Term Labels for a Prebuilt WSS Reference Set
#'
#' Some prebuilt reference sets ship a companion table of human-readable
#' term labels for each module, summarized from the underlying pathway
#' enrichment used to build the module.
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A data.table of module term labels, or `NULL` if the given
#'   reference set has no associated terms file.
#'
#' @examples
#' head(load_prebuilt_terms())
#'
#' @importFrom data.table fread
#' @export
load_prebuilt_terms <-function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  if (is.null(entry$terms)) {
    return(NULL)
  }

  fread(.wss_prebuilt_path(entry$terms))
}

#' Load a Prebuilt WSS Prediction Model
#'
#' Some prebuilt reference sets ship a fitted `joinet` model (a two-layer
#' elastic net) that predicts all of that reference set's module scores from
#' a small set of covariates, without needing to measure the full proteomic
#' panel. See [predict_WSS()] to predict from this model directly.
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A list with components `model` (the fitted `joinet` object),
#'   `outcomes` (the module names it predicts), `x_cols` (the required
#'   predictor columns, in order), and `variant`; or `NULL` if the given
#'   reference set has no associated prediction model.
#'
#' @seealso [predict_WSS()], [list_prebuilt_wss()]
#'
#' @examples
#' model_info <- load_prebuilt_wss_model()
#' model_info$variant
#' model_info$x_cols
#' head(model_info$outcomes)
#'
#' @export
load_prebuilt_wss_model <-function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  if (is.null(entry$model)) {
    return(NULL)
  }

  readRDS(.wss_prebuilt_path(entry$model))
}

#' Load a Prebuilt Biomarker Reference Distribution
#'
#' Some prebuilt reference sets ship a baseline biomarker reference
#' distribution used by [normalize_wss_biomarkers()] (`method = "project"`)
#' to rank-normalize new raw biomarker values without needing a sizeable
#' cohort of your own. Each biomarker's reference values are drawn from the
#' real baseline cohort, but **independently randomly shuffled** so that no
#' two values across different biomarkers can be attributed to the same
#' original subject. This preserves each biomarker's marginal distribution
#' exactly (every value, tie, and extreme is real) while destroying the
#' joint/multivariate fingerprint that carries re-identification risk.
#'
#' Because of this, **the returned reference has no valid joint structure**
#' — inter-biomarker correlations are deliberately destroyed. Only use it
#' for independent per-biomarker operations (rank/quantile lookups), never
#' for anything requiring multivariate relationships (e.g. correlation or
#' joint modeling).
#'
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#'
#' @return A data.table with columns `biomarker` and `value` (long format,
#'   one row per reference observation), or `NULL` if the given reference set
#'   has no associated biomarker reference.
#'
#' @seealso [normalize_wss_biomarkers()], [project_rank_norm()]
#'
#' @examples
#' ref <- load_prebuilt_biomarker_reference()
#' table(ref$biomarker)
#'
#' @importFrom data.table fread
#' @export
load_prebuilt_biomarker_reference <- function(prebuilt = "core_AD_plasma_biomarkers") {
  entry <- .wss_prebuilt_entry(prebuilt)

  if (is.null(entry$biomarker_reference)) {
    return(NULL)
  }

  fread(.wss_prebuilt_path(entry$biomarker_reference))
}
