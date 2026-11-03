#' Simulated Expression Matrix for WSS Examples
#'
#' A small, fully simulated samples-by-proteins matrix for trying out
#' [calculate_WSS()]. Its columns are the 887 proteins (gene symbols) that
#' survive the default deduplication and module-size filter of the bundled
#' `"core_AD_plasma_biomarkers"` reference set, so
#' `calculate_WSS(wss_example_expression)` scores all 75 modules without
#' warnings. Values are on a z-score-like scale, as expected for normalized
#' input.
#'
#' Each module was given a latent per-sample factor, and its proteins load on
#' that factor in the direction of their `mean_beta`, plus independent noise.
#' The data are random draws and contain no participant data; they carry no
#' biological meaning and are for demonstration only.
#'
#' @format A numeric matrix with 50 rows (samples, named `sample01` to
#'   `sample50`) and 887 columns (proteins, named by gene symbol).
#'
#' @source Simulated; see `data-raw/simulate_examples.R` in the package
#'   source.
#'
#' @seealso [calculate_WSS()], [wss_example_biomarkers]
#'
#' @examples
#' dim(wss_example_expression)
#' res <- calculate_WSS(wss_example_expression)
#' head(res$scores[, 1:3])
"wss_example_expression"

#' Simulated Plasma Biomarkers for WSS Prediction Examples
#'
#' A small, fully simulated table of raw plasma biomarker values for trying
#' out [normalize_wss_biomarkers()] and [predict_WSS()]. Biomarkers are drawn
#' from log-normal distributions at plausible scales. The data are random
#' draws and contain no participant data; they carry no biological meaning
#' and are for demonstration only.
#'
#' @format A data.frame with 20 rows (named `patient01` to `patient20`) and
#'   10 columns:
#'   \describe{
#'     \item{Age}{Age in years.}
#'     \item{Gender}{`1` = male, `0` = not male.}
#'     \item{PlasmaPTau181, PlasmaAB142P, PlasmaAB140P, PlasmaABRatio,
#'       PlasmapTau217, PlasmapTau217_AB42Ratio, PlasmaGFAP, PlasmaNfL}{Raw
#'       (un-normalized) plasma biomarker values.}
#'   }
#'
#' @source Simulated; see `data-raw/simulate_examples.R` in the package
#'   source.
#'
#' @seealso [normalize_wss_biomarkers()], [predict_WSS()],
#'   [wss_example_expression]
#'
#' @examples
#' head(wss_example_biomarkers)
"wss_example_biomarkers"
