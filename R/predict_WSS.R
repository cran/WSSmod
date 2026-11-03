#' Predict WSS Module Scores from Core Biomarkers
#'
#' Predicts module scores for a prebuilt reference set's `joinet` model (a
#' two-layer elastic net stack) without needing to measure the full
#' proteomic panel. For `"core_AD_plasma_biomarkers"`, this predicts all 75
#' module scores from Age, Gender, and the 8 core plasma biomarkers used in
#' the original analysis.
#'
#' @param newdata A data.frame or matrix with one row per sample and (at
#'   least) the columns required by the model; see `x_cols` in
#'   [load_prebuilt_wss_model()] for the exact set and order. For
#'   `"core_AD_plasma_biomarkers"` this is `Age` (raw years), `Gender`
#'   (numeric, `1` = male, `0` = not male), and the 8 core biomarkers
#'   (`PlasmaPTau181`, `PlasmaAB142P`, `PlasmaAB140P`, `PlasmaABRatio`,
#'   `PlasmapTau217`, `PlasmapTau217_AB42Ratio`, `PlasmaGFAP`, `PlasmaNfL`).
#'
#'   **The biomarker columns must already be rank-based inverse-normal
#'   transformed before calling this function** -- see
#'   [normalize_wss_biomarkers()] to build them from raw biomarker values,
#'   including for a single new patient.
#' @param prebuilt Name of a prebuilt reference set. See
#'   [list_prebuilt_wss()] for available options.
#' @param type Prediction type passed to `predict.joinet()`: `"response"`
#'   (the default) or `"link"`.
#'
#' @return A list with components `base` (first-layer-only predictions) and
#'   `meta` (final stacked predictions), each a matrix with one row per
#'   sample (row names taken from `newdata`, if present) and one column per
#'   module (named by the model's `outcomes`).
#'
#' @seealso [normalize_wss_biomarkers()], [load_prebuilt_wss_model()],
#'   [calculate_WSS()]
#'
#' @examples
#' if (requireNamespace("glmnet", quietly = TRUE) &&
#'     requireNamespace("joinet", quietly = TRUE)) {
#'   # Simulated raw biomarkers; works for a single patient too
#'   raw_biomarkers <- wss_example_biomarkers[, -(1:2)]
#'   normalized <- normalize_wss_biomarkers(raw_biomarkers)
#'
#'   newdata <- cbind(
#'     Age = wss_example_biomarkers$Age,
#'     Gender = wss_example_biomarkers$Gender,
#'     normalized
#'   )
#'   pred <- predict_WSS(newdata)
#'   pred$meta[1:3, 1:4]
#' }
#'
#' @export
predict_WSS <- function(newdata, prebuilt = "core_AD_plasma_biomarkers", type = "response") {
  model_info <- load_prebuilt_wss_model(prebuilt)
  if (is.null(model_info)) {
    stop(sprintf(
      "Prebuilt reference set '%s' has no associated prediction model", prebuilt
    ), call. = FALSE)
  }

  missing_cols <- setdiff(model_info$x_cols, colnames(newdata))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "newdata is missing required column(s): %s",
      paste(missing_cols, collapse = ", ")
    ), call. = FALSE)
  }

  newdata <- newdata[, model_info$x_cols, drop = FALSE]
  non_numeric <- names(newdata)[!vapply(newdata, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(sprintf(
      "newdata column(s) must be numeric: %s",
      paste(non_numeric, collapse = ", ")
    ), call. = FALSE)
  }

  for (pkg in c("glmnet", "joinet")) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(sprintf(
        "Package '%s' is required for predict_WSS() but is not installed", pkg
      ), call. = FALSE)
    }
    loadNamespace(pkg)
  }

  newx <- as.matrix(newdata)
  pred <- stats::predict(model_info$model, newx = newx, type = type)

  colnames(pred$base) <- colnames(pred$meta) <- model_info$outcomes
  rownames(pred$base) <- rownames(pred$meta) <- rownames(newx)

  pred
}
