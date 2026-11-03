#' Normalize Raw Biomarker Values for predict_WSS()
#'
#' Builds the rank-based inverse-normal-transformed biomarker columns
#' [predict_WSS()] requires, from raw biomarker values, using one of two
#' methods. Output columns keep the same names as the input columns (e.g.
#' `PlasmaPTau181` stays `PlasmaPTau181`) -- [predict_WSS()] reads whatever
#' column names its underlying model actually expects.
#'
#' @param raw_biomarkers A data.frame or matrix, one row per sample, with
#'   columns named by the raw (unsuffixed) biomarker names (e.g.
#'   `PlasmaPTau181`, not `PlasmaPTau181_norm`).
#' @param prebuilt Name of a prebuilt reference set to use when
#'   `method = "project"`. See [list_prebuilt_wss()] for available options.
#'   Ignored when `method = "self"`.
#' @param method Normalization method:
#'   \describe{
#'     \item{`"project"`}{(default) Project each raw value onto a bundled
#'       reference distribution via [project_rank_norm()]. Works for a
#'       single sample or a small/differently-distributed cohort, since
#'       each value is scored independently against the fixed reference.}
#'     \item{`"self"`}{Rank-normalize (via [RNOmni::RankNorm()]) within
#'       `raw_biomarkers` itself. Appropriate if you have a sizeable cohort
#'       of your own with a distribution you're comfortable normalizing
#'       against directly, instead of the bundled reference. Requires at
#'       least two samples, since a rank transform needs multiple values to
#'       rank against.}
#'   }
#'
#' @return A data.frame with one normalized column per input biomarker
#'   column (same names as `raw_biomarkers`), row names preserved. Combine
#'   with `Age`/`Gender` columns and pass to [predict_WSS()].
#'
#' @seealso [predict_WSS()], [project_rank_norm()],
#'   [load_prebuilt_biomarker_reference()]
#'
#' @examples
#' raw_biomarkers <- data.frame(PlasmaPTau181 = 1.5, PlasmaNfL = 20)
#' normalize_wss_biomarkers(raw_biomarkers)
#'
#' @export
normalize_wss_biomarkers <- function(raw_biomarkers, prebuilt = "core_AD_plasma_biomarkers",
                                      method = c("project", "self")) {
  method <- match.arg(method)

  raw_biomarkers <- as.data.frame(raw_biomarkers)
  biomarkers <- names(raw_biomarkers)
  non_numeric <- biomarkers[!vapply(raw_biomarkers, is.numeric, logical(1))]
  if (length(non_numeric) > 0) {
    stop(sprintf(
      "raw_biomarkers column(s) must be numeric: %s",
      paste(non_numeric, collapse = ", ")
    ), call. = FALSE)
  }

  if (method == "self") {
    if (nrow(raw_biomarkers) < 2) {
      stop(
        "method = \"self\" requires multiple samples to rank against; ",
        "use method = \"project\" for a single sample or small cohort.",
        call. = FALSE
      )
    }
    if (!requireNamespace("RNOmni", quietly = TRUE)) {
      stop("Package 'RNOmni' is required for method = \"self\" but is not installed", call. = FALSE)
    }
    normalized <- lapply(raw_biomarkers, RNOmni::RankNorm)
  } else {
    reference <- load_prebuilt_biomarker_reference(prebuilt)
    if (is.null(reference)) {
      stop(sprintf(
        "Prebuilt reference set '%s' has no associated biomarker reference", prebuilt
      ), call. = FALSE)
    }

    missing_biomarkers <- setdiff(biomarkers, unique(reference$biomarker))
    if (length(missing_biomarkers) > 0) {
      stop(sprintf(
        "No reference distribution for biomarker(s): %s",
        paste(missing_biomarkers, collapse = ", ")
      ), call. = FALSE)
    }

    normalized <- lapply(biomarkers, function(b) {
      train_vals <- reference$value[reference$biomarker == b]
      project_rank_norm(train_vals, raw_biomarkers[[b]])
    })
    names(normalized) <- biomarkers
  }

  normalized <- as.data.frame(normalized)
  names(normalized) <- biomarkers
  rownames(normalized) <- rownames(raw_biomarkers)

  normalized
}
