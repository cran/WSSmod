#' Project New Values onto a Frozen Rank-Normal Reference
#'
#' Projects `new_vals` onto the rank-based inverse-normal transform defined
#' by `train_vals`, without ever refitting or observing `new_vals` when
#' determining the transform. Each new value is treated as a single
#' insertion into the frozen, sorted `train_vals` reference: its rank is the
#' count of reference values at or below it (ties in `train_vals` share the
#' same rank), and that rank is converted to a normal quantile using the
#' same Blom-type constant (`k = 0.375`) as [RNOmni::RankNorm()].
#'
#' This is useful for scoring new samples (including a single new sample)
#' against a fixed reference cohort without leaking information from the
#' new samples into the reference's distribution -- e.g. avoiding the
#' train/test leakage that results from rank-normalizing a pooled
#' reference+new-sample dataset together.
#'
#' Note this does **not** reproduce `RNOmni::RankNorm(train_vals)`
#' bit-for-bit when `new_vals` happens to equal `train_vals`: it always
#' treats each new value as an *additional* observation inserted into the
#' `n`-sample reference (using `n + 1` in the denominator), rather than as
#' one of the original `n` samples. The two are close but not identical.
#'
#' @param train_vals Numeric vector, the frozen reference distribution
#'   (e.g. a baseline cohort's biomarker values). Must not contain `NA`.
#' @param new_vals Numeric vector of new values to project onto
#'   `train_vals`. May be a single value.
#' @param k Blom-type constant, passed through to match
#'   [RNOmni::RankNorm()]'s default.
#'
#' @return A numeric vector the same length as `new_vals`, giving each
#'   value's projected normal quantile.
#'
#' @seealso [normalize_wss_biomarkers()]
#'
#' @examples
#' train_vals <- rnorm(200)
#' project_rank_norm(train_vals, new_vals = c(-1, 0, 1))
#'
#' @export
project_rank_norm <- function(train_vals, new_vals, k = 0.375) {
  if (anyNA(train_vals)) {
    stop("train_vals must not contain NA", call. = FALSE)
  }

  n <- length(train_vals)
  sorted_train <- sort(train_vals)
  r_new <- findInterval(new_vals, sorted_train) + 1
  denom <- (n + 1) - 2 * k + 1
  stats::qnorm((r_new - k) / denom)
}
