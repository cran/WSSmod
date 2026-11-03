#' Coerce a module table to a data.table
#'
#' Accepts a data.table as is and converts a data.frame (including a tibble).
#' Anything else is rejected with a message naming the class that was received.
#'
#' @param result A data.table or data.frame of module information.
#'
#' @return A data.table.
#'
#' @noRd
.as_result_table <- function(result) {
  if (data.table::is.data.table(result)) {
    return(result)
  }

  if (is.data.frame(result)) {
    return(data.table::as.data.table(result))
  }

  stop(
    sprintf(
      paste(
        "`result` must be a data.table or data.frame with columns module,",
        "symbol, mean_beta and mean_alpha_scaled, not an object of class '%s'"
      ),
      paste(class(result), collapse = "/")
    ),
    call. = FALSE
  )
}

#' Check that the proteins used by the modules have no missing values
#'
#' Only columns of `expr_matrix` that a retained module uses are checked, so
#' missing values in proteins no module uses are harmless.
#'
#' @param expr_matrix Samples-by-proteins matrix passed to [calculate_WSS()].
#' @param result Module table after filtering by `.filter_wss_result()`.
#'
#' @return `NULL`, invisibly. Stops with an informative message if any protein
#'   used by a module has a missing value.
#'
#' @noRd
.check_no_missing_expr <- function(expr_matrix, result) {
  used <- intersect(unique(result$symbol), colnames(expr_matrix))
  has_na <- colSums(is.na(expr_matrix[, used, drop = FALSE])) > 0
  na_proteins <- used[has_na]

  if (length(na_proteins) == 0) {
    return(invisible(NULL))
  }

  examples <- na_proteins[seq_len(min(5, length(na_proteins)))]
  stop(
    sprintf(
      paste(
        "`expr_matrix` has missing values (NA) in %d protein(s) used by the",
        "modules, e.g. %s. calculate_WSS() does not handle missing values;",
        "impute them, or drop the affected proteins or samples, before",
        "scoring."
      ),
      length(na_proteins),
      paste(examples, collapse = ", ")
    ),
    call. = FALSE
  )
}
