#' Calculate Weighted Sum Scores Using Signed Protein Values
#'
#' This function calculates module scores by sign-flipping proteins based on
#' their biomarker associations (mean_beta) and weighting by alpha centrality.
#' This approach aligns proteins in the direction of core biomarker progression
#' and weights by network importance and biomarker association strength.
#'
#' @param expr_matrix A numeric matrix with samples as rows and genes/proteins
#'   as columns, already normalized and scaled (see Details). Column names
#'   should match the 'symbol' values in the result table. Proteins used by a
#'   module must not have missing values (see Details).
#' @param result A data.table containing module information with required columns:
#'   module, symbol (gene identifiers), mean_beta (for sign-flipping),
#'   mean_alpha_scaled (for weighting). A data.frame is converted to a
#'   data.table; any other class is an error. If `NULL` (the default), a
#'   prebuilt reference set shipped with the package is used instead; see
#'   `prebuilt`.
#' @param prebuilt Name of a prebuilt reference set to use when `result` is
#'   `NULL`. Defaults to `"core_AD_plasma_biomarkers"`. See
#'   [list_prebuilt_wss()] for available options. Ignored if `result` is
#'   supplied.
#' @param dedup_by Name of the column used to rank duplicate gene/module
#'   assignments when a gene appears in `result` more than once; only the
#'   highest-ranked assignment per gene is kept. Defaults to
#'   `"mean_alpha_scaled"`.
#' @param min_module_size Minimum number of genes a module must retain
#'   (after deduplication) to be included in the output. Defaults to `3`,
#'   matching the module size filter (`N > 2`) used to build the original
#'   `core_AD_plasma_biomarkers` reference set.
#'
#' @return A list containing:
#'   \item{scores}{Matrix of Weighted Sum Scores (samples x modules)}
#'   \item{summary_stats}{data.table with summary statistics for each module score}
#'   \item{module_sd}{Named vector of standard deviations for each module}
#'   \item{module_proteins}{List of protein names in each module}
#'
#'   A module with none of its proteins in `expr_matrix` cannot be scored. It
#'   triggers a warning, has all-`NA` scores and a `module_sd` of 0, and is left
#'   out of `summary_stats` and `module_proteins`.
#'
#' @details
#' For each module, proteins are sign-flipped based on sign(mean_beta) to align
#' all proteins in the direction of disease progression. The signed values are
#' then weighted by mean_alpha_scaled (network centrality + association strength)
#' and averaged to produce a module score.
#'
#' Higher module scores indicate greater module activity in the direction of
#' disease progression.
#'
#' **Input scale.** The function does no normalization or scaling of
#' `expr_matrix`; supply values that are already normalized and scaled. Each
#' module score is a weighted average of sign-aligned protein values, so the
#' scores inherit the scale of the input, whether that is a rank-based inverse
#' normal transformed (INT) scale, z-scaled log2 ratios, or another choice
#' appropriate to the platform. The bundled reference set was built from
#' SomaScan data on an INT scale. Scores computed on differently normalized or
#' scaled input are not directly comparable with each other, and
#' [predict_WSS()] predicts scores on the INT-based scale of the bundled
#' reference.
#'
#' **Missing values.** Modules are scored from the proteins that are present in
#' `expr_matrix` (see the return value for modules with no proteins present),
#' but missing values (`NA`) in a protein that a module uses are not supported:
#' the function stops with an error that names the affected proteins. Impute
#' the values, or drop the affected proteins or samples, before scoring.
#' Missing values in proteins that no module uses are ignored.
#'
#' @examples
#' expr_matrix <- matrix(
#'   c(1, 2, 3, 4, 5, 6, 7, 8),
#'   nrow = 2,
#'   dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
#' )
#' cluster_result <- data.table::data.table(
#'   module = c("mod1", "mod1", "mod1", "mod2"),
#'   symbol = c("geneA", "geneB", "geneC", "geneD"),
#'   mean_beta = c(1, -1, 1, 1),
#'   mean_alpha_scaled = c(1, 1, 2, 1)
#' )
#'
#' # Calculate Weighted Sum Scores (min_module_size = 1 to keep this toy example)
#' module_results <- calculate_WSS(expr_matrix, cluster_result, min_module_size = 1)
#'
#' # Access score matrix
#' scores <- module_results$scores
#'
#' # View summary statistics
#' print(module_results$summary_stats)
#'
#' # Check standard deviations
#' print(module_results$module_sd)
#'
#' # With no `result`, the default bundled reference set is used. The
#' # simulated example matrix has columns for all of its modules.
#' default_results <- calculate_WSS(wss_example_expression)
#' dim(default_results$scores)
#'
#' @seealso [list_prebuilt_wss()], [load_prebuilt_wss()], [load_prebuilt_terms()]
#'
#' @import data.table
#' @importFrom stats sd
#' @export
calculate_WSS <- function(expr_matrix, result = NULL, prebuilt = "core_AD_plasma_biomarkers",
                           dedup_by = "mean_alpha_scaled", min_module_size = 3) {

  if (is.null(result)) {
    result <- load_prebuilt_wss(prebuilt)
  } else {
    result <- .as_result_table(result)
  }

  result <- .filter_wss_result(result, dedup_by = dedup_by, min_module_size = min_module_size)

  # Stop early, with a clear message, rather than failing on a missing value
  .check_no_missing_expr(expr_matrix, result)

  # Get unique modules
  modules <- unique(result$module)

  # Initialize output structures
  score_matrix <- matrix(NA,
                         nrow = nrow(expr_matrix),
                         ncol = length(modules),
                         dimnames = list(rownames(expr_matrix), modules))

  summary_list <- list()
  sd_vector <- numeric(length(modules))
  names(sd_vector) <- modules
  protein_list <- list()

  # Calculate scores for each module
  for (mod in modules) {

    # Get proteins in this module
    mod_proteins <- result[module == mod]
    protein_ids <- mod_proteins$symbol

    # Check which proteins are in expression matrix
    available_proteins <- protein_ids[protein_ids %in% colnames(expr_matrix)]

    if (length(available_proteins) == 0) {
      warning(sprintf("Module %s has no proteins in expression matrix", mod))
      next
    }

    # Store protein list
    protein_list[[mod]] <- available_proteins

    # Subset expression matrix and module info
    mod_expr <- expr_matrix[, available_proteins, drop = FALSE]
    mod_info <- mod_proteins[symbol %in% available_proteins]

    # Ensure order matches
    mod_info <- mod_info[match(available_proteins, symbol)]

    # Get sign from mean_beta
    sign_vector <- sign(mod_info$mean_beta)

    # Get weights from mean_alpha_scaled
    weights <- mod_info$mean_alpha_scaled

    # Handle missing or zero weights
    if (any(is.na(weights)) || sum(weights) == 0) {
      warning(sprintf("Module %s has missing or zero weights, using equal weights", mod))
      weights <- rep(1, length(weights))
    }

    # Sign-flip proteins to align with disease progression
    signed_expr <- sweep(mod_expr, 2, sign_vector, "*")

    # Calculate weighted sum
    weighted_sum <- signed_expr %*% weights

    # Normalize by sum of weights
    module_score <- as.vector(weighted_sum / sum(weights))

    # Store in matrix
    score_matrix[, mod] <- module_score

    # Calculate summary statistics
    score_summary <- summary(module_score)
    summary_list[[mod]] <- as.list(score_summary)

    # Calculate SD
    sd_vector[mod] <- sd(module_score)
  }

  # Convert summary list to data.table
  summary_dt <- rbindlist(summary_list, idcol = "module")

  # Clean up column names if needed
  setnames(summary_dt,
           old = c("Min.", "1st Qu.", "3rd Qu.", "Max."),
           new = c("Min", "Q1", "Q3", "Max"),
           skip_absent = TRUE)

  # Return results
  return(list(
    scores = score_matrix,
    summary_stats = summary_dt,
    module_sd = sd_vector,
    module_proteins = protein_list
  ))
}
