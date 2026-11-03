# Deduplicate a WSS module/weighting table and drop undersized modules.
#
# For genes assigned to more than one module, only the row with the highest
# value of `dedup_by` is kept. Modules with fewer than `min_module_size`
# remaining genes are dropped entirely.
.filter_wss_result <- function(result, dedup_by = "mean_alpha_scaled", min_module_size = 3) {
  if (!dedup_by %in% names(result)) {
    stop(sprintf(
      "dedup_by column '%s' not found in result table. Available columns: %s",
      dedup_by, paste(names(result), collapse = ", ")
    ), call. = FALSE)
  }
  if (!"module" %in% names(result) || !"symbol" %in% names(result)) {
    stop("result table must contain 'module' and 'symbol' columns", call. = FALSE)
  }

  result <- data.table::copy(result)
  data.table::setorderv(result, dedup_by, order = -1L)
  result <- result[, .SD[1], by = symbol]

  keep <- result[, .N, by = module][N >= min_module_size, module]
  result <- result[module %in% keep]

  result[]
}
