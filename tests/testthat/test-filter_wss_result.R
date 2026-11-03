test_that(".filter_wss_result keeps the highest-ranked duplicate per symbol", {
  result <- data.table::data.table(
    module = c("mod1", "mod2"),
    symbol = c("geneA", "geneA"),
    mean_beta = c(1, 1),
    mean_alpha_scaled = c(1, 5)
  )

  out <- .filter_wss_result(result, min_module_size = 1)

  expect_equal(nrow(out), 1)
  expect_equal(out$module, "mod2")
})

test_that(".filter_wss_result can rank duplicates by a different column", {
  result <- data.table::data.table(
    module = c("mod1", "mod2"),
    symbol = c("geneA", "geneA"),
    mean_beta = c(1, 1),
    mean_alpha_scaled = c(5, 1),
    custom_rank = c(1, 9)
  )

  out <- .filter_wss_result(result, dedup_by = "custom_rank", min_module_size = 1)

  expect_equal(out$module, "mod2")
})

test_that(".filter_wss_result drops modules below min_module_size", {
  result <- data.table::data.table(
    module = c("mod1", "mod1", "mod1", "mod2"),
    symbol = c("geneA", "geneB", "geneC", "geneD"),
    mean_beta = c(1, 1, 1, 1),
    mean_alpha_scaled = c(1, 1, 1, 1)
  )

  out <- .filter_wss_result(result, min_module_size = 3)

  expect_equal(unique(out$module), "mod1")
  expect_equal(nrow(out), 3)
})

test_that(".filter_wss_result errors when dedup_by column is missing", {
  result <- data.table::data.table(
    module = "mod1",
    symbol = "geneA",
    mean_beta = 1,
    mean_alpha_scaled = 1
  )

  expect_error(
    .filter_wss_result(result, dedup_by = "not_a_column"),
    "not found in result table"
  )
})

test_that(".filter_wss_result errors when required columns are missing", {
  result <- data.table::data.table(
    symbol = "geneA",
    mean_alpha_scaled = 1
  )

  expect_error(
    .filter_wss_result(result),
    "must contain 'module' and 'symbol'"
  )
})

test_that(".filter_wss_result does not mutate the caller's table", {
  result <- data.table::data.table(
    module = c("mod1", "mod1"),
    symbol = c("geneA", "geneB"),
    mean_beta = c(1, 1),
    mean_alpha_scaled = c(2, 1)
  )
  original_order <- data.table::copy(result)

  invisible(.filter_wss_result(result, min_module_size = 1))

  expect_equal(result, original_order)
})
