test_that("calculate_WSS computes expected weighted averages", {
  set.seed(1)
  expr_matrix <- matrix(
    c(1, 2, 3, 4, 5, 6, 7, 8),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
  )

  result <- data.table::data.table(
    module = c("mod1", "mod1", "mod1", "mod2"),
    symbol = c("geneA", "geneB", "geneC", "geneD"),
    mean_beta = c(1, -1, 1, 1),
    mean_alpha_scaled = c(1, 1, 2, 1)
  )

  out <- calculate_WSS(expr_matrix, result, min_module_size = 1)

  expect_named(out, c("scores", "summary_stats", "module_sd", "module_proteins"))
  expect_equal(colnames(out$scores), c("mod1", "mod2"))

  # sample1: geneA=1, geneB=3 -> -3 (flipped), geneC=5, weights 1,1,2 -> (1 - 3 + 5*2)/4 = 2
  expect_equal(out$scores["sample1", "mod1"], 2)
  expect_equal(out$scores["sample1", "mod2"], 7)

  expect_equal(sort(out$module_proteins$mod1), c("geneA", "geneB", "geneC"))
  expect_equal(out$module_proteins$mod2, "geneD")
})

test_that("calculate_WSS falls back to equal weights when all weights are zero", {
  expr_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB"))
  )

  result <- data.table::data.table(
    module = c("mod1", "mod1"),
    symbol = c("geneA", "geneB"),
    mean_beta = c(1, 1),
    mean_alpha_scaled = c(0, 0)
  )

  expect_warning(
    out <- calculate_WSS(expr_matrix, result, min_module_size = 1),
    "missing or zero weights"
  )

  expect_equal(out$scores["sample1", "mod1"], mean(c(1, 3)))
})

test_that("calculate_WSS warns when a module has no matching proteins", {
  expr_matrix <- matrix(
    c(1, 2),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), "geneA")
  )

  result <- data.table::data.table(
    module = "mod1",
    symbol = "geneMissing",
    mean_beta = 1,
    mean_alpha_scaled = 1
  )

  expect_warning(
    calculate_WSS(expr_matrix, result, min_module_size = 1),
    "no proteins in expression matrix"
  )
})

test_that("calculate_WSS uses the default prebuilt reference set when result is NULL", {
  # HYOU1 and DNAJB11 are both members of module Merged.M7 in the default
  # core_AD_plasma_biomarkers reference set.
  expr_matrix <- matrix(
    c(1, 2, 3, 4),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("HYOU1", "DNAJB11"))
  )

  # Every other module in the reference set will warn about missing
  # proteins, since expr_matrix only supplies genes for one module.
  out <- suppressWarnings(calculate_WSS(expr_matrix))

  expect_true("Merged.M7" %in% colnames(out$scores))
  expect_true(all(c("HYOU1", "DNAJB11") %in% out$module_proteins$Merged.M7))
})

test_that("calculate_WSS errors on an unknown prebuilt name", {
  expr_matrix <- matrix(1, nrow = 1, dimnames = list("sample1", "geneA"))
  expect_error(
    calculate_WSS(expr_matrix, prebuilt = "not_a_real_set"),
    "Unknown prebuilt WSS reference set"
  )
})

test_that("calculate_WSS drops modules smaller than min_module_size", {
  expr_matrix <- matrix(
    1:6,
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC"))
  )

  result <- data.table::data.table(
    module = c("mod1", "mod1", "mod2"),
    symbol = c("geneA", "geneB", "geneC"),
    mean_beta = c(1, 1, 1),
    mean_alpha_scaled = c(1, 1, 1)
  )

  # Default min_module_size = 3 drops both the 2-gene and 1-gene modules
  out_default <- calculate_WSS(expr_matrix, result)
  expect_equal(ncol(out_default$scores), 0)

  # min_module_size = 2 keeps mod1 (2 genes) but drops mod2 (1 gene)
  out_relaxed <- calculate_WSS(expr_matrix, result, min_module_size = 2)
  expect_equal(colnames(out_relaxed$scores), "mod1")
})

test_that("calculate_WSS dedups genes assigned to multiple modules by dedup_by", {
  expr_matrix <- matrix(
    1:8,
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("geneA", "geneB", "geneC", "geneD"))
  )

  # geneA is assigned to both mod1 and mod2; mod2's assignment has the
  # higher mean_alpha_scaled and should win.
  result <- data.table::data.table(
    module = c("mod1", "mod1", "mod2", "mod2", "mod2"),
    symbol = c("geneA", "geneB", "geneA", "geneC", "geneD"),
    mean_beta = c(1, 1, 1, 1, 1),
    mean_alpha_scaled = c(1, 1, 5, 1, 1)
  )

  out <- calculate_WSS(expr_matrix, result, min_module_size = 1)

  expect_equal(sort(colnames(out$scores)), c("mod1", "mod2"))
  expect_equal(out$module_proteins$mod1, "geneB")
  expect_equal(sort(out$module_proteins$mod2), c("geneA", "geneC", "geneD"))
})

test_that("calculate_WSS errors when dedup_by column is missing", {
  expr_matrix <- matrix(1, nrow = 1, dimnames = list("sample1", "geneA"))
  result <- data.table::data.table(
    module = "mod1",
    symbol = "geneA",
    mean_beta = 1,
    mean_alpha_scaled = 1
  )

  expect_error(
    calculate_WSS(expr_matrix, result, dedup_by = "not_a_column"),
    "not found in result table"
  )
})
