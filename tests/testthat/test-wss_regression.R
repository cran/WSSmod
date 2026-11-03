# Lightweight synthetic regression test for the full calculate_WSS() pipeline
# (dedup -> module-size filter -> sign-flip weighted average), using default
# arguments. This mirrors, on synthetic data, an end-to-end comparison against
# the original 09-WSS_construction.R analysis output that was verified by hand
# (max abs diff ~5e-15) to confirm calculate_WSS()'s defaults (dedup_by =
# "mean_alpha_scaled", min_module_size = 3) reproduce the paper's results.

test_that("calculate_WSS() defaults match the original analysis semantics (N > 2 module filter)", {
  expect_equal(formals(calculate_WSS)$min_module_size, 3)
  expect_equal(formals(calculate_WSS)$dedup_by, "mean_alpha_scaled")
})

test_that("calculate_WSS() end-to-end pipeline matches a hand-computed synthetic example", {
  # G3 is assigned to both modA (alpha 0.8) and modC (alpha 0.3); dedup should
  # keep the modA assignment and drop modC entirely (it has no genes left).
  # modB has only 2 genes, below the default min_module_size of 3, and should
  # be dropped.
  result <- data.table::data.table(
    module = c("modA", "modA", "modA", "modC", "modB", "modB"),
    symbol = c("G1", "G2", "G3", "G3", "G4", "G5"),
    mean_beta = c(1, -1, 1, 1, 1, 1),
    mean_alpha_scaled = c(0.5, 1.0, 0.8, 0.3, 1.0, 1.0)
  )

  expr_matrix <- matrix(
    c(2, 1, 4, 3, 6, 5, 8, 7, 10, 9),
    nrow = 2,
    dimnames = list(c("sample1", "sample2"), c("G1", "G2", "G3", "G4", "G5"))
  )

  out <- calculate_WSS(expr_matrix, result)

  # Only modA survives dedup + the min_module_size filter.
  expect_equal(colnames(out$scores), "modA")
  expect_equal(sort(out$module_proteins$modA), c("G1", "G2", "G3"))

  # Hand-computed weighted average: sign-flip by mean_beta, weight by
  # mean_alpha_scaled, normalize by the sum of weights.
  #   sample1: G1=+2*0.5, G2=-4*1.0, G3=+6*0.8 -> (1 - 4 + 4.8) / 2.3
  #   sample2: G1=+1*0.5, G2=-3*1.0, G3=+5*0.8 -> (0.5 - 3 + 4) / 2.3
  expected_sample1 <- (2 * 0.5 - 4 * 1.0 + 6 * 0.8) / (0.5 + 1.0 + 0.8)
  expected_sample2 <- (1 * 0.5 - 3 * 1.0 + 5 * 0.8) / (0.5 + 1.0 + 0.8)

  expect_equal(out$scores["sample1", "modA"], expected_sample1)
  expect_equal(out$scores["sample2", "modA"], expected_sample2)
})
