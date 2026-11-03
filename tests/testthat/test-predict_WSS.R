test_that("load_prebuilt_wss_model returns the bundled joinet model for the default set", {
  model_info <- load_prebuilt_wss_model("core_AD_plasma_biomarkers")

  expect_type(model_info, "list")
  expect_true(all(c("model", "outcomes", "variant", "x_cols") %in% names(model_info)))
  expect_equal(model_info$variant, "Biomarkers_Age_Sex")
  expect_equal(
    model_info$x_cols,
    c(
      "Age", "Gender", "PlasmaPTau181", "PlasmaAB142P",
      "PlasmaAB140P", "PlasmaABRatio", "PlasmapTau217",
      "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
    )
  )
  expect_length(model_info$outcomes, 75)
})

test_that("predict_WSS errors on an unknown prebuilt name", {
  newdata <- data.frame(Age = 72, Gender = 1)
  expect_error(
    predict_WSS(newdata, prebuilt = "not_a_real_set"),
    "Unknown prebuilt WSS reference set"
  )
})

test_that("predict_WSS errors when required columns are missing", {
  newdata <- data.frame(Age = 72, Gender = 1)
  expect_error(
    predict_WSS(newdata),
    "missing required column"
  )
})

test_that("predict_WSS errors when a required column is non-numeric", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("joinet")
  skip_on_cran()

  newdata <- data.frame(
    Age = 72, Gender = "M",
    PlasmaPTau181 = 0, PlasmaAB142P = 0, PlasmaAB140P = 0,
    PlasmaABRatio = 0, PlasmapTau217 = 0,
    PlasmapTau217_AB42Ratio = 0, PlasmaGFAP = 0, PlasmaNfL = 0
  )

  expect_error(predict_WSS(newdata), "must be numeric")
})

test_that("predict_WSS returns base and meta predictions for all outcomes", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("joinet")
  skip_on_cran()

  newdata <- data.frame(
    Age = c(72, 65),
    Gender = c(1, 0),
    PlasmaPTau181 = c(0, 0.5),
    PlasmaAB142P = c(0, -0.5),
    PlasmaAB140P = c(0, 0.2),
    PlasmaABRatio = c(0, -0.2),
    PlasmapTau217 = c(0, 0.1),
    PlasmapTau217_AB42Ratio = c(0, -0.1),
    PlasmaGFAP = c(0, 0.3),
    PlasmaNfL = c(0, -0.3),
    row.names = c("subject1", "subject2")
  )

  pred <- predict_WSS(newdata)

  expect_named(pred, c("base", "meta"))
  expect_equal(dim(pred$base), c(2, 75))
  expect_equal(dim(pred$meta), c(2, 75))
  expect_equal(rownames(pred$meta), c("subject1", "subject2"))
  expect_true(all(load_prebuilt_wss_model()$outcomes %in% colnames(pred$meta)))
  expect_true(all(is.finite(pred$meta)))

  # Predictions for the two (different) synthetic subjects should differ
  expect_false(isTRUE(all.equal(pred$meta["subject1", ], pred$meta["subject2", ])))
})

test_that("predict_WSS accepts a matrix in addition to a data.frame", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("joinet")
  skip_on_cran()

  x_cols <- load_prebuilt_wss_model()$x_cols
  newdata_mat <- matrix(0, nrow = 1, ncol = length(x_cols), dimnames = list("subject1", x_cols))
  newdata_mat[, "Age"] <- 72
  newdata_mat[, "Gender"] <- 1

  pred <- predict_WSS(newdata_mat)
  expect_equal(dim(pred$meta), c(1, 75))
})

test_that("predict_WSS with RankNorm applied within a synthetic batch is deterministic", {
  # This codifies the workflow validated against the paper's held-out
  # follow-up cohort: RNOmni::RankNorm() applied to each raw biomarker
  # *within the batch being predicted on* (not against the training
  # cohort, and not left untransformed) exactly reproduced the paper's
  # published follow-up predictions (max abs diff ~5e-15). This test uses
  # synthetic biomarker values, not real subject data, so it only checks
  # the pipeline is wired correctly and deterministic, not the paper's
  # actual numbers.
  skip_if_not_installed("glmnet")
  skip_if_not_installed("joinet")
  skip_on_cran()
  skip_if_not_installed("RNOmni")

  set.seed(42)
  biomarkers <- c(
    "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
    "PlasmapTau217", "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
  )
  raw <- as.data.frame(matrix(
    rlnorm(8 * 20), nrow = 20, dimnames = list(NULL, biomarkers)
  ))

  normalized <- as.data.frame(lapply(raw, RNOmni::RankNorm))

  newdata <- cbind(
    Age = round(runif(20, 55, 90)),
    Gender = rbinom(20, 1, 0.5),
    normalized
  )

  pred1 <- predict_WSS(newdata)
  pred2 <- predict_WSS(newdata)

  expect_equal(dim(pred1$meta), c(20, 75))
  expect_true(all(is.finite(pred1$meta)))
  expect_equal(pred1$meta, pred2$meta)
})
