test_that("load_prebuilt_biomarker_reference returns the bundled reference for the default set", {
  reference <- load_prebuilt_biomarker_reference("core_AD_plasma_biomarkers")

  expect_true(data.table::is.data.table(reference))
  expect_named(reference, c("biomarker", "value"))
  expect_setequal(
    unique(reference$biomarker),
    c(
      "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
      "PlasmapTau217", "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
    )
  )
  expect_true(all(is.finite(reference$value)))
})

test_that("normalize_wss_biomarkers(method = 'project') works on a single sample", {
  raw_biomarkers <- data.frame(
    PlasmaPTau181 = 1.5, PlasmaAB142P = 25, PlasmaAB140P = 300, PlasmaABRatio = 0.09,
    PlasmapTau217 = 0.3, PlasmapTau217_AB42Ratio = 0.012, PlasmaGFAP = 60, PlasmaNfL = 22
  )

  out <- normalize_wss_biomarkers(raw_biomarkers, method = "project")

  expect_equal(nrow(out), 1)
  expect_equal(sort(names(out)), sort(names(raw_biomarkers)))
  expect_true(all(vapply(out, is.finite, logical(1))))
})

test_that("normalize_wss_biomarkers(method = 'project') is deterministic and preserves row names", {
  raw_biomarkers <- data.frame(
    PlasmaPTau181 = c(1.5, 3.0),
    PlasmaNfL = c(22, 40),
    row.names = c("subject1", "subject2")
  )

  out1 <- normalize_wss_biomarkers(raw_biomarkers, method = "project")
  out2 <- normalize_wss_biomarkers(raw_biomarkers, method = "project")

  expect_equal(out1, out2)
  expect_equal(rownames(out1), c("subject1", "subject2"))
})

test_that("normalize_wss_biomarkers(method = 'self') requires multiple samples", {
  raw_biomarkers <- data.frame(PlasmaPTau181 = 1.5)
  expect_error(
    normalize_wss_biomarkers(raw_biomarkers, method = "self"),
    "requires multiple samples"
  )
})

test_that("normalize_wss_biomarkers(method = 'self') rank-normalizes within the batch", {
  skip_if_not_installed("RNOmni")

  raw_biomarkers <- data.frame(
    PlasmaPTau181 = c(1, 2, 3, 4, 5),
    PlasmaNfL = c(10, 20, 30, 40, 50)
  )

  out <- normalize_wss_biomarkers(raw_biomarkers, method = "self")

  expect_equal(nrow(out), 5)
  expect_equal(names(out), c("PlasmaPTau181", "PlasmaNfL"))
  # Monotonic increasing input should give monotonic increasing rank-norm output
  expect_true(all(diff(out$PlasmaPTau181) > 0))
})

test_that("normalize_wss_biomarkers errors on non-numeric columns", {
  raw_biomarkers <- data.frame(PlasmaPTau181 = "high", PlasmaNfL = 22)
  expect_error(
    normalize_wss_biomarkers(raw_biomarkers, method = "project"),
    "must be numeric"
  )
})

test_that("normalize_wss_biomarkers(method = 'project') errors on an unknown biomarker column", {
  raw_biomarkers <- data.frame(NotARealBiomarker = 1.5)
  expect_error(
    normalize_wss_biomarkers(raw_biomarkers, method = "project"),
    "No reference distribution"
  )
})

test_that("normalize_wss_biomarkers(method = 'project') output feeds directly into predict_WSS", {
  skip_if_not_installed("glmnet")
  skip_if_not_installed("joinet")
  skip_on_cran()

  raw_biomarkers <- data.frame(
    PlasmaPTau181 = 1.5, PlasmaAB142P = 25, PlasmaAB140P = 300, PlasmaABRatio = 0.09,
    PlasmapTau217 = 0.3, PlasmapTau217_AB42Ratio = 0.012, PlasmaGFAP = 60, PlasmaNfL = 22
  )
  normalized <- normalize_wss_biomarkers(raw_biomarkers, method = "project")

  newdata <- cbind(Age = 72, Gender = 1, normalized)
  pred <- predict_WSS(newdata)

  expect_equal(dim(pred$meta), c(1, 75))
  expect_true(all(is.finite(pred$meta)))
})
