test_that("wss_example_expression is a numeric samples-by-proteins matrix", {
  expect_true(is.matrix(wss_example_expression))
  expect_true(is.numeric(wss_example_expression))
  expect_equal(dim(wss_example_expression), c(50, 887))
  expect_false(anyNA(wss_example_expression))
  expect_false(is.null(rownames(wss_example_expression)))
  expect_false(is.null(colnames(wss_example_expression)))
})

test_that("calculate_WSS scores every default module from wss_example_expression", {
  expect_no_warning(res <- calculate_WSS(wss_example_expression))

  expected_modules <- unique(.filter_wss_result(load_prebuilt_wss())$module)
  expect_setequal(colnames(res$scores), expected_modules)
  expect_equal(nrow(res$scores), nrow(wss_example_expression))
  expect_true(all(is.finite(res$scores)))
})

test_that("wss_example_biomarkers has the columns normalize_wss_biomarkers/predict_WSS need", {
  biomarkers <- c(
    "PlasmaPTau181", "PlasmaAB142P", "PlasmaAB140P", "PlasmaABRatio",
    "PlasmapTau217", "PlasmapTau217_AB42Ratio", "PlasmaGFAP", "PlasmaNfL"
  )

  expect_s3_class(wss_example_biomarkers, "data.frame")
  expect_equal(nrow(wss_example_biomarkers), 20)
  expect_named(wss_example_biomarkers, c("Age", "Gender", biomarkers))
  expect_true(all(wss_example_biomarkers$Gender %in% c(0, 1)))
  expect_false(anyNA(wss_example_biomarkers))

  normalized <- normalize_wss_biomarkers(wss_example_biomarkers[, biomarkers])
  expect_equal(dim(normalized), c(20, 8))
  expect_true(all(is.finite(as.matrix(normalized))))
})
