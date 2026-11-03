# A small custom module table over proteins that exist in the example matrix
proteins <- colnames(wss_example_expression)[1:10]

make_custom <- function() {
  data.table::data.table(
    module = c(rep("modA", 5), rep("modB", 5)),
    symbol = proteins,
    mean_beta = c(rep(1, 5), rep(-1, 5)),
    mean_alpha_scaled = rep(c(0.9, 0.8, 0.7, 0.6, 0.5), 2)
  )
}

test_that("a data.frame result gives the same scores as a data.table", {
  custom <- make_custom()

  from_table <- calculate_WSS(wss_example_expression, custom)
  from_frame <- calculate_WSS(wss_example_expression, as.data.frame(custom))

  expect_equal(from_frame$scores, from_table$scores)
  expect_equal(from_frame$module_proteins, from_table$module_proteins)
})

test_that("a result that is not a table is rejected with an informative error", {
  expect_error(
    calculate_WSS(wss_example_expression, as.matrix(make_custom())),
    "must be a data.table or data.frame.*'matrix/array'"
  )
  expect_error(
    calculate_WSS(wss_example_expression, as.list(make_custom())),
    "must be a data.table or data.frame.*'list'"
  )
})

test_that("a missing value in a protein a module uses stops with a clear error", {
  expr <- wss_example_expression
  expr[1, proteins[2]] <- NA

  expect_error(
    calculate_WSS(expr, make_custom()),
    "`expr_matrix` has missing values.*1 protein"
  )
  expect_error(calculate_WSS(expr, make_custom()), proteins[2], fixed = TRUE)
})

test_that("the default reference set also stops on a missing value", {
  expr <- wss_example_expression
  expr[3, 5] <- NA

  expect_error(calculate_WSS(expr), "`expr_matrix` has missing values")
})

test_that("a missing value in a protein no module uses is ignored", {
  expr <- wss_example_expression
  unused <- colnames(expr)[20]
  expr[1, unused] <- NA

  res <- calculate_WSS(expr, make_custom())

  expect_true(all(is.finite(res$scores)))
})

test_that("a module with no proteins in the matrix warns and is documented as 0 sd", {
  custom <- make_custom()
  custom[module == "modB", symbol := paste0("NOT_A_PROTEIN", 1:5)]

  expect_warning(
    res <- calculate_WSS(wss_example_expression, custom),
    "Module modB has no proteins in expression matrix"
  )

  expect_true(all(is.na(res$scores[, "modB"])))
  expect_equal(unname(res$module_sd["modB"]), 0)
  expect_false("modB" %in% names(res$module_proteins))
  expect_false("modB" %in% res$summary_stats$module)
  expect_true(all(is.finite(res$scores[, "modA"])))
})
