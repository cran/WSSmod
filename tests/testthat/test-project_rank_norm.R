test_that("project_rank_norm is deterministic and monotonic", {
  set.seed(1)
  train_vals <- rnorm(200)

  new_vals <- c(-3, -1, 0, 1, 3)
  out1 <- project_rank_norm(train_vals, new_vals)
  out2 <- project_rank_norm(train_vals, new_vals)

  expect_equal(out1, out2)
  expect_true(all(diff(out1) > 0))
})

test_that("project_rank_norm handles a single new value", {
  set.seed(2)
  train_vals <- rnorm(100)

  out <- project_rank_norm(train_vals, 0.5)
  expect_length(out, 1)
  expect_true(is.finite(out))
})

test_that("project_rank_norm places a value far below the reference near the low tail", {
  train_vals <- 1:100
  out_low <- project_rank_norm(train_vals, -1000)
  out_high <- project_rank_norm(train_vals, 1000)

  expect_lt(out_low, -2)
  expect_gt(out_high, 2)
})

test_that("project_rank_norm assigns tied reference values the same projected quantile", {
  train_vals <- c(1, 2, 2, 2, 3, 4, 5)
  out <- project_rank_norm(train_vals, c(2, 2, 2))
  expect_equal(out[1], out[2])
  expect_equal(out[2], out[3])
})

test_that("project_rank_norm errors on NA in train_vals", {
  expect_error(
    project_rank_norm(c(1, 2, NA, 4), 1),
    "must not contain NA"
  )
})
