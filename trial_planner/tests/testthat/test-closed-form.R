# n_per_arm against base R's power.prop.test, an independent implementation (it solves the power equation
# numerically with uniroot instead of using the closed form).

test_that("n_per_arm agrees with stats::power.prop.test over a 128-design grid", {
  g <- expand.grid(p = c(0.02, 0.05, 0.10, 0.14, 0.19, 0.30, 0.50, 0.80), rrr = c(0.10, 0.20, 0.35, 0.50),
                   alpha = c(0.01, 0.05), power = c(0.80, 0.90))
  for (i in seq_len(nrow(g))) {
    p <- g$p[i]; rrr <- g$rrr[i]; alpha <- g$alpha[i]; pw <- g$power[i]
    n <- n_per_arm(p, rrr, alpha, pw)
    ppt <- power.prop.test(p1 = p, p2 = p * (1 - rrr), sig.level = alpha, power = pw)
    expect_lt(abs(n - ppt$n), 1e-3)
    back <- power.prop.test(n = n, p1 = p, p2 = p * (1 - rrr), sig.level = alpha)$power
    expect_equal(back, pw, tolerance = 1e-10)
  }
})

test_that("degenerate designs return Inf / NaN like the Python", {
  expect_identical(n_per_arm(0), Inf)
  expect_identical(n_per_arm(1), Inf)
  expect_identical(n_per_arm(0.2, rrr = 0), Inf)
  expect_true(is.nan(ppv_from_accuracy(0, 1, 0.3)))
})
