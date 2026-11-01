# The R port against the REAL ards-endotypes script: every value it computed (reference/python_reference.json,
# captured by tools/export_reference.py from an unmodified run that reproduced the committed outputs).

obs <- REF$observed
k <- REF$constants

test_that("n_per_arm reproduces all 14 calls the real script made", {
  expect_length(REF$calls$n_per_arm, 14)
  for (cl in REF$calls$n_per_arm) {
    a <- cl$args
    expect_equal(n_per_arm(a$p_control, a$rrr, a$alpha, a$power), cl$out, tolerance = 1e-12)
  }
})

test_that("ppv_from_accuracy reproduces all 12 calls", {
  expect_length(REF$calls$ppv_from_accuracy, 12)
  for (cl in REF$calls$ppv_from_accuracy) {
    a <- cl$args
    expect_equal(ppv_from_accuracy(a$sens, a$spec, a$prevalence), cl$out, tolerance = 1e-14)
  }
})

test_that("the accuracy curve equals L2_enrichment_curve.csv after the same rounding", {
  r <- round_curve(enrichment_curve(obs$prevalence, obs$p_hyper, obs$p_hypo, rrr = k$rrr, alpha = k$alpha,
                                    power = k$power, p_all = obs$p_all))
  py <- rows_to_df(REF$curve_rows)
  expect_equal(r, py[names(r)], tolerance = 0)
})

test_that("break-even accuracy matches, including the rounded-vs-unrounded comparison", {
  r <- round_curve(enrichment_curve(obs$prevalence, obs$p_hyper, obs$p_hypo, p_all = obs$p_all))
  expect_equal(break_even_accuracy(r, REF$n_unstrat), REF$break_even_accuracy)
})

test_that("the site table equals L2_enrichment_by_site.csv (the README's 12,494 and 647,654)", {
  r <- site_table(obs$prevalence, obs$p_hyper, obs$p_hypo, REF$scoreable, sens = k$acc_assumed,
                  spec = k$acc_assumed, rrr = k$rrr, alpha = k$alpha, power = k$power)
  expect_equal(r, rows_to_df(REF$site_rows), tolerance = 0)
  expect_equal(r$patients_to_screen[r$scenario == "restricted_setB"], 12494)
  expect_equal(r$patients_to_screen[r$scenario == "naive_setB"], 647654)
})

test_that("the site design's unrounded internals match the script's locals", {
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, k$acc_assumed, k$acc_assumed)
  expect_equal(d$ppv, REF$site$ppv, tolerance = 1e-14)
  expect_equal(d$p_enrolled, REF$site$p_enrolled, tolerance = 1e-14)
  expect_equal(d$n_per_arm, REF$site$n_per_arm, tolerance = 1e-12)
  expect_equal(d$frac_positive, REF$site$frac_positive, tolerance = 1e-14)
})

test_that("invariant: the two classes partition the cohort's overall event rate", {
  expect_equal(obs$prevalence * obs$p_hyper + (1 - obs$prevalence) * obs$p_hypo, obs$p_all, tolerance = 1e-14)
})

test_that("invariant: an uninformative classifier (accuracy 0.5) needs exactly the unstratified N", {
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, 0.5, 0.5, p_all = obs$p_all)
  expect_equal(d$n_per_arm, d$n_unstrat, tolerance = 1e-12)
})

test_that("app defaults carry the same inputs as the reference", {
  expect_identical(DEFAULTS$observed, REF$observed)
  expect_identical(DEFAULTS$scoreable, REF$scoreable)
})
