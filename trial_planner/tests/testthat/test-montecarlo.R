# The closed form against patient-by-patient simulation (app/R/simulate.R), which never uses the PPV, dilution or
# screening formulas. Tolerances are 3.5 Monte Carlo standard errors; seeds fix the draws.

obs <- REF$observed
mc_se <- function(p, n) sqrt(p * (1 - p) / n)

test_that("power at the formula's N is the nominal 80% (0.85 classifier, full coverage)", {
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, 0.85, 0.85)
  sim <- simulate_trials(d$n_per_arm, obs$prevalence, obs$p_hyper, obs$p_hypo, 0.85, 0.85, n_sims = 2000, seed = 1)
  expect_lt(abs(sim$power - POWER), 3.5 * mc_se(POWER, 2000))
  expect_lt(abs(sim$screen_per_enrolled / d$screen_per_enrolled - 1), 0.01)
  expect_lt(abs(sim$ppv - d$ppv), 0.002)
  expect_lt(abs(sim$control_event_rate - d$p_enrolled), 3.5 * mc_se(d$p_enrolled, 2000 * ceiling(d$n_per_arm)))
})

test_that("asymmetric classifier, partial coverage, RRR 25%: power and screening still match", {
  sens <- 0.92; spec <- 0.70; fs <- REF$scoreable$restricted_setB; rrr <- 0.25
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, sens, spec, fs, rrr = rrr)
  sim <- simulate_trials(d$n_per_arm, obs$prevalence, obs$p_hyper, obs$p_hypo, sens, spec, fs, rrr = rrr,
                         n_sims = 1000, seed = 2)
  expect_lt(abs(sim$power - POWER), 3.5 * mc_se(POWER, 1000))
  expect_lt(abs(sim$screen_per_enrolled / d$screen_per_enrolled - 1), 0.01)
})

test_that("negative control: sizing for PERFECT enrichment with a 0.70 classifier is visibly underpowered", {
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, 0.70, 0.70)
  n_perfect <- n_per_arm(obs$p_hyper)
  expected <- power.prop.test(n = n_perfect, p1 = d$p_enrolled, p2 = d$p_enrolled * (1 - RRR))$power
  sim <- simulate_trials(n_perfect, obs$prevalence, obs$p_hyper, obs$p_hypo, 0.70, 0.70, n_sims = 1000, seed = 3)
  expect_lt(sim$power, POWER - 3 * mc_se(POWER, 1000))
  expect_lt(abs(sim$power - expected), 3.5 * mc_se(expected, 1000))
})

test_that("with no treatment effect the rejection rate is alpha", {
  d <- design(obs$prevalence, obs$p_hyper, obs$p_hypo, 0.85, 0.85)
  sim <- simulate_trials(d$n_per_arm, obs$prevalence, obs$p_hyper, obs$p_hypo, 0.85, 0.85, rrr = 0,
                         n_sims = 1500, seed = 4)
  expect_lt(abs(sim$power - ALPHA), 3.5 * mc_se(ALPHA, 1500))
})
