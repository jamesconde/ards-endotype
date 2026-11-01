# Patient-by-patient Monte Carlo of an enrichment trial. Independent of the closed form in enrichment.R: it never
# uses the PPV, the dilution formula or the screening ratio, so agreement between the two is evidence for both.
#
# Each screened patient: scoreable ~ Bern(scoreable); true class hyper ~ Bern(prevalence); classifier positive ~
# Bern(sens) if hyper, else Bern(1 - spec). Scoreable positives are enrolled until 2n, then randomised 1:1.
# Events ~ Bern(p_class) in control and Bern(p_class * (1 - rrr)) in treatment. Analysis: the two-sided pooled
# two-proportion z-test, the test the sample-size formula is derived for.

simulate_trials <- function(n_per_arm, prevalence, p_hyper, p_hypo, sens, spec, scoreable = 1, rrr = RRR,
                            alpha = ALPHA, n_sims = 1000, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- ceiling(n_per_arm)                    # whole patients; at most +1 per arm vs the formula
  need <- 2L * n
  # Only sizes each screening batch; a wrong guess costs extra loops, never a different answer.
  q_guess <- max(scoreable * (sens * prevalence + (1 - spec) * (1 - prevalence)), 1e-6)
  out <- matrix(NA_real_, n_sims, 4, dimnames = list(NULL, c("screened", "e0", "e1", "hyper_frac")))
  for (s in seq_len(n_sims)) {
    screened <- 0
    cls <- logical(0)
    while (length(cls) < need) {
      m <- max(1000L, ceiling(1.1 * (need - length(cls)) / q_guess))
      hyper <- runif(m) < prevalence
      pos <- runif(m) < hyper * sens + (!hyper) * (1 - spec)    # P(positive | true class)
      if (scoreable < 1) pos <- pos & runif(m) < scoreable
      idx <- which(pos)
      take <- min(length(idx), need - length(cls))
      if (take > 0) cls <- c(cls, hyper[idx[seq_len(take)]])
      screened <- screened + if (length(cls) == need) idx[take] else m
    }
    arm <- sample(rep(c(FALSE, TRUE), each = n))          # TRUE = treatment
    p <- (p_hypo + cls * (p_hyper - p_hypo)) * (1 - rrr * arm)
    ev <- runif(need) < p
    out[s, ] <- c(screened, sum(ev[!arm]), sum(ev[arm]), mean(cls))
  }
  pbar <- (out[, "e0"] + out[, "e1"]) / (2 * n)
  z <- (out[, "e0"] - out[, "e1"]) / n / sqrt(pbar * (1 - pbar) * 2 / n)
  reject <- abs(z) > qnorm(1 - alpha / 2)
  reject[is.na(reject)] <- FALSE                          # no events in either arm: z undefined, not rejected
  pw <- mean(reject)
  se <- sqrt(pw * (1 - pw) / n_sims)
  list(n_per_arm_simulated = n, n_sims = n_sims, power = pw, power_se = se,
       power_ci = c(max(0, pw - 1.96 * se), min(1, pw + 1.96 * se)),
       screen_per_enrolled = mean(out[, "screened"]) / need,
       ppv = mean(out[, "hyper_frac"]), control_event_rate = mean(out[, "e0"]) / n)
}

# Expected screened patients per simulated trial, used only to cap the app's run time.
expected_screened <- function(n_per_arm, prevalence, sens, spec, scoreable) {
  2 * ceiling(n_per_arm) / max(scoreable * (sens * prevalence + (1 - spec) * (1 - prevalence)), 1e-9)
}
