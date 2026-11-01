# Closed-form prognostic-enrichment calculations.
#
# Port of ards-endotypes/layer2_clinical/09_enrichment_simulation.py (commit d5f5f66). Names, argument order and
# rounding mirror the Python; tests/testthat/test-parity.R checks every value against the real script's own
# outputs (reference/python_reference.json). The app generalises one thing: sensitivity and specificity can
# differ (the Python curve sets both to "accuracy").
#
# PROGNOSTIC enrichment only (which patients have more events), never predictive enrichment (who benefits more).

ALPHA <- 0.05
POWER <- 0.80
RRR <- 0.20
ACC_GRID <- c(1.00, 0.95, 0.90, 0.85, 0.80, 0.75, 0.70, 0.65, 0.60, 0.55, 0.50)

# Two-proportion sample size per arm, no continuity correction. Python: n_per_arm().
n_per_arm <- function(p_control, rrr = RRR, alpha = ALPHA, power = POWER) {
  p1 <- p_control
  p2 <- p_control * (1 - rrr)
  if (p1 <= 0 || p1 >= 1 || p1 == p2) return(Inf)
  pbar <- (p1 + p2) / 2
  z_a <- qnorm(1 - alpha / 2)
  z_b <- qnorm(power)
  num <- (z_a * sqrt(2 * pbar * (1 - pbar)) + z_b * sqrt(p1 * (1 - p1) + p2 * (1 - p2)))^2
  num / (p1 - p2)^2
}

# Fraction of classifier positives that truly are hyper. Python: ppv_from_accuracy().
ppv_from_accuracy <- function(sens, spec, prevalence) {
  tp <- sens * prevalence
  fp <- (1 - spec) * (1 - prevalence)
  if ((tp + fp) > 0) tp / (tp + fp) else NaN
}

# Fraction of screened (scoreable) patients the classifier calls positive.
frac_positive <- function(sens, spec, prevalence) sens * prevalence + (1 - spec) * (1 - prevalence)

# Event rate among enrolled patients: misclassification dilutes it toward the hypo rate.
enrolled_event_rate <- function(ppv, p_hyper, p_hypo) ppv * p_hyper + (1 - ppv) * p_hypo

# Everything for one design. scoreable = fraction of screened patients the classifier can score at the site
# (a patient who cannot be scored cannot be enrolled). p_all defaults to the law-of-total-probability mixture,
# which equals the observed overall rate on the study cohort (checked in test-parity.R).
design <- function(prevalence, p_hyper, p_hypo, sens, spec, scoreable = 1, rrr = RRR, alpha = ALPHA,
                   power = POWER, p_all = prevalence * p_hyper + (1 - prevalence) * p_hypo) {
  ppv <- ppv_from_accuracy(sens, spec, prevalence)
  p_enrolled <- enrolled_event_rate(ppv, p_hyper, p_hypo)
  n <- n_per_arm(p_enrolled, rrr, alpha, power)
  n_unstrat <- n_per_arm(p_all, rrr, alpha, power)
  fpos <- frac_positive(sens, spec, prevalence)
  eff <- fpos * scoreable
  list(ppv = ppv, p_enrolled = p_enrolled, n_per_arm = n, n_total = 2 * n, n_unstrat = n_unstrat,
       vs_unstratified = n / n_unstrat, frac_positive = fpos, scoreable = scoreable,
       screen_per_enrolled = if (eff > 0) 1 / eff else Inf,
       patients_to_screen = if (eff > 0) 2 * n / eff else Inf)
}

# Required N vs classifier accuracy, sens = spec = accuracy (the Python loop). Unrounded.
enrichment_curve <- function(prevalence, p_hyper, p_hypo, acc = ACC_GRID, rrr = RRR, alpha = ALPHA,
                             power = POWER, p_all = prevalence * p_hyper + (1 - prevalence) * p_hypo) {
  rows <- lapply(acc, function(a) {
    d <- design(prevalence, p_hyper, p_hypo, a, a, 1, rrr, alpha, power, p_all)
    data.frame(classifier_accuracy = a, ppv = d$ppv, event_rate_enrolled = d$p_enrolled,
               n_per_arm = d$n_per_arm, n_total = d$n_total, screen_per_enrolled = d$screen_per_enrolled,
               screened_total = d$n_total * d$screen_per_enrolled, vs_unstratified = d$vs_unstratified)
  })
  do.call(rbind, rows)
}

# The rounding the Python applies before writing L2_enrichment_curve.csv.
round_curve <- function(curve) {
  transform(curve, ppv = round(ppv, 4), event_rate_enrolled = round(event_rate_enrolled, 4),
            n_per_arm = round(n_per_arm), n_total = round(n_total),
            screen_per_enrolled = round(screen_per_enrolled, 2), screened_total = round(screened_total),
            vs_unstratified = round(vs_unstratified, 3))
}

# Largest accuracy whose ROUNDED N/arm reaches the UNROUNDED unstratified N, as the Python computes it.
# #Consult with James: at accuracy 0.5 the two N are equal in exact arithmetic (an uninformative classifier
# enrols the base rate), so this returns 0.5 or NA depending only on which way n_unstrat rounds.
break_even_accuracy <- function(curve_rounded, n_unstrat) {
  hit <- curve_rounded$classifier_accuracy[curve_rounded$n_per_arm >= n_unstrat]
  if (length(hit)) max(hit) else NA_real_
}

# Screening feasibility by site at one assumed accuracy (Python: acc_assumed = 0.85), rounded like
# L2_enrichment_by_site.csv.
site_table <- function(prevalence, p_hyper, p_hypo, scoreable, sens = 0.85, spec = sens, rrr = RRR,
                       alpha = ALPHA, power = POWER) {
  rows <- lapply(names(scoreable), function(nm) {
    fs <- scoreable[[nm]]
    d <- design(prevalence, p_hyper, p_hypo, sens, spec, fs, rrr, alpha, power)
    data.frame(scenario = nm, frac_scoreable = round(fs, 4), n_total_enrolled = round(d$n_total),
               patients_to_screen = round(d$patients_to_screen),
               screen_multiplier_vs_full_coverage = if (fs != 0) round(1 / fs, 2) else NA_real_)
  })
  do.call(rbind, rows)
}
