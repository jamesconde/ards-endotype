### What this tool answers

A **prognostically enriched** trial enrols only patients a classifier flags as high-risk, because they have more
events and so need fewer patients for the same power. Two things erode that saving, and this tool puts numbers on
both:

1. **Classifier accuracy.** Every misclassified low-risk patient dilutes the enrolled event rate toward the
   low-risk rate, so the required sample size climbs back toward an unstratified trial.
2. **Whether the site can score patients at all.** A patient missing the classifier's inputs cannot be scored and
   so cannot be enrolled. In the study cohort (hypoxemic respiratory failure) the best-performing model needs a blood-gas bicarbonate that one
   hospital system records for 2.1% of patients, so the same trial needs about 52 times more screening there.

### What it is not

**Prognostic, not predictive, enrichment.** It asks which patients have more events, which is directly observable.
It does not ask which patients benefit more from treatment: that needs a treatment-effect estimate, and in
observational data it is confounded by indication. The relative risk reduction is a design assumption you enter.

### Method

- Sample size: two-sided two-proportion test without continuity correction (Fleiss), identical to
  `stats::power.prop.test`.
- Enrolled event rate = PPV × high-risk rate + (1 − PPV) × low-risk rate.
- Patients to screen = enrolled ÷ (share screening positive × share scoreable).
- **Check by simulation** runs whole trials patient by patient and never uses these formulas.

### Caveats on the default numbers

The defaults come from a hypoxemic respiratory failure cohort (oximetry surrogate, not ARDS-adjudicated) with
**sepsis onset, not mortality,** as the endpoint. The method is the deliverable; the numbers are illustrative.
