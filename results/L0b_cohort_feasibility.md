# L0b — does the §4.3 cohort exist?

*Descriptive feasibility count on Challenge 2019. No models. SPEC.md §1.5, §4.3, §6.2.*

## Verdict

**PRIMARY THRESHOLD VIABLE (5687 patients)**

## FiO2 scale and validity

Observed min -50.000, median 0.500, max 4000.000 over 129,365 recorded values. 100.0% are <= 1.0, so the column is recorded as a **fraction**. A mixed scale would make the ratio meaningless; this one is not mixed.

But the raw range is impossible: room air is 0.209 and pure oxygen is 1.0. **395 values (0.31%) fall outside [0.21, 1.0] and are rejected as S/F candidates** (DECISIONS.md D10). They include zeros, which send S/F to infinity, and values above 1.0, which *deflate* S/F and would falsely admit a patient as severely hypoxemic. The patient is not excluded — only the impossible hour is.

## Attrition

|   step | label                                           |     n |   excluded | reason                         |
|-------:|:------------------------------------------------|------:|-----------:|:-------------------------------|
|      0 | set A patients                                  | 40336 |        nan |                                |
|      1 | age >= 18                                       | 40314 |         22 | paediatric or missing age      |
|      2 | ICULOS >= 24h                                   | 30646 |       9668 | short record                   |
|      3 | has >=1 hour with FiO2 AND O2Sat                | 14471 |      16175 | S/F not computable             |
|      4 | after O2Sat<97 and hour<=168 at candidate hours |  7451 |       7020 | oximeter saturated or too late |
|      5 | S/F <= 235 (primary)                            |  5687 |       1764 |                                |


At the relaxed threshold (S/F <= 315, a surrogate for P/F <= 300): **7010** patients.

## Minimum S/F per patient, among those where it is computable

| statistic | value |
|---|---|
| count | 7451.0 |
| mean | 183.3 |
| std | 80.1 |
| min | 24.0 |
| 5% | 89.0 |
| 25% | 120.0 |
| 50% | 186.0 |
| 75% | 235.0 |
| 95% | 320.0 |
| max | 459.5 |


## Who the cohort selects (vs. eligible but excluded)

| Variable | In cohort | Excluded | Std. diff |
|---|---|---|---|
| age | 64.698 | 61.592 | 0.197 |
| n_hours | 52.261 | 43.379 | 0.329 |
| sepsis_any | 0.154 | 0.039 | 0.397 |


Standardised differences above 0.1 indicate imbalance. This table belongs in Table 1 and in the README: the cohort is a selected subpopulation, not a random sample of ICU patients, and the selection is induced by the measurement itself (§6.2).
