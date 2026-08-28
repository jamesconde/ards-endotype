# L0b — does the §4.3 cohort exist?

*Descriptive feasibility count on Challenge 2019 set A. No models. SPEC.md §1.5, §4.3, §6.2.*

## Verdict

**PRIMARY THRESHOLD VIABLE (4682 patients)**

## FiO2 scale and validity

Observed min 0.000, median 0.500, max 10.000 over 112,155 recorded values. 100.0% are <= 1.0, so the column is recorded as a **fraction**. A mixed scale would make the ratio meaningless; this one is not mixed.

But the raw range is impossible: room air is 0.209 and pure oxygen is 1.0. **321 values (0.29%) fall outside [0.21, 1.0] and are rejected as S/F candidates** (DECISIONS.md D10). They include zeros, which send S/F to infinity, and values above 1.0, which *deflate* S/F and would falsely admit a patient as severely hypoxemic. The patient is not excluded — only the impossible hour is.

## Attrition

|   step | label                                           |     n |   excluded | reason                         |
|-------:|:------------------------------------------------|------:|-----------:|:-------------------------------|
|      0 | set A patients                                  | 20336 |        nan |                                |
|      1 | age >= 18                                       | 20336 |          0 | paediatric or missing age      |
|      2 | ICULOS >= 24h                                   | 15837 |       4499 | short record                   |
|      3 | has >=1 hour with FiO2 AND O2Sat                |  9739 |       6098 | S/F not computable             |
|      4 | after O2Sat<97 and hour<=168 at candidate hours |  5621 |       4118 | oximeter saturated or too late |
|      5 | S/F <= 235 (primary)                            |  4682 |        939 |                                |


At the relaxed threshold (S/F <= 315, a surrogate for P/F <= 300): **5552** patients.

## Minimum S/F per patient, among those where it is computable

| statistic | value |
|---|---|
| count | 5621.0 |
| mean | 167.7 |
| std | 59.6 |
| min | 24.0 |
| 5% | 89.0 |
| 25% | 115.0 |
| 50% | 182.0 |
| 75% | 193.0 |
| 95% | 254.3 |
| max | 457.1 |


## Who the cohort selects (vs. eligible but excluded)

| Variable | In cohort | Excluded | Std. diff |
|---|---|---|---|
| age | 65.175 | 62.156 | 0.191 |
| n_hours | 50.666 | 42.366 | 0.338 |
| sepsis_any | 0.139 | 0.044 | 0.333 |


Standardised differences above 0.1 indicate imbalance. This table belongs in Table 1 and in the README: the cohort is a selected subpopulation, not a random sample of ICU patients, and the selection is induced by the measurement itself (§6.2).
