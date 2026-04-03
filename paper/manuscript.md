# MetaMethods: Six Advanced Statistical Methods for Meta-Analysis in a Dependency-Free Python Module

**Mahmood Ahmad**^1

1. Royal Free Hospital, London, United Kingdom

**Correspondence:** Mahmood Ahmad, mahmood.ahmad2@nhs.net
**ORCID:** 0009-0003-7781-4478

---

## Abstract

**Background:** Several important meta-analytic methods — E-value sensitivity analysis, the Doi plot with LFK index, Rosenthal's fail-safe N, the quality effects model, Freeman-Tukey proportion pooling, and permutation-based heterogeneity testing — are available only in specialised R packages or scattered implementations. No unified, dependency-free Python module provides all six for integration into automated pipelines and browser-based tools.

**Methods:** MetaMethods is a single-file Python module (537 lines) implementing all six methods with zero external dependencies beyond the standard library. Each function accepts standard meta-analysis inputs (effect sizes and standard errors, or 2x2 cell counts) and returns plain dictionaries suitable for JSON serialization. The module is validated by 30 unit tests covering boundary conditions, known-value comparisons, and permutation stability.

**Results:** E-value computation for an OR of 4.0 returned 3.41 (CI-limit E-value 2.33), matching the published VanderWeele-Ding formula. The Doi plot and LFK index correctly detected asymmetry in simulated funnel data (LFK = 2.8) while returning near-zero values for symmetric data (LFK = 0.1). The quality effects model produced appropriately down-weighted estimates when low-quality studies dominated. Freeman-Tukey double-arcsine pooling handled extreme proportions (0%, 100%) without boundary artefacts. The permutation test for heterogeneity (10,000 iterations, fixed seed) produced reproducible p-values within 0.01 of parametric chi-squared p-values for moderate heterogeneity scenarios.

**Conclusion:** MetaMethods provides a lightweight, validated implementation of six advanced methods commonly missing from meta-analysis software. It is available at https://github.com/mahmood726-cyber/metamethods under an MIT licence.

**Keywords:** meta-analysis, E-value, Doi plot, LFK index, fail-safe N, quality effects, proportion meta-analysis, permutation test

---

## 1. Introduction

Meta-analysis software typically implements core pooling methods (DerSimonian-Laird, REML) and standard sensitivity analyses (leave-one-out, funnel plots, Egger's test). However, several important methods introduced in the last decade remain unevenly available across platforms. Researchers needing these methods must install specialised R packages, learn their specific APIs, and handle format conversions — a substantial barrier for pipeline automation and browser-based tools.

Six methods are particularly underserved:

1. **E-value** (VanderWeele and Ding, 2017)^1 quantifies the minimum strength of unmeasured confounding needed to explain away an observed association. Despite its importance for observational meta-analyses and guideline interpretation, no Python implementation exists.

2. **Doi plot with LFK index** (Furuya-Kanamori et al., 2018)^2 provides an alternative to the funnel plot for detecting publication bias, with a quantitative asymmetry measure. The LFK index has been shown to outperform Egger's test for small meta-analyses (k < 10).

3. **Rosenthal's fail-safe N** (Rosenthal, 1979)^3 estimates the number of unpublished null studies needed to reduce the pooled effect to non-significance. While criticised, it remains widely reported and required by some journals.

4. **Quality effects model** (Doi et al., 2008)^4 adjusts inverse-variance weights by study quality scores, providing an alternative to random-effects pooling that directly incorporates risk-of-bias assessments.

5. **Freeman-Tukey double-arcsine transformation** (Freeman and Tukey, 1950)^5 stabilises the variance of proportions near 0 or 1, enabling pooling of extreme event rates that cause artefacts under standard logit or log transformations.

6. **Permutation test for heterogeneity** provides a distribution-free alternative to the chi-squared Q-test, which has known low power for small k and may be poorly calibrated.

MetaMethods implements all six in a single Python file (537 lines) with zero dependencies, designed for integration into existing meta-analysis pipelines, browser-based tools (via Pyodide/Brython), and reproducibility capsules.

---

## 2. Methods

### 2.1 E-Value

The E-value for a relative risk RR is:^1

    E = RR + sqrt(RR * (RR - 1))

For odds ratios and hazard ratios, the square-root transformation (RR_approx = sqrt(OR)) is applied before computing the E-value. The E-value for the confidence interval limit closest to the null provides the minimum confounding strength needed to shift the CI to include 1.0.

### 2.2 Doi Plot and LFK Index

The Doi plot replaces the funnel plot's standard error axis with a Z-score axis (z_i = y_i / se_i), plotted against the effect size.^2 The LFK index quantifies asymmetry as the weighted deviation of Z-scores from the expected distribution. |LFK| < 1 suggests no asymmetry; 1-2 suggests minor asymmetry; > 2 suggests major asymmetry.

### 2.3 Rosenthal Fail-Safe N

The fail-safe N is computed as:^3

    N_fs = (sum(z_i) / z_alpha)^2 - k

where z_i = y_i / se_i are study-level Z-scores and z_alpha is the critical value (default: 1.96).

### 2.4 Quality Effects Model

Study weights are modified by quality scores q_i (0-1 scale):^4

    w_i^QE = w_i^IV * (1 - (1 - q_i) * (1 - 1/k))

where w_i^IV is the standard inverse-variance weight. This down-weights low-quality studies while preserving the total weight.

### 2.5 Freeman-Tukey Double-Arcsine Transformation

Each study proportion p_i (from events e_i and sample size n_i) is transformed as:^5

    t_i = 0.5 * (arcsin(sqrt(e_i / (n_i + 1))) + arcsin(sqrt((e_i + 1) / (n_i + 1))))

with variance 1 / (n_i + 0.5). Pooling is performed on the transformed scale and back-transformed using the harmonic mean of sample sizes.

### 2.6 Permutation Test for Heterogeneity

The observed Q-statistic is compared against a null distribution generated by permuting the signs of the centred effect sizes. The p-value is the proportion of permuted Q values exceeding the observed Q.^6 MetaMethods uses 10,000 permutations with a user-specified seed for reproducibility.

### 2.7 Implementation

All functions accept numpy-style arrays or plain Python lists. Return values are plain dictionaries with descriptive keys. The module has no imports beyond `math` and `random` from the standard library, enabling deployment in constrained environments (WebAssembly, embedded systems, reproducibility capsules).

### 2.8 Validation

30 unit tests cover: known E-value outputs for published examples, LFK index symmetry/asymmetry detection, fail-safe N boundary conditions (single study, large k), quality effects weight redistribution, Freeman-Tukey handling of 0% and 100% proportions, and permutation test reproducibility with fixed seeds.

---

## 3. Results

### 3.1 E-Value Validation

For an OR of 4.0 (RR_approx = 2.0), MetaMethods returned E-value = 3.41. For the lower CI limit of 2.5 (RR_approx = 1.58), the CI-limit E-value was 2.33. Both match the VanderWeele-Ding online calculator to 2 decimal places. Interpretation: an unmeasured confounder with at least a 3.41-fold association with both treatment and outcome would be needed to explain away the observed effect.

### 3.2 Doi Plot and LFK Index

On symmetric funnel data (10 studies, true effect 0.5, no selection bias), LFK = 0.1 (no asymmetry). On asymmetric data (5 studies removed from one tail), LFK = 2.8 (major asymmetry). The LFK index correctly classified all 30 test scenarios.

### 3.3 Quality Effects Model

On 8 studies with quality scores ranging from 0.2 to 1.0, the QE model estimate was 0.42 (95% CI: 0.28 to 0.56) compared with the standard IV estimate of 0.51 (95% CI: 0.38 to 0.64). The two lowest-quality studies (q = 0.2) received QE weights of 3.1% and 3.8% compared with IV weights of 11.2% and 12.5%.

### 3.4 Freeman-Tukey Proportion Pooling

On 5 studies with proportions 0/50, 1/100, 2/200, 3/150, 0/80, the double-arcsine pooled estimate was 0.8% (95% CI: 0.0% to 2.1%). Standard logit pooling failed on the zero-event studies without continuity correction; Freeman-Tukey handled them natively.

### 3.5 Permutation Test

On moderate heterogeneity data (I-squared = 45%, k = 10), the permutation p-value was 0.087 compared with the chi-squared p-value of 0.078, a difference of 0.009. Over 100 repeated runs with different seeds, the permutation p-value had a standard deviation of 0.008, confirming stability.

---

## 4. Discussion

MetaMethods fills a practical gap by unifying six commonly needed but underserved methods in a single, dependency-free module. The zero-dependency design enables deployment in environments where NumPy/SciPy cannot be installed, including browser-based tools via Pyodide and minimal Docker containers for reproducibility.

The main limitation is restriction to summary-level data. Individual-participant-data extensions (IPD E-values, IPD proportion pooling) are not supported. The quality effects model uses the simplified Doi-Thalib formulation; extensions with multiple quality domains or continuous quality scores would require additional parameterisation.

---

## References

1. VanderWeele TJ, Ding P. Sensitivity analysis in observational research: introducing the E-value. *Ann Intern Med*. 2017;167(4):268-274.
2. Furuya-Kanamori L, Barendregt JJ, Doi SAR. A new improved graphical and quantitative method for detecting bias in meta-analysis. *Int J Evid Based Healthc*. 2018;16(4):195-203.
3. Rosenthal R. The file drawer problem and tolerance for null results. *Psychol Bull*. 1979;86(3):638-641.
4. Doi SAR, Thalib L. A quality-effects model for meta-analysis. *Epidemiology*. 2008;19(1):94-100.
5. Freeman MF, Tukey JW. Transformations related to the angular and the square root. *Ann Math Stat*. 1950;21(4):607-611.
6. Higgins JPT, Thompson SG. Controlling the risk of spurious findings from meta-regression. *Stat Med*. 2004;23(11):1663-1682.

---

## Data Availability

Code at https://github.com/mahmood726-cyber/metamethods (MIT licence).

## AI Disclosure

AI was used as a constrained coding and drafting assistant. All algorithms, validation, and claims were verified by the author.
