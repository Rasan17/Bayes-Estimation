# Bayes-Estimation — Bayesian A/B Testing & Parameter Estimation Suite

**Bayes-Estimation** is an ultra-premium, interactive, zero-dependency clinical biostatistics and scientific Bayesian estimation suite. It is designed for biomedical researchers, clinicians, data scientists, and neurosurgeons to perform robust two-group comparisons, Bayesian A/B testing, credible interval estimation, Highest Density Interval (HDI) computation, Region of Practical Equivalence (ROPE) decision-making, and clinical diagnostic inference.

Conceived, supervised design and testing: **Dr G Narenthiran** MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN).  
Copyright © Dr G Narenthiran, `g_narenthiran@hotmail.com`.  
Dedicated to **Mrs Nirmaladevy Ganesalingam BSc (mother)**.

> ⚠️ **Notice**: *Bayes-Estimation* is an interactive Bayesian analytical suite intended for clinical research, data exploration, prospective studies, and teaching. Calculations run entirely client-side with zero external dependencies.

---

## 🚀 Key Modules & Capabilities

### 1. Bayesian Two-Group Comparison: Proportions & A/B Testing
- **Conjugate Beta-Binomial Updating**: Exact conjugate posterior updating for Group A (Control / Standard of Care) and Group B (Treatment / Experimental) with customizable $\text{Beta}(\alpha, \beta)$ priors.
- **Joint Posterior Simulation (40,000+ draws)**:
  - **Absolute Risk Difference (ARR)**: $\delta = p_B - p_A$.
  - **Relative Risk (RR)**: $p_B / p_A$ with median and 95% HDI.
  - **Odds Ratio (OR)**: $[p_B / (1 - p_B)] / [p_A / (1 - p_A)]$.
  - **Number Needed to Treat (NNT)**: $1 / |\delta|$.
- **Probability of Superiority**: Direct calculation of $P(p_B > p_A \mid \text{data})$.
- **Highest Density Interval (HDI)**: Gold-standard shortest interval containing 95% of posterior probability mass.
- **ROPE (Region of Practical Equivalence)**: Clinical equivalence bounds $[-r, +r]$ with automated Bayesian decision logic (*Superiority Established*, *Practical Equivalence Established*, *Inferiority Established*, or *Inconclusive*).
- **Savage-Dickey Density Ratio**: Bayes Factor ($BF_{10}$ and $BF_{01}$) quantifying evidence for difference vs null.

### 2. Bayesian Two-Group Comparison: Continuous Means (BEST t-Test)
- **Unequal Variances Bayesian Estimation**: Robust estimation allowing distinct variances between groups.
- **Dual Input Modes**: Summary statistics ($n, \bar{x}, s$) or raw numerical vector input (comma/space/newline separated).
- **Difference in Means ($\Delta\mu = \mu_B - \mu_A$)**: Posterior distribution with 95% HDI and ROPE analysis.
- **Bayesian Cohen's $d$**: Standardized effect size distribution $(\mu_B - \mu_A) / \sqrt{(\sigma_A^2 + \sigma_B^2)/2}$.
- **Variance Ratio ($\sigma_B / \sigma_A$)**: Quantifies heteroscedasticity.
- **JZS Bayes Factor**: Rouder et al. (2009) Cauchy-prior model comparison against point null.

### 3. Single Proportion Estimation (Beta-Binomial)
- **Bayesian Tri-Plot**: Live overlaid visualization of Prior (amber dashed), Normalized Likelihood (green dotted), and Posterior (blue solid).
- **Moments & Intervals**: Mean, median, mode, variance, 95% Equal-Tailed CrI, and 95% HDI.
- **Benchmark Hypothesis Testing**: Posterior probability that rate exceeds clinical threshold $p_0$ ($P(p > p_0 \mid \text{data})$).

### 4. Single Mean Estimation (Normal-Normal)
- **Normal-Normal Conjugate Model**: Precision-weighted conjugate updating combining prior beliefs with experimental data.
- **Continuous Tri-Plot**: Visual representation of Prior, Data Likelihood, and Posterior distributions.
- **Benchmark Exceedance**: Probability of meeting or exceeding benchmark target $\mu_0$.

### 5. Clinical Diagnostic Bayes & Fagan Nomogram
- **Likelihood Ratios**: Positive ($LR^+$) and Negative ($LR^-$) Likelihood Ratios from Sensitivity and Specificity.
- **Bayes' Rule in Odds Form**: $\text{Post-test Odds} = \text{Pre-test Odds} \times LR$.
- **Interactive SVG Fagan Nomogram**: Visual alignment ray connecting Pre-test Probability $\rightarrow$ Likelihood Ratio $\rightarrow$ Post-test Probability.
- **Uncertainty Propagation**: Beta-distributed sampling of test sensitivity and specificity providing 95% HDIs for post-test probabilities.

---

## 💻 Technical Highlights

- **Zero Dependencies**: Pure native HTML5, modern CSS3 (glassmorphic dark/light design system), and modular ES6 JavaScript.
- **High-Precision Special Functions**: Lanczos 6-term Gamma / lnGamma ($10^{-16}$ precision), Lentz continued fraction Incomplete Beta ($I_x(a, b)$), Acklam rational inverse normal CDF, and Abramowitz-Stegun error function.
- **High-DPI Canvas & SVG Graphics**: Responsive, crisp, retina-ready charts with shaded HDI and ROPE bands, mode/mean indicators, and 1-click high-res PNG export.
- **Medical Journal / APA Reporting**: 1-click clinical report generator formatted for direct inclusion in manuscripts, protocols, and study reports.
- **macOS Launcher**: Convenient double-clickable `Launch_Bayes_Estimation.command` script.

---

## 🏃 Quick Start

### Running Locally
Double-click `Launch_Bayes_Estimation.command` or simply open `index.html` in Safari, Chrome, Firefox, or Edge.

```bash
# Open in default browser
open index.html

# Run mathematical verification test suite
node verify_math.js
```

---

## 📄 License & Attribution
- Software Conceived and Designed by **Dr G Narenthiran** MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN), Southampton, UK (`g_narenthiran@hotmail.com`).
- Dedicated to **Mrs Nirmaladevy Ganesalingam BSc (mother)**.
