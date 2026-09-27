/**
 * Bayes-Estimation: Probability Distributions Module
 * Analytical PDFs, CDFs, Quantiles, and Moments for Bayesian Prior/Posterior Modeling.
 */

import { SpecialFunctions } from "./special-functions.js";

export const Distributions = {
  /**
   * Beta Distribution PDF: f(x; a, b)
   */
  betaPDF(x, a, b) {
    if (x < 0 || x > 1) return 0;
    if (x === 0) return a < 1 ? Infinity : (a === 1 ? b : 0);
    if (x === 1) return b < 1 ? Infinity : (b === 1 ? a : 0);
    if (a <= 0 || b <= 0) return 0;

    const logP = (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - SpecialFunctions.lnBeta(a, b);
    return Math.exp(logP);
  },

  /**
   * Beta Distribution CDF: F(x; a, b)
   */
  betaCDF(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    return SpecialFunctions.incBeta(x, a, b);
  },

  /**
   * Beta Distribution Quantile (Inverse CDF)
   */
  betaQuantile(p, a, b) {
    return SpecialFunctions.invIncBeta(p, a, b);
  },

  /**
   * Normal Distribution PDF: f(x; μ, σ)
   */
  normalPDF(x, mean = 0, sd = 1) {
    if (sd <= 0) return 0;
    const z = (x - mean) / sd;
    return Math.exp(-0.5 * z * z) / (sd * Math.sqrt(2 * Math.PI));
  },

  /**
   * Normal Distribution CDF: F(x; μ, σ)
   */
  normalCDF(x, mean = 0, sd = 1) {
    if (sd <= 0) return x >= mean ? 1 : 0;
    return SpecialFunctions.normalCDF((x - mean) / sd);
  },

  /**
   * Normal Distribution Quantile (Inverse CDF)
   */
  normalQuantile(p, mean = 0, sd = 1) {
    if (sd <= 0) return mean;
    return mean + sd * SpecialFunctions.invNormalCDF(p);
  },

  /**
   * Student-t Distribution PDF: f(t; ν, loc, scale)
   */
  studentTPDF(x, df, loc = 0, scale = 1) {
    if (scale <= 0 || df <= 0) return 0;
    const t = (x - loc) / scale;
    const logNorm = SpecialFunctions.lnGamma((df + 1) / 2) -
                    (0.5 * Math.log(df * Math.PI) + SpecialFunctions.lnGamma(df / 2) + Math.log(scale));
    const logVal = logNorm - ((df + 1) / 2) * Math.log(1 + (t * t) / df);
    return Math.exp(logVal);
  },

  /**
   * Student-t Distribution CDF
   */
  studentTCDF(x, df, loc = 0, scale = 1) {
    if (scale <= 0 || df <= 0) return 0;
    const t = (x - loc) / scale;
    const xBeta = df / (df + t * t);
    const pBeta = 0.5 * SpecialFunctions.incBeta(xBeta, df / 2, 0.5);
    return t >= 0 ? 1 - pBeta : pBeta;
  },

  /**
   * Student-t Quantile (Inverse CDF) via bisection & Newton
   */
  studentTQuantile(p, df, loc = 0, scale = 1) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;
    if (p === 0.5) return loc;

    // Use normal quantile as initial guess
    const z = SpecialFunctions.invNormalCDF(p);
    let t = z * Math.sqrt((df - 0.67) / df); // approximate correction

    // Newton refinement
    for (let i = 0; i < 25; i++) {
      const cdf = this.studentTCDF(t, df, 0, 1);
      const err = cdf - p;
      if (Math.abs(err) < 1e-11) break;
      const pdf = this.studentTPDF(t, df, 0, 1);
      if (pdf <= 0) break;
      t -= err / pdf;
    }
    return loc + scale * t;
  },

  /**
   * Cauchy Distribution PDF (used for JZS Bayes factors)
   */
  cauchyPDF(x, loc = 0, scale = 1) {
    if (scale <= 0) return 0;
    const z = (x - loc) / scale;
    return 1 / (Math.PI * scale * (1 + z * z));
  },

  /**
   * Cauchy CDF
   */
  cauchyCDF(x, loc = 0, scale = 1) {
    const z = (x - loc) / scale;
    return 0.5 + Math.atan(z) / Math.PI;
  },

  /**
   * Inverse Gamma PDF: f(x; α, β) = β^α / Γ(α) * x^(-α - 1) * exp(-β/x)
   */
  invGammaPDF(x, shape, scale) {
    if (x <= 0 || shape <= 0 || scale <= 0) return 0;
    const logP = shape * Math.log(scale) - SpecialFunctions.lnGamma(shape) - (shape + 1) * Math.log(x) - (scale / x);
    return Math.exp(logP);
  }
};
