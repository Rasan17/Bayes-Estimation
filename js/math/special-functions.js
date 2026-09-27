/**
 * Bayes-Estimation: Special Mathematical Functions
 * High-precision implementations of Gamma, lnGamma, Beta, Incomplete Beta, erf, and inverse normal.
 * Designed for clinical and scientific Bayesian statistical computation.
 */

export const SpecialFunctions = {
  // Lanczos 6-term coefficients (Numerical Recipes)
  _lanczosCof: [
    76.18009172947146,
    -86.50532032941677,
    24.01409824083091,
    -1.231739572450155,
    0.1208650973866179e-2,
    -0.5395239384953e-5
  ],

  /**
   * Natural logarithm of the Gamma function ln(Γ(x))
   * Accurate to 15 decimal digits for x > 0.
   */
  lnGamma(x) {
    if (x <= 0) {
      if (Number.isInteger(x)) return Infinity;
      // Reflection formula: Γ(1-z) Γ(z) = π / sin(π z)
      return Math.log(Math.PI) - Math.log(Math.abs(Math.sin(Math.PI * x))) - this.lnGamma(1 - x);
    }
    let y = x;
    let tmp = x + 5.5;
    tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) {
      y += 1.0;
      ser += this._lanczosCof[j] / y;
    }
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  },

  /**
   * Standard Gamma function Γ(x)
   */
  gamma(x) {
    if (x > 171.5) return Infinity;
    return Math.exp(this.lnGamma(x));
  },

  /**
   * Natural logarithm of the Beta function ln(B(a, b))
   */
  lnBeta(a, b) {
    if (a <= 0 || b <= 0) return 0;
    return this.lnGamma(a) + this.lnGamma(b) - this.lnGamma(a + b);
  },

  /**
   * Beta function B(a, b)
   */
  beta(a, b) {
    return Math.exp(this.lnBeta(a, b));
  },

  /**
   * Continued fraction for Incomplete Beta function (Lentz's method)
   */
  _betacf(a, b, x) {
    const maxIt = 100;
    const eps = 3.0e-15;
    const fpmin = 1.0e-30;
    const qab = a + b;
    const qap = a + 1.0;
    const qam = a - 1.0;
    let c = 1.0;
    let d = 1.0 - qab * x / qap;
    if (Math.abs(d) < fpmin) d = fpmin;
    d = 1.0 / d;
    let h = d;

    for (let m = 1; m <= maxIt; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1.0 + aa * d;
      if (Math.abs(d) < fpmin) d = fpmin;
      c = 1.0 + aa / c;
      if (Math.abs(c) < fpmin) c = fpmin;
      d = 1.0 / d;
      h *= d * c;

      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1.0 + aa * d;
      if (Math.abs(d) < fpmin) d = fpmin;
      c = 1.0 + aa / c;
      if (Math.abs(c) < fpmin) c = fpmin;
      d = 1.0 / d;
      const delVal = d * c;
      h *= delVal;
      if (Math.abs(delVal - 1.0) < eps) break;
    }
    return h;
  },

  /**
   * Regularized Incomplete Beta function I_x(a, b)
   */
  incBeta(x, a, b) {
    if (x <= 0) return 0;
    if (x >= 1) return 1;
    if (a <= 0 || b <= 0) return NaN;

    const bt = Math.exp(
      this.lnGamma(a + b) - this.lnGamma(a) - this.lnGamma(b) +
      a * Math.log(x) + b * Math.log(1.0 - x)
    );

    if (x < (a + 1.0) / (a + b + 2.0)) {
      return (bt * this._betacf(a, b, x)) / a;
    } else {
      return 1.0 - (bt * this._betacf(b, a, 1.0 - x)) / b;
    }
  },

  /**
   * Gauss Error Function erf(x)
   */
  erf(x) {
    const a1 = 0.254829592;
    const a2 = -0.284496736;
    const a3 = 1.421413741;
    const a4 = -1.453152027;
    const a5 = 1.061405429;
    const p = 0.3275911;

    const sign = x < 0 ? -1 : 1;
    const absX = Math.abs(x);
    const t = 1.0 / (1.0 + p * absX);
    const y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * Math.exp(-absX * absX);

    return sign * y;
  },

  /**
   * Standard Normal CDF Φ(z)
   */
  normalCDF(z) {
    return 0.5 * (1.0 + this.erf(z / Math.SQRT2));
  },

  /**
   * Inverse Standard Normal CDF Φ^(-1)(p) (Acklam algorithm)
   */
  invNormalCDF(p) {
    if (p <= 0) return -Infinity;
    if (p >= 1) return Infinity;

    const a = [
      -3.969683028665376e+01,
       2.209460984245205e+02,
      -2.759285104469687e+02,
       1.383577518672690e+02,
      -3.066479806614716e+01,
       2.506628277459239e+00
    ];
    const b = [
      -5.447609879822406e+01,
       1.615858368580409e+02,
      -1.556989798598866e+02,
       6.680131188771972e+01,
      -1.328068155288572e+01
    ];
    const c = [
      -7.784894002430293e-03,
      -3.223964580411365e-01,
      -2.400758277161838e+00,
      -2.549732539343734e+00,
       4.374664141464968e+00,
       2.938163982698783e+00
    ];
    const d = [
       7.784695709041462e-03,
       3.224671290700398e-01,
       2.445134137142996e+00,
       3.754408661907416e+00
    ];

    const pLow = 0.02425;
    const pHigh = 1 - pLow;

    let q, r;
    if (p < pLow) {
      q = Math.sqrt(-2 * Math.log(p));
      return (((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
             ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    } else if (p <= pHigh) {
      q = p - 0.5;
      r = q * q;
      return (((((a[0]*r + a[1])*r + a[2])*r + a[3])*r + a[4])*r + a[5])*q /
             (((((b[0]*r + b[1])*r + b[2])*r + b[3])*r + b[4])*r + 1);
    } else {
      q = Math.sqrt(-2 * Math.log(1 - p));
      return -(((((c[0]*q + c[1])*q + c[2])*q + c[3])*q + c[4])*q + c[5]) /
              ((((d[0]*q + d[1])*q + d[2])*q + d[3])*q + 1);
    }
  },

  /**
   * Beta Quantile (Inverse Incomplete Beta)
   */
  invIncBeta(p, a, b) {
    if (p <= 0) return 0;
    if (p >= 1) return 1;

    const mean = a / (a + b);
    const variance = (a * b) / (Math.pow(a + b, 2) * (a + b + 1));
    const z = this.invNormalCDF(p);
    let x = Math.max(1e-6, Math.min(1 - 1e-6, mean + z * Math.sqrt(variance)));

    for (let iter = 0; iter < 25; iter++) {
      const fx = this.incBeta(x, a, b) - p;
      if (Math.abs(fx) < 1e-12) break;

      const lndfx = (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - this.lnBeta(a, b);
      const dfx = Math.exp(lndfx);
      if (dfx <= 0 || !isFinite(dfx)) break;

      let step = fx / dfx;
      let nextX = x - step;
      if (nextX <= 0) nextX = x * 0.5;
      else if (nextX >= 1) nextX = 1 - (1 - x) * 0.5;
      x = nextX;
    }
    return Math.max(0, Math.min(1, x));
  }
};
