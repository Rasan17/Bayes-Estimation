/**
 * Bayes-Estimation — Standalone Unified Clinical Bayesian Engine
 * Compatible with local file:// protocol (Chrome, Safari, Firefox, Edge) and web servers.
 * Conceived, supervised design and testing: Dr G Narenthiran MB ChB BSc(MedSci)(Hons) MRCS(Ed.) FEBNS FRCS(SN), g_narenthiran@hotmail.com
 * Dedicated to Mrs Nirmaladevy Ganesalingam BSc (mother).
 */

(function () {
  'use strict';

  // =========================================================================
  // 1. SPECIAL MATHEMATICAL FUNCTIONS
  // =========================================================================
  const SpecialFunctions = {
    _lanczosCof: [
      76.18009172947146,
      -86.50532032941677,
      24.01409824083091,
      -1.231739572450155,
      0.1208650973866179e-2,
      -0.5395239384953e-5
    ],

    lnGamma(x) {
      if (x <= 0) {
        if (Number.isInteger(x)) return Infinity;
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

    gamma(x) {
      if (x > 171.5) return Infinity;
      return Math.exp(this.lnGamma(x));
    },

    lnBeta(a, b) {
      if (a <= 0 || b <= 0) return 0;
      return this.lnGamma(a) + this.lnGamma(b) - this.lnGamma(a + b);
    },

    beta(a, b) {
      return Math.exp(this.lnBeta(a, b));
    },

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

    normalCDF(z) {
      return 0.5 * (1.0 + this.erf(z / Math.SQRT2));
    },

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

  // =========================================================================
  // 2. PROBABILITY DISTRIBUTIONS
  // =========================================================================
  const Distributions = {
    betaPDF(x, a, b) {
      if (x < 0 || x > 1) return 0;
      if (x === 0) return a < 1 ? Infinity : (a === 1 ? b : 0);
      if (x === 1) return b < 1 ? Infinity : (b === 1 ? a : 0);
      if (a <= 0 || b <= 0) return 0;
      const logP = (a - 1) * Math.log(x) + (b - 1) * Math.log(1 - x) - SpecialFunctions.lnBeta(a, b);
      return Math.exp(logP);
    },

    betaCDF(x, a, b) {
      if (x <= 0) return 0;
      if (x >= 1) return 1;
      return SpecialFunctions.incBeta(x, a, b);
    },

    betaQuantile(p, a, b) {
      return SpecialFunctions.invIncBeta(p, a, b);
    },

    normalPDF(x, mean = 0, sd = 1) {
      if (sd <= 0) return 0;
      const z = (x - mean) / sd;
      return Math.exp(-0.5 * z * z) / (sd * Math.sqrt(2 * Math.PI));
    },

    normalCDF(x, mean = 0, sd = 1) {
      if (sd <= 0) return x >= mean ? 1 : 0;
      return SpecialFunctions.normalCDF((x - mean) / sd);
    },

    normalQuantile(p, mean = 0, sd = 1) {
      if (sd <= 0) return mean;
      return mean + sd * SpecialFunctions.invNormalCDF(p);
    },

    studentTPDF(x, df, loc = 0, scale = 1) {
      if (scale <= 0 || df <= 0) return 0;
      const t = (x - loc) / scale;
      const logNorm = SpecialFunctions.lnGamma((df + 1) / 2) -
                      (0.5 * Math.log(df * Math.PI) + SpecialFunctions.lnGamma(df / 2) + Math.log(scale));
      const logVal = logNorm - ((df + 1) / 2) * Math.log(1 + (t * t) / df);
      return Math.exp(logVal);
    },

    studentTCDF(x, df, loc = 0, scale = 1) {
      if (scale <= 0 || df <= 0) return 0;
      const t = (x - loc) / scale;
      const xBeta = df / (df + t * t);
      const pBeta = 0.5 * SpecialFunctions.incBeta(xBeta, df / 2, 0.5);
      return t >= 0 ? 1 - pBeta : pBeta;
    },

    studentTQuantile(p, df, loc = 0, scale = 1) {
      if (p <= 0) return -Infinity;
      if (p >= 1) return Infinity;
      if (p === 0.5) return loc;

      const z = SpecialFunctions.invNormalCDF(p);
      let t = z * Math.sqrt((df - 0.67) / df);

      for (let i = 0; i < 25; i++) {
        const cdf = this.studentTCDF(t, df, 0, 1);
        const err = cdf - p;
        if (Math.abs(err) < 1e-11) break;
        const pdf = this.studentTPDF(t, df, 0, 1);
        if (pdf <= 0) break;
        t -= err / pdf;
      }
      return loc + scale * t;
    }
  };

  // =========================================================================
  // 3. MCMC SAMPLER & HDI ENGINE
  // =========================================================================
  const MCMCSampler = {
    randomNormal(mean = 0, sd = 1) {
      let u1 = 0, u2 = 0;
      while (u1 === 0) u1 = Math.random();
      while (u2 === 0) u2 = Math.random();
      const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
      return mean + z0 * sd;
    },

    randomGamma(shape, scale = 1) {
      if (shape < 1) {
        const u = Math.random();
        return this.randomGamma(shape + 1, scale) * Math.pow(u, 1 / shape);
      }
      const d = shape - 1 / 3;
      const c = 1 / Math.sqrt(9 * d);

      while (true) {
        let x, v;
        do {
          x = this.randomNormal(0, 1);
          v = 1 + c * x;
        } while (v <= 0);

        v = v * v * v;
        const u = Math.random();

        if (u < 1 - 0.0331 * x * x * x * x) return d * v * scale;
        if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) return d * v * scale;
      }
    },

    randomBeta(a, b) {
      if (a <= 0 || b <= 0) return 0.5;
      const ga = this.randomGamma(a, 1);
      const gb = this.randomGamma(b, 1);
      const sum = ga + gb;
      if (sum === 0) return 0.5;
      return ga / sum;
    },

    summarize(samples) {
      const n = samples.length;
      if (n === 0) return null;

      let sum = 0;
      for (let i = 0; i < n; i++) sum += samples[i];
      const mean = sum / n;

      let sqDiffSum = 0;
      for (let i = 0; i < n; i++) {
        const diff = samples[i] - mean;
        sqDiffSum += diff * diff;
      }
      const sd = Math.sqrt(sqDiffSum / (n - 1));

      const sorted = new Float64Array(samples).sort();
      const mid = Math.floor(n / 2);
      const median = n % 2 !== 0 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;

      const kde = this.estimateKDE(sorted, 100);
      let maxDensity = -1;
      let mode = mean;
      for (const pt of kde) {
        if (pt.density > maxDensity) {
          maxDensity = pt.density;
          mode = pt.x;
        }
      }

      return { mean, median, mode, sd, sorted };
    },

    computeHDI(sortedSamples, credibleMass = 0.95) {
      const n = sortedSamples.length;
      const span = Math.floor(credibleMass * n);
      if (span >= n) {
        return { low: sortedSamples[0], high: sortedSamples[n - 1], width: sortedSamples[n - 1] - sortedSamples[0] };
      }

      let minWidth = Infinity;
      let bestLow = sortedSamples[0];
      let bestHigh = sortedSamples[span];

      const numIntervals = n - span;
      for (let i = 0; i < numIntervals; i++) {
        const low = sortedSamples[i];
        const high = sortedSamples[i + span];
        const width = high - low;
        if (width < minWidth) {
          minWidth = width;
          bestLow = low;
          bestHigh = high;
        }
      }

      return { low: bestLow, high: bestHigh, width: minWidth, mass: credibleMass };
    },

    computeCrI(sortedSamples, credibleMass = 0.95) {
      const n = sortedSamples.length;
      const alpha = 1 - credibleMass;
      const lowIdx = Math.max(0, Math.floor((alpha / 2) * n));
      const highIdx = Math.min(n - 1, Math.floor((1 - alpha / 2) * n));
      return { low: sortedSamples[lowIdx], high: sortedSamples[highIdx], mass: credibleMass };
    },

    evaluateROPE(sortedSamples, ropeLow, ropeHigh, hdi95) {
      const n = sortedSamples.length;
      let countBelow = 0, countIn = 0, countAbove = 0;

      for (let i = 0; i < n; i++) {
        const val = sortedSamples[i];
        if (val < ropeLow) countBelow++;
        else if (val > ropeHigh) countAbove++;
        else countIn++;
      }

      const pctBelow = (countBelow / n) * 100;
      const pctIn = (countIn / n) * 100;
      const pctAbove = (countAbove / n) * 100;

      let decision = "Undecided / Data Inconclusive";
      let decisionDetail = "The 95% HDI overlaps the ROPE boundary. More clinical trial samples or narrower equivalence margins required.";
      let statusClass = "status-neutral";

      if (hdi95.low > ropeHigh) {
        decision = "Superiority Established";
        decisionDetail = "The 95% HDI is completely above the ROPE region. Clear statistical evidence of clinically meaningful superiority.";
        statusClass = "status-positive";
      } else if (hdi95.high < ropeLow) {
        decision = "Inferiority Established";
        decisionDetail = "The 95% HDI is completely below the ROPE region. Clear statistical evidence of clinical inferiority.";
        statusClass = "status-negative";
      } else if (hdi95.low >= ropeLow && hdi95.high <= ropeHigh) {
        decision = "Practical Equivalence Established";
        decisionDetail = "The 95% HDI lies entirely within the ROPE region. The two interventions are practically equivalent.";
        statusClass = "status-equivalent";
      }

      return { ropeLow, ropeHigh, pctBelow, pctIn, pctAbove, decision, decisionDetail, statusClass };
    },

    estimateKDE(sortedSamples, numGridPoints = 120, xMin = null, xMax = null) {
      const n = sortedSamples.length;
      if (n === 0) return [];

      let sum = 0;
      for (let i = 0; i < n; i++) sum += sortedSamples[i];
      const mean = sum / n;

      let sqDiff = 0;
      for (let i = 0; i < n; i++) sqDiff += Math.pow(sortedSamples[i] - mean, 2);
      const sd = Math.sqrt(sqDiff / (n - 1)) || 1e-4;

      const q25 = sortedSamples[Math.floor(0.25 * n)];
      const q75 = sortedSamples[Math.floor(0.75 * n)];
      const iqr = (q75 - q25) || sd;
      const h = 0.9 * Math.min(sd, iqr / 1.34) * Math.pow(n, -0.2) || 1e-3;

      if (xMin === null) xMin = sortedSamples[0] - 3 * h;
      if (xMax === null) xMax = sortedSamples[n - 1] + 3 * h;

      const step = (xMax - xMin) / (numGridPoints - 1);
      const grid = [];
      const invH = 1 / h;
      const normFactor = 1 / (n * h * Math.sqrt(2 * Math.PI));

      for (let i = 0; i < numGridPoints; i++) {
        const x = xMin + i * step;
        let density = 0;
        for (let j = 0; j < n; j++) {
          const u = (x - sortedSamples[j]) * invH;
          density += Math.exp(-0.5 * u * u);
        }
        density *= normFactor;
        grid.push({ x, density });
      }

      return grid;
    },

    savageDickeyBF(priorDensityAtNull, posteriorSortedSamples, nullVal = 0) {
      const n = posteriorSortedSamples.length;
      if (n === 0) return { bf10: 1, bf01: 1, interpretation: "Insufficient data" };

      let sum = 0;
      for (let i = 0; i < n; i++) sum += posteriorSortedSamples[i];
      const mean = sum / n;

      let sqDiff = 0;
      for (let i = 0; i < n; i++) sqDiff += Math.pow(posteriorSortedSamples[i] - mean, 2);
      const sd = Math.sqrt(sqDiff / (n - 1)) || 1e-4;
      const h = 1.06 * sd * Math.pow(n, -0.2);

      let posteriorDensityAtNull = 0;
      const invH = 1 / h;
      for (let i = 0; i < n; i++) {
        const u = (nullVal - posteriorSortedSamples[i]) * invH;
        posteriorDensityAtNull += Math.exp(-0.5 * u * u);
      }
      posteriorDensityAtNull /= (n * h * Math.sqrt(2 * Math.PI));

      if (posteriorDensityAtNull < 1e-12) {
        return {
          bf10: 999999,
          bf01: 0.000001,
          interpretation: "Extreme evidence for H₁ (difference ≠ 0)",
          posteriorDensityAtNull: 0
        };
      }

      const bf01 = posteriorDensityAtNull / priorDensityAtNull;
      const bf10 = 1 / bf01;

      let interpretation = "Inconclusive evidence";
      if (bf10 > 100) interpretation = "Extreme evidence for H₁ (difference ≠ 0)";
      else if (bf10 > 30) interpretation = "Very strong evidence for H₁";
      else if (bf10 > 10) interpretation = "Strong evidence for H₁";
      else if (bf10 > 3) interpretation = "Moderate evidence for H₁";
      else if (bf10 > 1) interpretation = "Anecdotal / weak evidence for H₁";
      else if (bf10 === 1) interpretation = "No evidence either way";
      else if (bf10 > 1 / 3) interpretation = "Anecdotal / weak evidence for H₀ (no difference)";
      else if (bf10 > 1 / 10) interpretation = "Moderate evidence for H₀ (no difference)";
      else if (bf10 > 1 / 30) interpretation = "Strong evidence for H₀ (no difference)";
      else if (bf10 > 1 / 100) interpretation = "Very strong evidence for H₀ (no difference)";
      else interpretation = "Extreme evidence for H₀ (no difference)";

      return { bf10, bf01, priorDensityAtNull, posteriorDensityAtNull, interpretation };
    }
  };

  // =========================================================================
  // 3b. BOOTSTRAP RESAMPLING ENGINE (ESTIMATION STATISTICS - GARDNER-ALTMAN / CUMMING)
  // =========================================================================
  const BootstrapEngine = {
    B: 5000,

    resampleProportions(kA, nA, kB, nB, B = 5000) {
      const pA = kA / nA;
      const pB = kB / nB;
      const diffs = new Float64Array(B);
      const hs = new Float64Array(B);

      for (let b = 0; b < B; b++) {
        let kA_star = 0;
        for (let i = 0; i < nA; i++) if (Math.random() < pA) kA_star++;
        let kB_star = 0;
        for (let i = 0; i < nB; i++) if (Math.random() < pB) kB_star++;

        const pA_star = kA_star / nA;
        const pB_star = kB_star / nB;
        diffs[b] = pB_star - pA_star;
        hs[b] = 2 * (Math.asin(Math.sqrt(pB_star)) - Math.asin(Math.sqrt(pA_star)));
      }

      diffs.sort();
      hs.sort();
      const lowIdx = Math.floor(0.025 * B);
      const highIdx = Math.min(B - 1, Math.floor(0.975 * B));

      return {
        diffCI: { low: diffs[lowIdx], high: diffs[highIdx] },
        cohensHCI: { low: hs[lowIdx], high: hs[highIdx] },
        iterations: B
      };
    },

    resampleMeans(meanA, sdA, nA, meanB, sdB, nB, B = 5000) {
      const sampleA = new Float64Array(nA);
      const sampleB = new Float64Array(nB);

      let sumA = 0, sumB = 0;
      for (let i = 0; i < nA; i++) { sampleA[i] = MCMCSampler.randomNormal(0, 1); sumA += sampleA[i]; }
      for (let i = 0; i < nB; i++) { sampleB[i] = MCMCSampler.randomNormal(0, 1); sumB += sampleB[i]; }

      const rawMeanA = sumA / nA;
      const rawMeanB = sumB / nB;

      let varA = 0, varB = 0;
      for (let i = 0; i < nA; i++) varA += Math.pow(sampleA[i] - rawMeanA, 2);
      for (let i = 0; i < nB; i++) varB += Math.pow(sampleB[i] - rawMeanB, 2);

      const rawSdA = Math.sqrt(varA / Math.max(1, nA - 1)) || 1;
      const rawSdB = Math.sqrt(varB / Math.max(1, nB - 1)) || 1;

      for (let i = 0; i < nA; i++) sampleA[i] = ((sampleA[i] - rawMeanA) / rawSdA) * sdA + meanA;
      for (let i = 0; i < nB; i++) sampleB[i] = ((sampleB[i] - rawMeanB) / rawSdB) * sdB + meanB;

      const diffs = new Float64Array(B);
      const cohenDs = new Float64Array(B);
      const hedgesGs = new Float64Array(B);
      const dfPooled = nA + nB - 2;
      const jCorrection = 1 - (3 / (4 * Math.max(1, dfPooled) - 1));

      for (let b = 0; b < B; b++) {
        let bSumA = 0, bSumB = 0;
        for (let i = 0; i < nA; i++) bSumA += sampleA[Math.floor(Math.random() * nA)];
        for (let i = 0; i < nB; i++) bSumB += sampleB[Math.floor(Math.random() * nB)];
        const bMeanA = bSumA / nA;
        const bMeanB = bSumB / nB;

        let bVarA = 0, bVarB = 0;
        for (let i = 0; i < nA; i++) bVarA += Math.pow(sampleA[Math.floor(Math.random() * nA)] - bMeanA, 2);
        for (let i = 0; i < nB; i++) bVarB += Math.pow(sampleB[Math.floor(Math.random() * nB)] - bMeanB, 2);

        const bSdA = Math.sqrt(bVarA / Math.max(1, nA - 1));
        const bSdB = Math.sqrt(bVarB / Math.max(1, nB - 1));
        const bDiff = bMeanB - bMeanA;
        const bPooledVar = ((nA - 1) * bSdA * bSdA + (nB - 1) * bSdB * bSdB) / Math.max(1, dfPooled);
        const bPooledSd = Math.sqrt(bPooledVar) || 1e-6;
        const bD = bDiff / bPooledSd;

        diffs[b] = bDiff;
        cohenDs[b] = bD;
        hedgesGs[b] = bD * jCorrection;
      }

      diffs.sort();
      cohenDs.sort();
      hedgesGs.sort();

      const lowIdx = Math.floor(0.025 * B);
      const highIdx = Math.min(B - 1, Math.floor(0.975 * B));

      return {
        diffCI: { low: diffs[lowIdx], high: diffs[highIdx] },
        cohensDCI: { low: cohenDs[lowIdx], high: cohenDs[highIdx] },
        hedgesGCI: { low: hedgesGs[lowIdx], high: hedgesGs[highIdx] },
        iterations: B
      };
    }
  };

  // =========================================================================
  // 4. STATISTICAL MODULES (ESTIMATION & BAYESIAN ANALYSIS)
  // =========================================================================
  const TwoProportions = {
    analyze({
      kA, nA, aA = 1, bA = 1,
      kB, nB, aB = 1, bB = 1,
      ropeLow = -0.03, ropeHigh = 0.03,
      numSamples = 40000
    }) {
      kA = Math.max(0, parseInt(kA, 10) || 0);
      nA = Math.max(1, parseInt(nA, 10) || 1);
      kB = Math.max(0, parseInt(kB, 10) || 0);
      nB = Math.max(1, parseInt(nB, 10) || 1);

      if (kA > nA) kA = nA;
      if (kB > nB) kB = nB;

      aA = Math.max(0.001, parseFloat(aA) || 1);
      bA = Math.max(0.001, parseFloat(bA) || 1);
      aB = Math.max(0.001, parseFloat(aB) || 1);
      bB = Math.max(0.001, parseFloat(bB) || 1);

      // 1. Estimation Statistics: Sample Rates & Effect Sizes
      const sampleRateA = kA / nA;
      const sampleRateB = kB / nB;
      const sampleDiff = sampleRateB - sampleRateA;
      const sampleRR = sampleRateA > 0 ? sampleRateB / sampleRateA : Infinity;
      const sampleRRR = sampleRateA > 0 ? ((sampleRateB - sampleRateA) / sampleRateA) * 100 : 0;
      const sampleOddsA = sampleRateA / Math.max(1e-7, 1 - sampleRateA);
      const sampleOddsB = sampleRateB / Math.max(1e-7, 1 - sampleRateB);
      const sampleOR = sampleOddsA > 0 ? sampleOddsB / sampleOddsA : Infinity;

      // Standardized Effect Size: Cohen's h = 2 * (arcsin(sqrt(pB)) - arcsin(sqrt(pA)))
      const cohensH = 2 * (Math.asin(Math.sqrt(sampleRateB)) - Math.asin(Math.sqrt(sampleRateA)));
      const absH = Math.abs(cohensH);
      let cohensHInterpretation = "Negligible effect";
      if (absH >= 0.80) cohensHInterpretation = "Large effect size (|h| ≥ 0.8)";
      else if (absH >= 0.50) cohensHInterpretation = "Medium effect size (0.5 ≤ |h| < 0.8)";
      else if (absH >= 0.20) cohensHInterpretation = "Small effect size (0.2 ≤ |h| < 0.5)";

      // Estimation Statistics: 5,000 Resample Bootstrap Confidence Intervals
      const bootProp = BootstrapEngine.resampleProportions(kA, nA, kB, nB, 5000);
      const cohensHCI = bootProp.cohensHCI;
      const diffCI = bootProp.diffCI;

      // 2. Bayesian Conjugate Posterior Updating
      const postAA = aA + kA;
      const postBA = bA + (nA - kA);
      const postAB = aB + kB;
      const postBB = bB + (nB - kB);

      const meanA = postAA / (postAA + postBA);
      const varA = (postAA * postBA) / (Math.pow(postAA + postBA, 2) * (postAA + postBA + 1));
      const sdA = Math.sqrt(varA);
      const modeA = (postAA > 1 && postBA > 1) ? (postAA - 1) / (postAA + postBA - 2) : meanA;
      const cri95A = {
        low: Distributions.betaQuantile(0.025, postAA, postBA),
        high: Distributions.betaQuantile(0.975, postAA, postBA)
      };

      const meanB = postAB / (postAB + postBB);
      const varB = (postAB * postBB) / (Math.pow(postAB + postBB, 2) * (postAB + postBB + 1));
      const sdB = Math.sqrt(varB);
      const modeB = (postAB > 1 && postBB > 1) ? (postAB - 1) / (postAB + postBB - 2) : meanB;
      const cri95B = {
        low: Distributions.betaQuantile(0.025, postAB, postBB),
        high: Distributions.betaQuantile(0.975, postAB, postBB)
      };

      // Monte Carlo joint sampling
      const diffSamples = new Float64Array(numSamples);
      const rrSamples = new Float64Array(numSamples);
      const orSamples = new Float64Array(numSamples);

      let superiorCount = 0;

      for (let i = 0; i < numSamples; i++) {
        const pA = MCMCSampler.randomBeta(postAA, postBA);
        const pB = MCMCSampler.randomBeta(postAB, postBB);
        const diff = pB - pA;
        diffSamples[i] = diff;
        if (diff > 0) superiorCount++;

        const safePA = Math.max(1e-7, Math.min(1 - 1e-7, pA));
        const safePB = Math.max(1e-7, Math.min(1 - 1e-7, pB));
        rrSamples[i] = safePB / safePA;
        orSamples[i] = (safePB / (1 - safePB)) / (safePA / (1 - safePA));
      }

      const probSuperiority = superiorCount / numSamples;
      const probInferiority = 1 - probSuperiority;

      const diffSummary = MCMCSampler.summarize(diffSamples);
      const diffHdi95 = MCMCSampler.computeHDI(diffSummary.sorted, 0.95);
      const diffHdi90 = MCMCSampler.computeHDI(diffSummary.sorted, 0.90);
      const diffCri95 = MCMCSampler.computeCrI(diffSummary.sorted, 0.95);

      const ropeAnalysis = MCMCSampler.evaluateROPE(diffSummary.sorted, ropeLow, ropeHigh, diffHdi95);

      const rrSummary = MCMCSampler.summarize(rrSamples);
      const rrHdi95 = MCMCSampler.computeHDI(rrSummary.sorted, 0.95);
      const rrCri95 = MCMCSampler.computeCrI(rrSummary.sorted, 0.95);

      const orSummary = MCMCSampler.summarize(orSamples);
      const orHdi95 = MCMCSampler.computeHDI(orSummary.sorted, 0.95);
      const orCri95 = MCMCSampler.computeCrI(orSummary.sorted, 0.95);

      const absDiff = Math.abs(diffSummary.mean);
      const nnt = absDiff > 1e-5 ? (1 / absDiff) : Infinity;

      // Savage-Dickey Bayes Factor
      const priorSimN = 20000;
      const priorDeltaSamples = new Float64Array(priorSimN);
      for (let i = 0; i < priorSimN; i++) {
        const pA = MCMCSampler.randomBeta(aA, bA);
        const pB = MCMCSampler.randomBeta(aB, bB);
        priorDeltaSamples[i] = pB - pA;
      }
      const priorDeltaSorted = priorDeltaSamples.sort();
      const priorKde = MCMCSampler.estimateKDE(priorDeltaSorted, 50, -0.1, 0.1);
      let priorDensityAtNull = 1.0;
      for (const pt of priorKde) {
        if (Math.abs(pt.x) < 0.01) {
          priorDensityAtNull = pt.density;
          break;
        }
      }

      const bayesFactor = MCMCSampler.savageDickeyBF(priorDensityAtNull, diffSummary.sorted, 0);

      return {
        sampleStats: {
          rateA: sampleRateA,
          rateB: sampleRateB,
          diff: sampleDiff,
          diffCI,
          cohensH,
          cohensHCI,
          rr: sampleRR,
          rrr: sampleRRR,
          or: sampleOR,
          cohensH,
          cohensHInterpretation
        },
        groupA: { kA, nA, prior: { a: aA, b: bA }, posterior: { a: postAA, b: postBA }, mean: meanA, mode: modeA, sd: sdA, cri95: cri95A },
        groupB: { kB, nB, prior: { a: aB, b: bB }, posterior: { a: postAB, b: postBB }, mean: meanB, mode: modeB, sd: sdB, cri95: cri95B },
        difference: {
          mean: diffSummary.mean,
          median: diffSummary.median,
          mode: diffSummary.mode,
          sd: diffSummary.sd,
          hdi95: diffHdi95,
          hdi90: diffHdi90,
          cri95: diffCri95,
          probSuperiority,
          probInferiority,
          rope: ropeAnalysis,
          kde: MCMCSampler.estimateKDE(diffSummary.sorted, 100)
        },
        relativeRisk: { mean: rrSummary.mean, median: rrSummary.median, sd: rrSummary.sd, hdi95: rrHdi95, cri95: rrCri95 },
        oddsRatio: { mean: orSummary.mean, median: orSummary.median, sd: orSummary.sd, hdi95: orHdi95, cri95: orCri95 },
        nnt: isFinite(nnt) ? nnt : null,
        bayesFactor,
        samplesCount: numSamples
      };
    }
  };

  const TwoMeans = {
    parseData(text) {
      if (!text || typeof text !== "string") return [];
      return text
        .split(/[\s,;\t\n]+/)
        .map(v => parseFloat(v.trim()))
        .filter(v => !isNaN(v) && isFinite(v));
    },

    getSampleStats(data) {
      const n = data.length;
      if (n < 2) return null;
      let sum = 0;
      for (let i = 0; i < n; i++) sum += data[i];
      const mean = sum / n;
      let sqDiffSum = 0;
      for (let i = 0; i < n; i++) sqDiffSum += Math.pow(data[i] - mean, 2);
      const variance = sqDiffSum / (n - 1);
      const sd = Math.sqrt(variance);
      const sem = sd / Math.sqrt(n);
      return { n, mean, variance, sd, sem };
    },

    analyze(params) {
      let {
        dataA, dataB,
        nA, meanA, sdA,
        nB, meanB, sdB,
        ropeLow = -0.5, ropeHigh = 0.5,
        priorScale = 0.707,
        numSamples = 40000
      } = params;

      if (dataA && dataA.length >= 2) {
        const stats = this.getSampleStats(dataA);
        nA = stats.n; meanA = stats.mean; sdA = stats.sd;
      }
      if (dataB && dataB.length >= 2) {
        const stats = this.getSampleStats(dataB);
        nB = stats.n; meanB = stats.mean; sdB = stats.sd;
      }

      nA = Math.max(2, parseInt(nA, 10) || 2);
      nB = Math.max(2, parseInt(nB, 10) || 2);
      meanA = parseFloat(meanA) || 0;
      meanB = parseFloat(meanB) || 0;
      sdA = Math.max(1e-6, parseFloat(sdA) || 1);
      sdB = Math.max(1e-6, parseFloat(sdB) || 1);

      // 1. Estimation Statistics: Sample Effect Sizes
      const sampleDiff = meanB - meanA;
      const pooledVar = ((nA - 1) * sdA * sdA + (nB - 1) * sdB * sdB) / (nA + nB - 2);
      const pooledSd = Math.sqrt(pooledVar);
      const sampleCohensD = sampleDiff / pooledSd;
      const glassDelta = sampleDiff / sdA;

      // Common Language Effect Size (CLES) = probability that random draw B > random draw A
      const cles = Distributions.normalCDF(sampleCohensD / Math.SQRT2, 0, 1);

      const absD = Math.abs(sampleCohensD);
      let cohensDInterpretation = "Negligible effect";
      if (absD >= 1.20) cohensDInterpretation = "Very large effect (|d| ≥ 1.2)";
      else if (absD >= 0.80) cohensDInterpretation = "Large effect size (|d| ≥ 0.8)";
      else if (absD >= 0.50) cohensDInterpretation = "Medium effect size (0.5 ≤ |d| < 0.8)";
      else if (absD >= 0.20) cohensDInterpretation = "Small effect size (0.2 ≤ |d| < 0.5)";

      // Estimation Statistics: 5,000 Resample Bootstrap Confidence Intervals (Gardner-Altman / Cumming)
      const dfPooled = nA + nB - 2;
      const jCorrection = 1 - (3 / (4 * Math.max(1, dfPooled) - 1));
      const hedgesG = sampleCohensD * jCorrection;

      const bootMeans = BootstrapEngine.resampleMeans(meanA, sdA, nA, meanB, sdB, nB, 5000);
      const cohensDCI = bootMeans.cohensDCI;
      const hedgesGCI = bootMeans.hedgesGCI;
      const diffCI = bootMeans.diffCI;

      // 2. Bayesian Estimation (BEST / t-test)
      const dfA = nA - 1;
      const dfB = nB - 1;
      const semA = sdA / Math.sqrt(nA);
      const semB = sdB / Math.sqrt(nB);

      const tCritA = Distributions.studentTQuantile(0.975, dfA, 0, 1);
      const tCritB = Distributions.studentTQuantile(0.975, dfB, 0, 1);

      const cri95A = { low: meanA - tCritA * semA, high: meanA + tCritA * semA };
      const cri95B = { low: meanB - tCritB * semB, high: meanB + tCritB * semB };

      const diffSamples = new Float64Array(numSamples);
      const cohenSamples = new Float64Array(numSamples);
      const sdRatioSamples = new Float64Array(numSamples);

      let superiorCount = 0;

      for (let i = 0; i < numSamples; i++) {
        const chiA = MCMCSampler.randomGamma(dfA / 2, 2);
        const varA_sample = (dfA * sdA * sdA) / chiA;
        const sigmaA_sample = Math.sqrt(varA_sample);

        const chiB = MCMCSampler.randomGamma(dfB / 2, 2);
        const varB_sample = (dfB * sdB * sdB) / chiB;
        const sigmaB_sample = Math.sqrt(varB_sample);

        const muA_sample = MCMCSampler.randomNormal(meanA, sigmaA_sample / Math.sqrt(nA));
        const muB_sample = MCMCSampler.randomNormal(meanB, sigmaB_sample / Math.sqrt(nB));

        const diff = muB_sample - muA_sample;
        diffSamples[i] = diff;
        if (diff > 0) superiorCount++;

        const currentPooledSd = Math.sqrt((varA_sample + varB_sample) / 2);
        cohenSamples[i] = diff / currentPooledSd;
        sdRatioSamples[i] = sigmaB_sample / sigmaA_sample;
      }

      const probSuperiority = superiorCount / numSamples;
      const probInferiority = 1 - probSuperiority;

      const diffSummary = MCMCSampler.summarize(diffSamples);
      const diffHdi95 = MCMCSampler.computeHDI(diffSummary.sorted, 0.95);
      const diffHdi90 = MCMCSampler.computeHDI(diffSummary.sorted, 0.90);
      const diffCri95 = MCMCSampler.computeCrI(diffSummary.sorted, 0.95);

      const ropeAnalysis = MCMCSampler.evaluateROPE(diffSummary.sorted, ropeLow, ropeHigh, diffHdi95);

      const cohenSummary = MCMCSampler.summarize(cohenSamples);
      const cohenHdi95 = MCMCSampler.computeHDI(cohenSummary.sorted, 0.95);
      const cohenCri95 = MCMCSampler.computeCrI(cohenSummary.sorted, 0.95);

      const sdRatioSummary = MCMCSampler.summarize(sdRatioSamples);
      const sdRatioHdi95 = MCMCSampler.computeHDI(sdRatioSummary.sorted, 0.95);

      // Rouder et al. (2009) JZS Bayes Factor
      const tClassic = sampleDiff / Math.sqrt(pooledVar * (1 / nA + 1 / nB));
      const dfTotal = nA + nB - 2;
      const nEff = (nA * nB) / (nA + nB);

      let bf10 = 1.0;
      try {
        const r = priorScale;
        let integral = 0;
        const gSteps = 1000;
        const gMax = 50.0;
        const dg = gMax / gSteps;

        for (let s = 1; s < gSteps; s++) {
          const g = s * dg;
          const scaleFactor = 1 + nEff * r * r * g;
          const term1 = 1 / Math.sqrt(scaleFactor);
          const term2 = Math.pow(1 + (tClassic * tClassic) / (dfTotal * scaleFactor), -(dfTotal + 1) / 2);
          const priorG = (1 / Math.sqrt(2 * Math.PI)) * Math.pow(g, -1.5) * Math.exp(-1 / (2 * g));
          const weight = (s % 2 === 0) ? 2 : 4;
          integral += weight * term1 * term2 * priorG;
        }
        integral = (integral + (1 / Math.sqrt(2 * Math.PI)) * Math.pow(gMax, -1.5) * Math.exp(-1 / (2 * gMax))) * (dg / 3);

        const nullTerm = Math.pow(1 + (tClassic * tClassic) / dfTotal, -(dfTotal + 1) / 2);
        bf10 = Math.max(1e-6, integral / nullTerm);
      } catch (e) {
        bf10 = 1.0;
      }

      const bf01 = 1 / bf10;
      let bfInterpretation = "Inconclusive evidence";
      if (bf10 > 100) bfInterpretation = "Extreme evidence for difference (H₁)";
      else if (bf10 > 30) bfInterpretation = "Very strong evidence for difference (H₁)";
      else if (bf10 > 10) bfInterpretation = "Strong evidence for difference (H₁)";
      else if (bf10 > 3) bfInterpretation = "Moderate evidence for difference (H₁)";
      else if (bf10 > 1) bfInterpretation = "Anecdotal evidence for difference (H₁)";
      else if (bf10 === 1) bfInterpretation = "No evidence either way";
      else if (bf10 > 1 / 3) bfInterpretation = "Anecdotal evidence for equality (H₀)";
      else if (bf10 > 1 / 10) bfInterpretation = "Moderate evidence for equality (H₀)";
      else if (bf10 > 1 / 30) bfInterpretation = "Strong evidence for equality (H₀)";
      else bfInterpretation = "Extreme evidence for equality (H₀)";

      return {
        sampleStats: {
          meanA, meanB, sdA, sdB,
          diff: sampleDiff,
          diffCI,
          cohensD: sampleCohensD,
          cohensDCI,
          hedgesG,
          hedgesGCI,
          glassDelta,
          cles,
          cohensDInterpretation
        },
        groupA: { n: nA, mean: meanA, sd: sdA, sem: semA, df: dfA, cri95: cri95A },
        groupB: { n: nB, mean: meanB, sd: sdB, sem: semB, df: dfB, cri95: cri95B },
        difference: {
          mean: diffSummary.mean,
          median: diffSummary.median,
          mode: diffSummary.mode,
          sd: diffSummary.sd,
          hdi95: diffHdi95,
          hdi90: diffHdi90,
          cri95: diffCri95,
          probSuperiority,
          probInferiority,
          rope: ropeAnalysis,
          kde: MCMCSampler.estimateKDE(diffSummary.sorted, 100)
        },
        cohensD: {
          mean: cohenSummary.mean,
          median: cohenSummary.median,
          sd: cohenSummary.sd,
          hdi95: cohenHdi95,
          cri95: cohenCri95,
          kde: MCMCSampler.estimateKDE(cohenSummary.sorted, 80)
        },
        sdRatio: { mean: sdRatioSummary.mean, median: sdRatioSummary.median, hdi95: sdRatioHdi95 },
        bayesFactor: { bf10, bf01, priorScale, tClassic, interpretation: bfInterpretation }
      };
    }
  };

  const SingleProportion = {
    analyze({ k, n, priorA = 1, priorB = 1, benchmark = 0.50, numSamples = 25000 }) {
      k = Math.max(0, parseInt(k, 10) || 0);
      n = Math.max(1, parseInt(n, 10) || 1);
      if (k > n) k = n;

      priorA = Math.max(0.001, parseFloat(priorA) || 1);
      priorB = Math.max(0.001, parseFloat(priorB) || 1);
      benchmark = Math.max(0, Math.min(1, parseFloat(benchmark) || 0.50));

      const sampleRate = k / n;
      const effectDiff = sampleRate - benchmark;
      const effectH = 2 * (Math.asin(Math.sqrt(sampleRate)) - Math.asin(Math.sqrt(benchmark)));

      const postA = priorA + k;
      const postB = priorB + (n - k);

      const postMean = postA / (postA + postB);
      const postVar = (postA * postB) / (Math.pow(postA + postB, 2) * (postA + postB + 1));
      const postSd = Math.sqrt(postVar);
      const postMode = (postA > 1 && postB > 1) ? (postA - 1) / (postA + postB - 2) : postMean;

      const cri95 = {
        low: Distributions.betaQuantile(0.025, postA, postB),
        high: Distributions.betaQuantile(0.975, postA, postB)
      };

      const samples = new Float64Array(numSamples);
      let countAboveBenchmark = 0;
      for (let i = 0; i < numSamples; i++) {
        const s = MCMCSampler.randomBeta(postA, postB);
        samples[i] = s;
        if (s > benchmark) countAboveBenchmark++;
      }
      const sorted = samples.sort();
      const hdi95 = MCMCSampler.computeHDI(sorted, 0.95);
      const probExceedBenchmark = countAboveBenchmark / numSamples;

      const gridPoints = 200;
      const triPlotData = [];
      const mP = k / n;
      const maxLikVal = Math.pow(mP, k) * Math.pow(1 - mP, n - k) || 1e-10;
      const maxPostDensity = Math.max(1, Distributions.betaPDF(postMode, postA, postB));

      for (let i = 0; i <= gridPoints; i++) {
        const p = i / gridPoints;
        const priorDens = Distributions.betaPDF(p, priorA, priorB);
        const postDens = Distributions.betaPDF(p, postA, postB);
        const rawLik = (p === 0 && k === 0) || (p === 1 && k === n) ? 1 :
                      (p > 0 && p < 1) ? Math.pow(p, k) * Math.pow(1 - p, n - k) : 0;
        const normLik = (rawLik / maxLikVal) * maxPostDensity;
        triPlotData.push({ p, prior: priorDens, likelihood: normLik, posterior: postDens });
      }

      return {
        k, n, sampleRate, effectDiff, effectH,
        prior: { a: priorA, b: priorB },
        posterior: { a: postA, b: postB, mean: postMean, mode: postMode, sd: postSd, cri95, hdi95 },
        benchmark: { value: benchmark, probExceed: probExceedBenchmark },
        triPlotData
      };
    }
  };

  const SingleMean = {
    analyze({ sampleMean, sampleSd, sampleN, priorMean = 0, priorSd = 10, benchmark = 0 }) {
      sampleMean = parseFloat(sampleMean) || 0;
      sampleSd = Math.max(1e-6, parseFloat(sampleSd) || 1);
      sampleN = Math.max(1, parseInt(sampleN, 10) || 1);
      priorMean = parseFloat(priorMean) || 0;
      priorSd = Math.max(1e-6, parseFloat(priorSd) || 10);
      benchmark = parseFloat(benchmark) || 0;

      const effectDiff = sampleMean - benchmark;
      const effectD = effectDiff / sampleSd;

      const sampleSem = sampleSd / Math.sqrt(sampleN);
      const tauPrior = 1 / (priorSd * priorSd);
      const tauData = 1 / (sampleSem * sampleSem);
      const tauPost = tauPrior + tauData;

      const postMean = (tauPrior * priorMean + tauData * sampleMean) / tauPost;
      const postSd = 1 / Math.sqrt(tauPost);

      const zCrit = 1.959963984540054;
      const cri95 = { low: postMean - zCrit * postSd, high: postMean + zCrit * postSd };
      const probExceed = 1 - Distributions.normalCDF(benchmark, postMean, postSd);

      const allMeans = [priorMean, sampleMean, postMean, benchmark];
      const minCenter = Math.min(...allMeans);
      const maxCenter = Math.max(...allMeans);
      const maxSd = Math.max(priorSd, sampleSem, postSd);

      const xMin = minCenter - 3.5 * maxSd;
      const xMax = maxCenter + 3.5 * maxSd;
      const steps = 150;
      const stepSize = (xMax - xMin) / steps;

      const triPlotData = [];
      for (let i = 0; i <= steps; i++) {
        const x = xMin + i * stepSize;
        triPlotData.push({
          x,
          prior: Distributions.normalPDF(x, priorMean, priorSd),
          likelihood: Distributions.normalPDF(x, sampleMean, sampleSem),
          posterior: Distributions.normalPDF(x, postMean, postSd)
        });
      }

      return {
        sample: { mean: sampleMean, sd: sampleSd, n: sampleN, sem: sampleSem, effectDiff, effectD },
        prior: { mean: priorMean, sd: priorSd },
        posterior: { mean: postMean, sd: postSd, cri95, hdi95: cri95 },
        benchmark: { value: benchmark, probExceed },
        triPlotData
      };
    }
  };

  const DiagnosticNomogram = {
    analyze({ preTestProb = 0.20, sensitivity = 0.85, specificity = 0.90, sampleN_disease = 100, sampleN_healthy = 100 }) {
      preTestProb = Math.max(0.001, Math.min(0.999, parseFloat(preTestProb) || 0.20));
      sensitivity = Math.max(0.001, Math.min(0.999, parseFloat(sensitivity) || 0.85));
      specificity = Math.max(0.001, Math.min(0.999, parseFloat(specificity) || 0.90));

      const plr = sensitivity / Math.max(1e-7, 1 - specificity);
      const nlr = (1 - sensitivity) / Math.max(1e-7, specificity);
      const dor = plr / Math.max(1e-7, nlr);

      const preTestOdds = preTestProb / (1 - preTestProb);
      const postOddsPos = preTestOdds * plr;
      const postProbPos = postOddsPos / (1 + postOddsPos);

      const postOddsNeg = preTestOdds * nlr;
      const postProbNeg = postOddsNeg / (1 + postOddsNeg);

      const simN = 10000;
      const postPosSamples = new Float64Array(simN);
      const postNegSamples = new Float64Array(simN);

      const nDis = Math.max(10, parseInt(sampleN_disease, 10) || 100);
      const nHlt = Math.max(10, parseInt(sampleN_healthy, 10) || 100);

      const sensA = sensitivity * nDis + 1;
      const sensB = (1 - sensitivity) * nDis + 1;
      const specA = specificity * nHlt + 1;
      const specB = (1 - specificity) * nHlt + 1;

      for (let i = 0; i < simN; i++) {
        const s = MCMCSampler.randomBeta(sensA, sensB);
        const sp = MCMCSampler.randomBeta(specA, specB);
        const simPlr = s / Math.max(1e-6, 1 - sp);
        const simNlr = (1 - s) / Math.max(1e-6, sp);

        const simOddsPos = preTestOdds * simPlr;
        postPosSamples[i] = simOddsPos / (1 + simOddsPos);

        const simOddsNeg = preTestOdds * simNlr;
        postNegSamples[i] = simOddsNeg / (1 + simOddsNeg);
      }

      const posSummary = MCMCSampler.summarize(postPosSamples);
      const posHdi95 = MCMCSampler.computeHDI(posSummary.sorted, 0.95);

      const negSummary = MCMCSampler.summarize(postNegSamples);
      const negHdi95 = MCMCSampler.computeHDI(negSummary.sorted, 0.95);

      return {
        preTestProb, preTestOdds, sensitivity, specificity, plr, nlr, dor,
        postTestPositive: { prob: postProbPos, hdi95: posHdi95 },
        postTestNegative: { prob: postProbNeg, hdi95: negHdi95 }
      };
    }
  };

  // =========================================================================
  // 5. THEME & PRESETS
  // =========================================================================
  const ThemeManager = {
    init() {
      const savedTheme = localStorage.getItem("bayes_estimation_theme");
      const prefersDark = window.matchMedia && window.matchMedia("(prefers-color-scheme: dark)").matches;
      const activeTheme = savedTheme || (prefersDark ? "dark" : "light");
      this.setTheme(activeTheme);

      const toggleBtn = document.getElementById("themeToggleBtn");
      if (toggleBtn) {
        toggleBtn.addEventListener("click", () => {
          const currentTheme = document.documentElement.getAttribute("data-theme") || "dark";
          const newTheme = currentTheme === "dark" ? "light" : "dark";
          this.setTheme(newTheme);
        });
      }
    },

    setTheme(theme) {
      document.documentElement.setAttribute("data-theme", theme);
      localStorage.setItem("bayes_estimation_theme", theme);

      const toggleBtn = document.getElementById("themeToggleBtn");
      if (toggleBtn) {
        toggleBtn.innerHTML = theme === "dark"
          ? `<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><circle cx="12" cy="12" r="5"></circle><line x1="12" y1="1" x2="12" y2="3"></line><line x1="12" y1="21" x2="12" y2="23"></line><line x1="4.22" y1="4.22" x2="5.64" y2="5.64"></line><line x1="18.36" y1="18.36" x2="19.78" y2="19.78"></line><line x1="1" y1="12" x2="3" y2="12"></line><line x1="21" y1="12" x2="23" y2="12"></line><line x1="4.22" y1="19.78" x2="5.64" y2="18.36"></line><line x1="18.36" y1="5.64" x2="19.78" y2="4.22"></line></svg><span>Light Mode</span>`
          : `<svg width="16" height="16" viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2"><path d="M21 12.79A9 9 0 1 1 11.21 3 7 7 0 0 0 21 12.79z"></path></svg><span>Dark Mode</span>`;
      }
    }
  };

  const Presets = {
    twoProportions: [
      {
        id: "stroke_trial",
        title: "Acute Ischemic Stroke: Thrombectomy vs Medical Therapy",
        kA: 58, nA: 267,
        kB: 114, nB: 267,
        priorA: 1, priorB: 1,
        ropeLow: -0.03, ropeHigh: 0.03
      },
      {
        id: "craniosynostosis",
        title: "Craniosynostosis: Endoscopic vs Open Transfusion Rate",
        kA: 42, nA: 50,
        kB: 7, nB: 50,
        priorA: 1, priorB: 1,
        ropeLow: -0.05, ropeHigh: 0.05
      },
      {
        id: "sah_vasospasm",
        title: "Aneurysmal SAH: Prophylactic Nimodipine vs Control",
        kA: 68, nA: 175,
        kB: 38, nB: 179,
        priorA: 1, priorB: 1,
        ropeLow: -0.04, ropeHigh: 0.04
      }
    ],

    twoMeans: [
      {
        id: "tbi_icp",
        title: "Severe TBI: Intracranial Pressure (mmHg)",
        meanA: 21.4, sdA: 4.8, nA: 64,
        meanB: 15.2, sdB: 3.9, nB: 62,
        ropeLow: -1.0, ropeHigh: 1.0,
        priorScale: 0.707
      },
      {
        id: "inph_gait",
        title: "iNPH Tap Test: Timed-Up-and-Go (TUG) Seconds",
        meanA: 24.8, sdA: 5.6, nA: 38,
        meanB: 18.2, sdB: 4.9, nB: 38,
        ropeLow: -1.5, ropeHigh: 1.5,
        priorScale: 0.707
      },
      {
        id: "glioma_eor",
        title: "Glioma Resection: Extent of Resection (% Volume)",
        meanA: 87.5, sdA: 8.4, nA: 52,
        meanB: 95.8, sdB: 4.1, nB: 48,
        ropeLow: -2.0, ropeHigh: 2.0,
        priorScale: 0.707
      }
    ],

    singleProportion: [
      {
        id: "microdiscectomy_recurrence",
        title: "Microdiscectomy: 2-Year Recurrent Disc Herniation",
        k: 14, n: 210,
        priorA: 1, priorB: 15,
        benchmark: 0.10
      },
      {
        id: "shunt_infection",
        title: "VP Shunt Infection Rate (Antibiotic Catheters)",
        k: 5, n: 160,
        priorA: 1, priorB: 1,
        benchmark: 0.05
      }
    ],

    diagnosticNomogram: [
      {
        id: "inph_tap_test",
        title: "iNPH Diagnosis: CSF Tap Test Accuracy",
        preTestProb: 0.40,
        sensitivity: 0.65,
        specificity: 0.88
      },
      {
        id: "carotid_stenosis",
        title: "Severe Carotid Stenosis (≥70%): Duplex Ultrasound",
        preTestProb: 0.25,
        sensitivity: 0.90,
        specificity: 0.92
      }
    ]
  };

  // =========================================================================
  // 6. CANVAS CHARTS
  // =========================================================================
  const CanvasCharts = {
    setupCanvas(canvas) {
      const dpr = window.devicePixelRatio || 1;
      const rect = canvas.getBoundingClientRect();
      const width = rect.width || canvas.width || 600;
      const height = rect.height || canvas.height || 300;

      canvas.width = width * dpr;
      canvas.height = height * dpr;

      const ctx = canvas.getContext("2d");
      ctx.scale(dpr, dpr);
      return { ctx, width, height, dpr };
    },

    formatNum(n, d = 3) {
      if (typeof n !== "number" || isNaN(n)) return "--";
      return Number(n.toFixed(d)).toString();
    },

    drawDifferenceDistribution(canvas, kde, hdi95, rope = null, pointEstimates = {}, options = {}) {
      if (!canvas || !kde || kde.length === 0) return;
      const { ctx, width, height } = this.setupCanvas(canvas);

      const isDark = document.documentElement.getAttribute("data-theme") === "dark";
      const textColor = isDark ? "#94a3b8" : "#475569";
      const textBright = isDark ? "#f8fafc" : "#0f172a";
      const gridColor = isDark ? "rgba(148, 163, 184, 0.12)" : "rgba(100, 116, 139, 0.15)";
      const lineColor = options.lineColor || (isDark ? "#38bdf8" : "#0284c7");
      const fillColor = options.fillColor || (isDark ? "rgba(56, 189, 248, 0.15)" : "rgba(2, 132, 199, 0.12)");
      const hdiFillColor = options.hdiFillColor || (isDark ? "rgba(56, 189, 248, 0.35)" : "rgba(2, 132, 199, 0.28)");
      const ropeColor = options.ropeColor || (isDark ? "rgba(52, 211, 153, 0.22)" : "rgba(16, 185, 129, 0.18)");
      const ropeBorder = options.ropeBorder || (isDark ? "#34d399" : "#059669");

      const padding = { top: 35, right: 30, bottom: 50, left: 55 };
      const plotW = width - padding.left - padding.right;
      const plotH = height - padding.top - padding.bottom;

      let xMin = kde[0].x;
      let xMax = kde[kde.length - 1].x;
      if (hdi95) {
        xMin = Math.min(xMin, hdi95.low);
        xMax = Math.max(xMax, hdi95.high);
      }
      if (rope) {
        xMin = Math.min(xMin, rope.ropeLow);
        xMax = Math.max(xMax, rope.ropeHigh);
      }
      xMin = Math.min(xMin, 0);
      xMax = Math.max(xMax, 0);

      const xPadding = (xMax - xMin) * 0.08 || 0.05;
      xMin -= xPadding;
      xMax += xPadding;

      let yMax = 0;
      for (const pt of kde) {
        if (pt.density > yMax) yMax = pt.density;
      }
      yMax = (yMax || 1) * 1.15;

      const getX = x => padding.left + ((x - xMin) / (xMax - xMin)) * plotW;
      const getY = y => padding.top + plotH - (y / yMax) * plotH;

      ctx.clearRect(0, 0, width, height);

      // ROPE band
      if (rope && rope.ropeLow !== undefined && rope.ropeHigh !== undefined) {
        const rx1 = Math.max(padding.left, getX(rope.ropeLow));
        const rx2 = Math.min(padding.left + plotW, getX(rope.ropeHigh));
        const rw = rx2 - rx1;

        if (rw > 0) {
          ctx.fillStyle = ropeColor;
          ctx.fillRect(rx1, padding.top, rw, plotH);

          ctx.beginPath();
          ctx.strokeStyle = ropeBorder;
          ctx.lineWidth = 1.5;
          ctx.setLineDash([4, 3]);
          ctx.moveTo(rx1, padding.top); ctx.lineTo(rx1, padding.top + plotH);
          ctx.moveTo(rx2, padding.top); ctx.lineTo(rx2, padding.top + plotH);
          ctx.stroke();
          ctx.setLineDash([]);

          ctx.fillStyle = ropeBorder;
          ctx.font = "600 11px system-ui, -apple-system, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText(`ROPE [${this.formatNum(rope.ropeLow, 2)}, ${this.formatNum(rope.ropeHigh, 2)}]`, (rx1 + rx2) / 2, padding.top + 16);
        }
      }

      // Null reference line
      const zeroX = getX(0);
      if (zeroX >= padding.left && zeroX <= padding.left + plotW) {
        ctx.beginPath();
        ctx.strokeStyle = isDark ? "#ef4444" : "#dc2626";
        ctx.lineWidth = 1.5;
        ctx.setLineDash([3, 3]);
        ctx.moveTo(zeroX, padding.top);
        ctx.lineTo(zeroX, padding.top + plotH);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = isDark ? "#f87171" : "#dc2626";
        ctx.font = "bold 10px monospace";
        ctx.textAlign = "center";
        ctx.fillText("Null (0)", zeroX, padding.top + plotH + 34);
      }

      // HDI region
      if (hdi95) {
        const hdiPts = kde.filter(pt => pt.x >= hdi95.low && pt.x <= hdi95.high);
        if (hdiPts.length > 1) {
          ctx.beginPath();
          ctx.moveTo(getX(hdiPts[0].x), getY(0));
          for (const pt of hdiPts) {
            ctx.lineTo(getX(pt.x), getY(pt.density));
          }
          ctx.lineTo(getX(hdiPts[hdiPts.length - 1].x), getY(0));
          ctx.closePath();
          ctx.fillStyle = hdiFillColor;
          ctx.fill();

          const hdiY = getY(0) - 2;
          const hx1 = getX(hdi95.low);
          const hx2 = getX(hdi95.high);

          ctx.beginPath();
          ctx.strokeStyle = textBright;
          ctx.lineWidth = 3;
          ctx.moveTo(hx1, hdiY); ctx.lineTo(hx2, hdiY);
          ctx.moveTo(hx1, hdiY - 6); ctx.lineTo(hx1, hdiY + 6);
          ctx.moveTo(hx2, hdiY - 6); ctx.lineTo(hx2, hdiY + 6);
          ctx.stroke();

          ctx.fillStyle = textBright;
          ctx.font = "bold 11px system-ui, -apple-system, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText(`95% HDI: [${this.formatNum(hdi95.low, 3)}, ${this.formatNum(hdi95.high, 3)}]`, (hx1 + hx2) / 2, hdiY - 12);
        }
      }

      // Fill full curve
      ctx.beginPath();
      ctx.moveTo(getX(kde[0].x), getY(0));
      for (const pt of kde) {
        ctx.lineTo(getX(pt.x), getY(pt.density));
      }
      ctx.lineTo(getX(kde[kde.length - 1].x), getY(0));
      ctx.closePath();
      ctx.fillStyle = fillColor;
      ctx.fill();

      // Stroke curve
      ctx.beginPath();
      ctx.moveTo(getX(kde[0].x), getY(kde[0].density));
      for (let i = 1; i < kde.length; i++) {
        ctx.lineTo(getX(kde[i].x), getY(kde[i].density));
      }
      ctx.strokeStyle = lineColor;
      ctx.lineWidth = 2.5;
      ctx.stroke();

      // Mean indicator
      if (pointEstimates.mean !== undefined) {
        const mx = getX(pointEstimates.mean);
        if (mx >= padding.left && mx <= padding.left + plotW) {
          ctx.beginPath();
          ctx.strokeStyle = lineColor;
          ctx.lineWidth = 1.5;
          ctx.setLineDash([4, 2]);
          ctx.moveTo(mx, padding.top + 20);
          ctx.lineTo(mx, padding.top + plotH);
          ctx.stroke();
          ctx.setLineDash([]);

          ctx.fillStyle = lineColor;
          ctx.font = "600 11px system-ui, -apple-system, sans-serif";
          ctx.textAlign = "center";
          ctx.fillText(`Mean: ${this.formatNum(pointEstimates.mean, 3)}`, mx, padding.top + 16);
        }
      }

      // Axis ticks
      ctx.beginPath();
      ctx.strokeStyle = gridColor;
      ctx.lineWidth = 1;
      ctx.moveTo(padding.left, padding.top + plotH);
      ctx.lineTo(padding.left + plotW, padding.top + plotH);
      ctx.stroke();

      const xTicks = 7;
      ctx.fillStyle = textColor;
      ctx.font = "11px system-ui, -apple-system, sans-serif";
      ctx.textAlign = "center";

      for (let i = 0; i < xTicks; i++) {
        const val = xMin + (i / (xTicks - 1)) * (xMax - xMin);
        const px = getX(val);
        ctx.beginPath();
        ctx.moveTo(px, padding.top + plotH);
        ctx.lineTo(px, padding.top + plotH + 5);
        ctx.stroke();
        ctx.fillText(this.formatNum(val, 2), px, padding.top + plotH + 18);
      }

      ctx.fillStyle = textBright;
      ctx.font = "600 12px system-ui, -apple-system, sans-serif";
      ctx.fillText(options.xLabel || "Parameter Difference (Group B − Group A)", padding.left + plotW / 2, padding.top + plotH + 46);
    },

    drawTriPlot(canvas, triPlotData) {
      if (!canvas || !triPlotData || triPlotData.length === 0) return;
      const { ctx, width, height } = this.setupCanvas(canvas);

      const isDark = document.documentElement.getAttribute("data-theme") === "dark";
      const textColor = isDark ? "#94a3b8" : "#475569";
      const textBright = isDark ? "#f8fafc" : "#0f172a";
      const gridColor = isDark ? "rgba(148, 163, 184, 0.12)" : "rgba(100, 116, 139, 0.15)";

      const priorColor = isDark ? "#f59e0b" : "#d97706";
      const likColor = isDark ? "#10b981" : "#059669";
      const postColor = isDark ? "#38bdf8" : "#0284c7";

      const padding = { top: 35, right: 30, bottom: 45, left: 55 };
      const plotW = width - padding.left - padding.right;
      const plotH = height - padding.top - padding.bottom;

      const xProp = triPlotData[0].p !== undefined ? "p" : "x";
      const xMin = triPlotData[0][xProp];
      const xMax = triPlotData[triPlotData.length - 1][xProp];

      let yMax = 0;
      for (const pt of triPlotData) {
        if (pt.prior > yMax) yMax = pt.prior;
        if (pt.likelihood > yMax) yMax = pt.likelihood;
        if (pt.posterior > yMax) yMax = pt.posterior;
      }
      yMax *= 1.15;

      const getX = x => padding.left + ((x - xMin) / (xMax - xMin)) * plotW;
      const getY = y => padding.top + plotH - (y / yMax) * plotH;

      ctx.clearRect(0, 0, width, height);

      // Posterior Fill
      ctx.beginPath();
      ctx.moveTo(getX(triPlotData[0][xProp]), getY(0));
      for (const pt of triPlotData) {
        ctx.lineTo(getX(pt[xProp]), getY(pt.posterior));
      }
      ctx.lineTo(getX(triPlotData[triPlotData.length - 1][xProp]), getY(0));
      ctx.closePath();
      ctx.fillStyle = isDark ? "rgba(56, 189, 248, 0.20)" : "rgba(2, 132, 199, 0.15)";
      ctx.fill();

      // Prior line
      ctx.beginPath();
      ctx.moveTo(getX(triPlotData[0][xProp]), getY(triPlotData[0].prior));
      for (let i = 1; i < triPlotData.length; i++) {
        ctx.lineTo(getX(triPlotData[i][xProp]), getY(triPlotData[i].prior));
      }
      ctx.strokeStyle = priorColor;
      ctx.lineWidth = 2;
      ctx.setLineDash([5, 4]);
      ctx.stroke();
      ctx.setLineDash([]);

      // Likelihood line
      ctx.beginPath();
      ctx.moveTo(getX(triPlotData[0][xProp]), getY(triPlotData[0].likelihood));
      for (let i = 1; i < triPlotData.length; i++) {
        ctx.lineTo(getX(triPlotData[i][xProp]), getY(triPlotData[i].likelihood));
      }
      ctx.strokeStyle = likColor;
      ctx.lineWidth = 2;
      ctx.setLineDash([2, 3]);
      ctx.stroke();
      ctx.setLineDash([]);

      // Posterior line
      ctx.beginPath();
      ctx.moveTo(getX(triPlotData[0][xProp]), getY(triPlotData[0].posterior));
      for (let i = 1; i < triPlotData.length; i++) {
        ctx.lineTo(getX(triPlotData[i][xProp]), getY(triPlotData[i].posterior));
      }
      ctx.strokeStyle = postColor;
      ctx.lineWidth = 2.5;
      ctx.stroke();

      // Axis
      ctx.beginPath();
      ctx.strokeStyle = gridColor;
      ctx.moveTo(padding.left, padding.top + plotH);
      ctx.lineTo(padding.left + plotW, padding.top + plotH);
      ctx.stroke();

      const xTicks = 6;
      ctx.fillStyle = textColor;
      ctx.font = "11px system-ui, -apple-system, sans-serif";
      ctx.textAlign = "center";
      for (let i = 0; i < xTicks; i++) {
        const val = xMin + (i / (xTicks - 1)) * (xMax - xMin);
        ctx.fillText(this.formatNum(val, 2), getX(val), padding.top + plotH + 18);
      }

      // Legend
      const legends = [
        { name: "Prior", color: priorColor, dash: [5, 4] },
        { name: "Likelihood (Scaled)", color: likColor, dash: [2, 3] },
        { name: "Posterior", color: postColor, dash: [] }
      ];
      let legX = padding.left + 10;
      ctx.font = "600 12px system-ui, -apple-system, sans-serif";
      for (const leg of legends) {
        ctx.beginPath();
        ctx.strokeStyle = leg.color;
        ctx.lineWidth = 2.5;
        ctx.setLineDash(leg.dash);
        ctx.moveTo(legX, padding.top - 16);
        ctx.lineTo(legX + 22, padding.top - 16);
        ctx.stroke();
        ctx.setLineDash([]);

        ctx.fillStyle = textBright;
        ctx.textAlign = "left";
        ctx.fillText(leg.name, legX + 28, padding.top - 12);
        legX += ctx.measureText(leg.name).width + 48;
      }
    },

    renderFaganNomogram(container, preTestProb, plr, nlr, postProbPos, postProbNeg) {
      if (!container) return;
      const isDark = document.documentElement.getAttribute("data-theme") === "dark";
      const strokeColor = isDark ? "#334155" : "#cbd5e1";
      const textColor = isDark ? "#94a3b8" : "#475569";
      const textBright = isDark ? "#f8fafc" : "#0f172a";
      const linePosColor = isDark ? "#38bdf8" : "#0284c7";
      const lineNegColor = isDark ? "#f43f5e" : "#e11d48";

      const width = 640;
      const height = 500;
      const margin = { top: 68, bottom: 44, left: 75, right: 88 };
      const axisH = height - margin.top - margin.bottom;

      const minProb = 0.001, maxProb = 0.999;
      const minOdds = Math.log(minProb / (1 - minProb));
      const maxOdds = Math.log(maxProb / (1 - maxProb));

      const probToY = p => {
        p = Math.max(minProb, Math.min(maxProb, p));
        const odds = Math.log(p / (1 - p));
        return margin.top + (1 - (odds - minOdds) / (maxOdds - minOdds)) * axisH;
      };

      const minLR = -3, maxLR = 3;
      const lrToY = lr => {
        const logLR = Math.log10(Math.max(0.001, Math.min(1000, lr)));
        return margin.top + (1 - (logLR - minLR) / (maxLR - minLR)) * axisH;
      };

      const xPre = margin.left;
      const xLR = width / 2;
      const xPost = width - margin.right;

      const yPre = probToY(preTestProb);
      const yPosLR = lrToY(plr);
      const yPostPos = probToY(postProbPos);
      const yNegLR = lrToY(nlr);
      const yPostNeg = probToY(postProbNeg);

      const probTicks = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99, 0.999];
      const lrTicks = [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001];

      let svg = `<svg viewBox="0 0 ${width} ${height}" class="fagan-nomogram-svg" style="width: 100%; height: auto; font-family: system-ui, -apple-system, sans-serif;">
        <text x="${width/2}" y="25" text-anchor="middle" font-size="14.5" font-weight="800" fill="${textBright}" letter-spacing="0.02em">Bayesian Clinical Fagan Nomogram</text>
        <line x1="${xPre}" y1="${margin.top}" x2="${xPre}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />
        <line x1="${xLR}" y1="${margin.top}" x2="${xLR}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />
        <line x1="${xPost}" y1="${margin.top}" x2="${xPost}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />

        <text x="${xPre}" y="${margin.top - 14}" text-anchor="middle" font-size="12" font-weight="700" fill="${textBright}">Pre-Test Prob</text>
        <text x="${xLR}" y="${margin.top - 14}" text-anchor="middle" font-size="12" font-weight="700" fill="${textBright}">Likelihood Ratio</text>
        <text x="${xPost}" y="${margin.top - 14}" text-anchor="middle" font-size="12" font-weight="700" fill="${textBright}">Post-Test Prob</text>
      `;

      for (const p of probTicks) {
        const y = probToY(p);
        const label = p >= 0.01 && p <= 0.99 ? `${(p * 100).toFixed(p < 0.1 ? 1 : 0)}%` : `${p * 100}%`;
        svg += `
          <line x1="${xPre - 4}" y1="${y}" x2="${xPre + 4}" y2="${y}" stroke="${strokeColor}" stroke-width="1.5" />
          <text x="${xPre - 8}" y="${y + 3}" text-anchor="end" font-size="10" fill="${textColor}">${label}</text>
          <line x1="${xPost - 4}" y1="${y}" x2="${xPost + 4}" y2="${y}" stroke="${strokeColor}" stroke-width="1.5" />
          <text x="${xPost + 8}" y="${y + 3}" text-anchor="start" font-size="10" fill="${textColor}">${label}</text>
        `;
      }

      for (const lr of lrTicks) {
        const y = lrToY(lr);
        svg += `
          <line x1="${xLR - 5}" y1="${y}" x2="${xLR + 5}" y2="${y}" stroke="${strokeColor}" stroke-width="1.5" />
          <text x="${xLR + 8}" y="${y + 3}" text-anchor="start" font-size="10" fill="${textColor}">${lr}</text>
        `;
      }

      svg += `
        <line x1="${xPre}" y1="${yPre}" x2="${xPost}" y2="${yPostPos}" stroke="${linePosColor}" stroke-width="3" stroke-linecap="round" opacity="0.9" />
        <circle cx="${xPre}" cy="${yPre}" r="5" fill="${linePosColor}" />
        <circle cx="${xLR}" cy="${yPosLR}" r="5" fill="${linePosColor}" />
        <circle cx="${xPost}" cy="${yPostPos}" r="6" fill="${linePosColor}" />
        <text x="${xPost + 14}" y="${yPostPos + 4}" font-size="11" font-weight="700" fill="${linePosColor}">+Test: ${(postProbPos * 100).toFixed(1)}%</text>

        <line x1="${xPre}" y1="${yPre}" x2="${xPost}" y2="${yPostNeg}" stroke="${lineNegColor}" stroke-width="2.5" stroke-dasharray="6,4" stroke-linecap="round" opacity="0.9" />
        <circle cx="${xLR}" cy="${yNegLR}" r="4" fill="${lineNegColor}" />
        <circle cx="${xPost}" cy="${yPostNeg}" r="5" fill="${lineNegColor}" />
        <text x="${xPost + 14}" y="${yPostNeg + 4}" font-size="11" font-weight="700" fill="${lineNegColor}">−Test: ${(postProbNeg * 100).toFixed(1)}%</text>
      </svg>`;

      container.innerHTML = svg;
    },

    exportPNG(canvas, filename = "bayes_chart.png") {
      if (!canvas) return;
      const link = document.createElement("a");
      link.download = filename;
      link.href = canvas.toDataURL("image/png");
      link.click();
    }
  };

  // =========================================================================
  // 7. MAIN APP CONTROLLER
  // =========================================================================
  const App = {
    activeTab: "twoProportions",
    lastResults: {},

    init() {
      ThemeManager.init();
      this.setupTabs();
      this.setupPresetSelectors();
      this.setupTwoProportionsModule();
      this.setupTwoMeansModule();
      this.setupSingleProportionModule();
      this.setupSingleMeanModule();
      this.setupDiagnosticModule();
      this.setupExportActions();

      this.runTwoProportions();
      this.runTwoMeans();
      this.runSingleProportion();
      this.runSingleMean();
      this.runDiagnosticNomogram();

      window.addEventListener("resize", () => {
        this.redrawActiveCharts();
      });

      const themeBtn = document.getElementById("themeToggleBtn");
      if (themeBtn) {
        themeBtn.addEventListener("click", () => {
          setTimeout(() => this.redrawActiveCharts(), 50);
        });
      }

      // Disclaimer modal interactions (Escape key or backdrop click)
      const modal = document.getElementById("disclaimerModal");
      if (modal) {
        modal.addEventListener("click", (e) => {
          if (e.target === modal) modal.classList.add("hidden");
        });
        document.addEventListener("keydown", (e) => {
          if (e.key === "Escape" && !modal.classList.contains("hidden")) {
            modal.classList.add("hidden");
          }
        });
      }
    },

    setupTabs() {
      const tabButtons = document.querySelectorAll(".tab-btn");
      tabButtons.forEach(btn => {
        btn.addEventListener("click", () => {
          const tabId = btn.getAttribute("data-tab");
          if (!tabId) return;

          tabButtons.forEach(b => b.classList.remove("active"));
          btn.classList.add("active");

          document.querySelectorAll(".module-panel").forEach(p => p.classList.remove("active"));
          const activePanel = document.getElementById(`panel-${tabId}`);
          if (activePanel) activePanel.classList.add("active");

          this.activeTab = tabId;
          setTimeout(() => this.redrawActiveCharts(), 30);
        });
      });
    },

    setupPresetSelectors() {
      const propPresetSelect = document.getElementById("presetTwoProp");
      if (propPresetSelect) {
        Presets.twoProportions.forEach(p => {
          const opt = document.createElement("option");
          opt.value = p.id;
          opt.textContent = p.title;
          propPresetSelect.appendChild(opt);
        });
        propPresetSelect.addEventListener("change", e => {
          const preset = Presets.twoProportions.find(p => p.id === e.target.value);
          if (preset) {
            document.getElementById("propKA").value = preset.kA;
            document.getElementById("propNA").value = preset.nA;
            document.getElementById("propKB").value = preset.kB;
            document.getElementById("propNB").value = preset.nB;
            document.getElementById("propRopeLow").value = preset.ropeLow;
            document.getElementById("propRopeHigh").value = preset.ropeHigh;
            this.runTwoProportions();
          }
        });
      }

      const meansPresetSelect = document.getElementById("presetTwoMeans");
      if (meansPresetSelect) {
        Presets.twoMeans.forEach(p => {
          const opt = document.createElement("option");
          opt.value = p.id;
          opt.textContent = p.title;
          meansPresetSelect.appendChild(opt);
        });
        meansPresetSelect.addEventListener("change", e => {
          const preset = Presets.twoMeans.find(p => p.id === e.target.value);
          if (preset) {
            document.getElementById("meansMeanA").value = preset.meanA;
            document.getElementById("meansSdA").value = preset.sdA;
            document.getElementById("meansNA").value = preset.nA;
            document.getElementById("meansMeanB").value = preset.meanB;
            document.getElementById("meansSdB").value = preset.sdB;
            document.getElementById("meansNB").value = preset.nB;
            document.getElementById("meansRopeLow").value = preset.ropeLow;
            document.getElementById("meansRopeHigh").value = preset.ropeHigh;
            this.runTwoMeans();
          }
        });
      }

      const singlePropPresetSelect = document.getElementById("presetSingleProp");
      if (singlePropPresetSelect) {
        Presets.singleProportion.forEach(p => {
          const opt = document.createElement("option");
          opt.value = p.id;
          opt.textContent = p.title;
          singlePropPresetSelect.appendChild(opt);
        });
        singlePropPresetSelect.addEventListener("change", e => {
          const preset = Presets.singleProportion.find(p => p.id === e.target.value);
          if (preset) {
            document.getElementById("spK").value = preset.k;
            document.getElementById("spN").value = preset.n;
            document.getElementById("spPriorA").value = preset.priorA;
            document.getElementById("spPriorB").value = preset.priorB;
            document.getElementById("spBenchmark").value = preset.benchmark;
            this.runSingleProportion();
          }
        });
      }

      const diagPresetSelect = document.getElementById("presetDiagnostic");
      if (diagPresetSelect) {
        Presets.diagnosticNomogram.forEach(p => {
          const opt = document.createElement("option");
          opt.value = p.id;
          opt.textContent = p.title;
          diagPresetSelect.appendChild(opt);
        });
        diagPresetSelect.addEventListener("change", e => {
          const preset = Presets.diagnosticNomogram.find(p => p.id === e.target.value);
          if (preset) {
            document.getElementById("diagPreProb").value = preset.preTestProb;
            document.getElementById("diagSens").value = preset.sensitivity;
            document.getElementById("diagSpec").value = preset.specificity;
            this.runDiagnosticNomogram();
          }
        });
      }
    },

    setupTwoProportionsModule() {
      const btn = document.getElementById("btnCalcTwoProp");
      if (btn) btn.addEventListener("click", () => this.runTwoProportions());

      const autoInputs = ["propKA", "propNA", "propPriorA", "propPriorA_b", "propKB", "propNB", "propPriorB", "propPriorB_b", "propRopeLow", "propRopeHigh"];
      autoInputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener("input", () => this.runTwoProportions());
      });
    },

    runTwoProportions() {
      const kA = parseInt(document.getElementById("propKA")?.value, 10) || 58;
      const nA = parseInt(document.getElementById("propNA")?.value, 10) || 267;
      const aA = parseFloat(document.getElementById("propPriorA")?.value) || 1;
      const bA = parseFloat(document.getElementById("propPriorA_b")?.value) || 1;

      const kB = parseInt(document.getElementById("propKB")?.value, 10) || 114;
      const nB = parseInt(document.getElementById("propNB")?.value, 10) || 267;
      const aB = parseFloat(document.getElementById("propPriorB")?.value) || 1;
      const bB = parseFloat(document.getElementById("propPriorB_b")?.value) || 1;

      const ropeLow = parseFloat(document.getElementById("propRopeLow")?.value) || -0.03;
      const ropeHigh = parseFloat(document.getElementById("propRopeHigh")?.value) || 0.03;

      const res = TwoProportions.analyze({
        kA, nA, aA, bA,
        kB, nB, aB, bB,
        ropeLow, ropeHigh,
        numSamples: 40000
      });
      this.lastResults.twoProportions = res;

      const pctSup = (res.difference.probSuperiority * 100).toFixed(1);
      const diffPct = (res.difference.mean * 100).toFixed(2);
      const hdiLow = (res.difference.hdi95.low * 100).toFixed(2);
      const hdiHigh = (res.difference.hdi95.high * 100).toFixed(2);

      const elProbSup = document.getElementById("propKpiProbSup");
      if (elProbSup) elProbSup.textContent = `${pctSup}%`;

      const elDiff = document.getElementById("propKpiDiff");
      if (elDiff) {
        elDiff.textContent = `${diffPct > 0 ? "+" : ""}${diffPct}%`;
        elDiff.className = `kpi-value ${res.difference.mean >= 0 ? "text-success" : "text-danger"}`;
      }

      const elHdi = document.getElementById("propKpiHdi");
      if (elHdi) elHdi.textContent = `[${hdiLow}%, ${hdiHigh}%]`;

      // Effect Size: Cohen's h & Confidence Interval
      const elDiffCI = document.getElementById("propKpiDiffCI");
      if (elDiffCI && res.sampleStats.diffCI) {
        elDiffCI.textContent = `95% Bootstrap CI: [${(res.sampleStats.diffCI.low * 100).toFixed(1)}%, ${(res.sampleStats.diffCI.high * 100).toFixed(1)}%]`;
      }
      const elCohensH = document.getElementById("propKpiCohensH");
      if (elCohensH) elCohensH.textContent = `h = ${res.sampleStats.cohensH.toFixed(3)}`;
      const elCohensHCI = document.getElementById("propKpiCohensHCI");
      if (elCohensHCI && res.sampleStats.cohensHCI) {
        elCohensHCI.textContent = `95% Bootstrap CI: [${res.sampleStats.cohensHCI.low.toFixed(3)}, ${res.sampleStats.cohensHCI.high.toFixed(3)}]`;
      }
      const elCohensHLabel = document.getElementById("propKpiCohensHLabel");
      if (elCohensHLabel) elCohensHLabel.textContent = res.sampleStats.cohensHInterpretation;

      const elRR = document.getElementById("propKpiRR");
      if (elRR) elRR.textContent = `${res.relativeRisk.median.toFixed(2)} (95% CrI ${res.relativeRisk.cri95.low.toFixed(2)}–${res.relativeRisk.cri95.high.toFixed(2)})`;

      const elOR = document.getElementById("propKpiOR");
      if (elOR) elOR.textContent = `${res.oddsRatio.median.toFixed(2)} (95% CrI ${res.oddsRatio.cri95.low.toFixed(2)}–${res.oddsRatio.cri95.high.toFixed(2)})`;

      const elNNT = document.getElementById("propKpiNNT");
      if (elNNT) elNNT.textContent = res.nnt ? Math.round(res.nnt).toString() : "N/A";

      const elBF = document.getElementById("propKpiBF");
      if (elBF) elBF.textContent = res.bayesFactor.bf10.toFixed(2);

      const badge = document.getElementById("propRopeBadge");
      if (badge) {
        badge.textContent = res.difference.rope.decision;
        badge.className = `status-badge ${res.difference.rope.statusClass}`;
      }
      const ropeDetail = document.getElementById("propRopeDetail");
      if (ropeDetail) {
        ropeDetail.textContent = `${res.difference.rope.decisionDetail} (${res.difference.rope.pctIn.toFixed(1)}% in ROPE, ${res.difference.rope.pctAbove.toFixed(1)}% above, ${res.difference.rope.pctBelow.toFixed(1)}% below).`;
      }

      const bfEl = document.getElementById("propBayesFactor");
      if (bfEl) {
        bfEl.textContent = `BF₁₀ = ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
      }

      // Update dynamic interpretation guide text
      const elIntArr = document.getElementById("interpPropArr");
      if (elIntArr) elIntArr.textContent = `${diffPct > 0 ? "+" : ""}${diffPct}% (95% Bootstrap CI: [${(res.sampleStats.diffCI.low * 100).toFixed(1)}%, ${(res.sampleStats.diffCI.high * 100).toFixed(1)}%])`;
      const elIntBF = document.getElementById("interpPropBF");
      if (elIntBF) elIntBF.textContent = res.bayesFactor.bf10.toFixed(2);
      const elIntSup = document.getElementById("interpPropSup");
      if (elIntSup) elIntSup.textContent = `${pctSup}%`;
      const elIntHdi = document.getElementById("interpPropHdi");
      if (elIntHdi) elIntHdi.textContent = `${hdiLow}% to ${hdiHigh}%`;
      const elIntH = document.getElementById("interpPropH");
      if (elIntH) elIntH.textContent = `h = ${res.sampleStats.cohensH.toFixed(3)}, ${res.sampleStats.cohensHInterpretation}`;
      const elIntNnt = document.getElementById("interpPropNnt");
      if (elIntNnt) elIntNnt.textContent = res.nnt ? Math.round(res.nnt).toString() : "N/A";
      const elIntOR = document.getElementById("interpPropOR");
      if (elIntOR) elIntOR.textContent = res.oddsRatio.median.toFixed(2);
      const elIntRR = document.getElementById("interpPropRR");
      if (elIntRR) elIntRR.textContent = res.relativeRisk.median.toFixed(2);

      const canvasDiff = document.getElementById("chartPropDiff");
      if (canvasDiff) {
        CanvasCharts.drawDifferenceDistribution(
          canvasDiff,
          res.difference.kde,
          res.difference.hdi95,
          res.difference.rope,
          { mean: res.difference.mean, mode: res.difference.mode },
          { xLabel: "Difference in Success Rates (p_Treatment − p_Control)" }
        );
      }

      const reportBox = document.getElementById("reportTwoProp");
      if (reportBox) {
        reportBox.textContent = `CLINICAL & STATISTICAL BAYESIAN REPORT: TWO PROPORTIONS (A/B)
================================================================================
ESTIMATION STATISTICS (Sample Data & Empirical Effect Sizes):
  • Group A (Control): ${res.groupA.kA} / ${res.groupA.nA} events (${(res.sampleStats.rateA*100).toFixed(1)}%)
  • Group B (Treatment): ${res.groupB.kB} / ${res.groupB.nB} events (${(res.sampleStats.rateB*100).toFixed(1)}%)
  • Sample Risk Difference (ARR): ${(res.sampleStats.diff*100).toFixed(2)}%
  • Standardized Effect Size (Cohen's h): ${res.sampleStats.cohensH.toFixed(3)} (${res.sampleStats.cohensHInterpretation})
  • Relative Risk (RR): ${res.sampleStats.rr.toFixed(3)}
  • Relative Risk Reduction (RRR): ${res.sampleStats.rrr.toFixed(1)}%
  • Sample Odds Ratio (OR): ${res.sampleStats.or.toFixed(3)}
  • Number Needed to Treat (NNT): ${res.nnt ? res.nnt.toFixed(1) : "Undefined"}

BAYESIAN ANALYSIS (Conjugate Beta Updating & Posterior Estimation):
  • Group A Posterior Beta(${res.groupA.posterior.a}, ${res.groupA.posterior.b}): Mean = ${(res.groupA.mean*100).toFixed(1)}%, 95% CrI [${(res.groupA.cri95.low*100).toFixed(1)}%, ${(res.groupA.cri95.high*100).toFixed(1)}%]
  • Group B Posterior Beta(${res.groupB.posterior.a}, ${res.groupB.posterior.b}): Mean = ${(res.groupB.mean*100).toFixed(1)}%, 95% CrI [${(res.groupB.cri95.low*100).toFixed(1)}%, ${(res.groupB.cri95.high*100).toFixed(1)}%]
  • Posterior Mean Difference (ARR): ${(res.difference.mean*100).toFixed(2)}% (SD = ${(res.difference.sd*100).toFixed(2)}%)
  • 95% Highest Density Interval (HDI): [${hdiLow}%, ${hdiHigh}%]
  • Probability of Superiority P(p_B > p_A | data): ${pctSup}%
  • Posterior Relative Risk (RR): Median ${res.relativeRisk.median.toFixed(3)} [95% HDI: ${res.relativeRisk.hdi95.low.toFixed(3)} to ${res.relativeRisk.hdi95.high.toFixed(3)}]
  • Posterior Odds Ratio (OR): Median ${res.oddsRatio.median.toFixed(3)} [95% HDI: ${res.oddsRatio.hdi95.low.toFixed(3)} to ${res.oddsRatio.hdi95.high.toFixed(3)}]

ROPE Clinical Equivalence [${(res.difference.rope.ropeLow*100).toFixed(1)}%, ${(res.difference.rope.ropeHigh*100).toFixed(1)}%]:
  • Verdict: ${res.difference.rope.decision}
  • Mass in ROPE: ${res.difference.rope.pctIn.toFixed(1)}% | Superior: ${res.difference.rope.pctAbove.toFixed(1)}% | Inferior: ${res.difference.rope.pctBelow.toFixed(1)}%
  • Savage-Dickey Bayes Factor BF₁₀: ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
      }
    },

    setupTwoMeansModule() {
      const btn = document.getElementById("btnCalcTwoMeans");
      if (btn) btn.addEventListener("click", () => this.runTwoMeans());

      const inputs = ["meansMeanA", "meansSdA", "meansNA", "meansMeanB", "meansSdB", "meansNB", "meansRopeLow", "meansRopeHigh"];
      inputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener("input", () => this.runTwoMeans());
      });

      const rawA = document.getElementById("meansRawDataA");
      const rawB = document.getElementById("meansRawDataB");
      if (rawA) rawA.addEventListener("input", () => this.runTwoMeans());
      if (rawB) rawB.addEventListener("input", () => this.runTwoMeans());
    },

    runTwoMeans() {
      const rawTextA = document.getElementById("meansRawDataA")?.value.trim();
      const rawTextB = document.getElementById("meansRawDataB")?.value.trim();

      const dataA = rawTextA ? TwoMeans.parseData(rawTextA) : null;
      const dataB = rawTextB ? TwoMeans.parseData(rawTextB) : null;

      const meanA = parseFloat(document.getElementById("meansMeanA")?.value) || 21.4;
      const sdA = parseFloat(document.getElementById("meansSdA")?.value) || 4.8;
      const nA = parseInt(document.getElementById("meansNA")?.value, 10) || 64;

      const meanB = parseFloat(document.getElementById("meansMeanB")?.value) || 15.2;
      const sdB = parseFloat(document.getElementById("meansSdB")?.value) || 3.9;
      const nB = parseInt(document.getElementById("meansNB")?.value, 10) || 62;

      const ropeLow = parseFloat(document.getElementById("meansRopeLow")?.value) || -1.0;
      const ropeHigh = parseFloat(document.getElementById("meansRopeHigh")?.value) || 1.0;

      const res = TwoMeans.analyze({
        dataA, dataB,
        meanA, sdA, nA,
        meanB, sdB, nB,
        ropeLow, ropeHigh,
        numSamples: 40000
      });
      this.lastResults.twoMeans = res;

      const pctSup = (res.difference.probSuperiority * 100).toFixed(1);
      const elProbSup = document.getElementById("meansKpiProbSup");
      if (elProbSup) elProbSup.textContent = `${pctSup}%`;

      const elDiff = document.getElementById("meansKpiDiff");
      if (elDiff) elDiff.textContent = `${res.difference.mean > 0 ? "+" : ""}${res.difference.mean.toFixed(2)}`;

      const elHdi = document.getElementById("meansKpiHdi");
      if (elHdi) elHdi.textContent = `[${res.difference.hdi95.low.toFixed(2)}, ${res.difference.hdi95.high.toFixed(2)}]`;

      const elMeansDiffCI = document.getElementById("meansKpiDiffCI");
      if (elMeansDiffCI && res.sampleStats.diffCI) {
        elMeansDiffCI.textContent = `95% Bootstrap CI: [${res.sampleStats.diffCI.low.toFixed(2)}, ${res.sampleStats.diffCI.high.toFixed(2)}]`;
      }
      const elCohen = document.getElementById("meansKpiCohensD");
      if (elCohen) elCohen.textContent = `d = ${res.sampleStats.cohensD.toFixed(3)}`;
      const elCohenCI = document.getElementById("meansKpiCohensDCI");
      if (elCohenCI && res.sampleStats.cohensDCI) {
        elCohenCI.textContent = `95% Bootstrap CI: [${res.sampleStats.cohensDCI.low.toFixed(3)}, ${res.sampleStats.cohensDCI.high.toFixed(3)}] (Hedges g: ${res.sampleStats.hedgesG.toFixed(3)})`;
      }
      const elCohenLabel = document.getElementById("meansKpiCohensDLabel");
      if (elCohenLabel) elCohenLabel.textContent = res.sampleStats.cohensDInterpretation;

      const elCLES = document.getElementById("meansKpiCLES");
      if (elCLES) elCLES.textContent = `${(res.sampleStats.cles * 100).toFixed(1)}%`;

      const elSdRatio = document.getElementById("meansKpiSdRatio");
      if (elSdRatio) elSdRatio.textContent = `${res.sdRatio.median.toFixed(2)} (95% HDI ${res.sdRatio.hdi95.low.toFixed(2)}–${res.sdRatio.hdi95.high.toFixed(2)})`;

      const elBF = document.getElementById("meansKpiBF");
      if (elBF) elBF.textContent = `BF₁₀ = ${res.bayesFactor.bf10.toFixed(2)}`;

      const badge = document.getElementById("meansRopeBadge");
      if (badge) {
        badge.textContent = res.difference.rope.decision;
        badge.className = `status-badge ${res.difference.rope.statusClass}`;
      }
      const ropeDetail = document.getElementById("meansRopeDetail");
      if (ropeDetail) {
        ropeDetail.textContent = `${res.difference.rope.decisionDetail} (${res.difference.rope.pctIn.toFixed(1)}% in ROPE, ${res.difference.rope.pctAbove.toFixed(1)}% above, ${res.difference.rope.pctBelow.toFixed(1)}% below).`;
      }

      // Update dynamic interpretation guide
      const elIntDiff = document.getElementById("interpMeansDiff");
      if (elIntDiff) elIntDiff.textContent = `${res.difference.mean > 0 ? "+" : ""}${res.difference.mean.toFixed(2)} units (95% Bootstrap CI: [${res.sampleStats.diffCI.low.toFixed(2)}, ${res.sampleStats.diffCI.high.toFixed(2)}])`;
      const elIntD = document.getElementById("interpMeansD");
      if (elIntD) elIntD.textContent = `d = ${res.sampleStats.cohensD.toFixed(3)} [95% Bootstrap CI: ${res.sampleStats.cohensDCI.low.toFixed(3)}, ${res.sampleStats.cohensDCI.high.toFixed(3)}] (${res.sampleStats.cohensDInterpretation})`;
      const elIntCles = document.getElementById("interpMeansCles");
      if (elIntCles) elIntCles.textContent = `${(res.sampleStats.cles * 100).toFixed(1)}%`;
      const elIntBF = document.getElementById("interpMeansBF");
      if (elIntBF) elIntBF.textContent = res.bayesFactor.bf10.toFixed(2);

      const canvasDiff = document.getElementById("chartMeansDiff");
      if (canvasDiff) {
        CanvasCharts.drawDifferenceDistribution(
          canvasDiff,
          res.difference.kde,
          res.difference.hdi95,
          res.difference.rope,
          { mean: res.difference.mean, mode: res.difference.mode },
          { xLabel: "Difference in Continuous Means (μ_Treatment − μ_Control)" }
        );
      }

      const canvasCohen = document.getElementById("chartMeansCohen");
      if (canvasCohen) {
        CanvasCharts.drawDifferenceDistribution(
          canvasCohen,
          res.cohensD.kde,
          res.cohensD.hdi95,
          null,
          { mean: res.cohensD.mean, mode: res.cohensD.median },
          {
            xLabel: "Bayesian Cohen's d (Standardized Effect Size)",
            lineColor: "#8b5cf6",
            fillColor: "rgba(139, 92, 246, 0.15)",
            hdiFillColor: "rgba(139, 92, 246, 0.3)"
          }
        );
      }

      const reportBox = document.getElementById("reportTwoMeans");
      if (reportBox) {
        reportBox.textContent = `CLINICAL & STATISTICAL BAYESIAN REPORT: TWO CONTINUOUS MEANS (BEST)
================================================================================
ESTIMATION STATISTICS (Sample Data & Empirical Effect Sizes):
  • Group A (Control): N = ${res.groupA.n}, Mean = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)} (SEM = ${res.groupA.sem.toFixed(2)})
  • Group B (Treatment): N = ${res.groupB.n}, Mean = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)} (SEM = ${res.groupB.sem.toFixed(2)})
  • Sample Mean Difference: ${res.sampleStats.diff.toFixed(3)}
  • Standardized Effect Size (Cohen's d): ${res.sampleStats.cohensD.toFixed(3)} (${res.sampleStats.cohensDInterpretation})
  • Glass's Delta (relative to control): ${res.sampleStats.glassDelta.toFixed(3)}
  • Common Language Effect Size (CLES): ${(res.sampleStats.cles*100).toFixed(1)}% probability that a random patient from B scores higher than A

BAYESIAN ANALYSIS (Unequal Variances BEST Updating & Posterior Estimation):
  • Posterior Mean Difference: ${res.difference.mean.toFixed(3)} (SD = ${res.difference.sd.toFixed(3)})
  • 95% Highest Density Interval (HDI): [${res.difference.hdi95.low.toFixed(3)}, ${res.difference.hdi95.high.toFixed(3)}]
  • Probability of Superiority P(μ_B > μ_A | data): ${pctSup}%
  • Posterior Bayesian Cohen's d: Median ${res.cohensD.median.toFixed(3)} [95% HDI: ${res.cohensD.hdi95.low.toFixed(3)} to ${res.cohensD.hdi95.high.toFixed(3)}]
  • Standard Deviation Ratio (σ_B / σ_A): Median ${res.sdRatio.median.toFixed(3)} [95% HDI: ${res.sdRatio.hdi95.low.toFixed(3)} to ${res.sdRatio.hdi95.high.toFixed(3)}]

ROPE Clinical Equivalence Verdict [${res.difference.rope.ropeLow}, ${res.difference.rope.ropeHigh}]:
  • Verdict: ${res.difference.rope.decision}
  • JZS Bayes Factor BF₁₀: ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
      }
    },

    setupSingleProportionModule() {
      const inputs = ["spK", "spN", "spPriorA", "spPriorB", "spBenchmark"];
      inputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener("input", () => this.runSingleProportion());
      });
    },

    runSingleProportion() {
      const k = parseInt(document.getElementById("spK")?.value, 10) || 14;
      const n = parseInt(document.getElementById("spN")?.value, 10) || 210;
      const priorA = parseFloat(document.getElementById("spPriorA")?.value) || 1;
      const priorB = parseFloat(document.getElementById("spPriorB")?.value) || 15;
      const benchmark = parseFloat(document.getElementById("spBenchmark")?.value) || 0.10;

      const res = SingleProportion.analyze({ k, n, priorA, priorB, benchmark });
      this.lastResults.singleProportion = res;

      const elMean = document.getElementById("spKpiMean");
      if (elMean) elMean.textContent = `${(res.posterior.mean * 100).toFixed(2)}%`;

      const elHdi = document.getElementById("spKpiHdi");
      if (elHdi) elHdi.textContent = `[${(res.posterior.hdi95.low * 100).toFixed(2)}%, ${(res.posterior.hdi95.high * 100).toFixed(2)}%]`;

      const elEffectH = document.getElementById("spKpiEffectH");
      if (elEffectH) elEffectH.textContent = `h = ${res.effectH.toFixed(3)}`;

      const elEffectHLabel = document.getElementById("spKpiEffectHLabel");
      if (elEffectHLabel) {
        const absH = Math.abs(res.effectH);
        elEffectHLabel.textContent = absH >= 0.5 ? "Medium/Large Effect" : absH >= 0.2 ? "Small Effect" : "Negligible Effect";
      }

      const elBench = document.getElementById("spKpiProbBench");
      if (elBench) elBench.textContent = `${(res.benchmark.probExceed * 100).toFixed(1)}%`;

      const canvas = document.getElementById("chartSingleProp");
      if (canvas) {
        CanvasCharts.drawTriPlot(canvas, res.triPlotData);
      }
    },

    setupSingleMeanModule() {
      const inputs = ["smMean", "smSd", "smN", "smPriorMean", "smPriorSd", "smBenchmark"];
      inputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener("input", () => this.runSingleMean());
      });
    },

    runSingleMean() {
      const sampleMean = parseFloat(document.getElementById("smMean")?.value) || 120;
      const sampleSd = parseFloat(document.getElementById("smSd")?.value) || 15;
      const sampleN = parseInt(document.getElementById("smN")?.value, 10) || 50;
      const priorMean = parseFloat(document.getElementById("smPriorMean")?.value) || 125;
      const priorSd = parseFloat(document.getElementById("smPriorSd")?.value) || 10;
      const benchmark = parseFloat(document.getElementById("smBenchmark")?.value) || 120;

      const res = SingleMean.analyze({ sampleMean, sampleSd, sampleN, priorMean, priorSd, benchmark });
      this.lastResults.singleMean = res;

      const elMean = document.getElementById("smKpiMean");
      if (elMean) elMean.textContent = res.posterior.mean.toFixed(2);

      const elHdi = document.getElementById("smKpiHdi");
      if (elHdi) elHdi.textContent = `[${res.posterior.cri95.low.toFixed(2)}, ${res.posterior.cri95.high.toFixed(2)}]`;

      const elEffectD = document.getElementById("smKpiEffectD");
      if (elEffectD) elEffectD.textContent = `d = ${res.sample.effectD.toFixed(3)}`;

      const elEffectDLabel = document.getElementById("smKpiEffectDLabel");
      if (elEffectDLabel) {
        const absD = Math.abs(res.sample.effectD);
        elEffectDLabel.textContent = absD >= 0.8 ? "Large Effect" : absD >= 0.5 ? "Medium Effect" : absD >= 0.2 ? "Small Effect" : "Negligible";
      }

      const elBench = document.getElementById("smKpiProbBench");
      if (elBench) elBench.textContent = `${(res.benchmark.probExceed * 100).toFixed(1)}%`;

      const canvas = document.getElementById("chartSingleMean");
      if (canvas) {
        CanvasCharts.drawTriPlot(canvas, res.triPlotData);
      }
    },

    setupDiagnosticModule() {
      const inputs = ["diagPreProb", "diagSens", "diagSpec"];
      inputs.forEach(id => {
        const el = document.getElementById(id);
        if (el) el.addEventListener("input", () => this.runDiagnosticNomogram());
      });
    },

    runDiagnosticNomogram() {
      const preTestProb = parseFloat(document.getElementById("diagPreProb")?.value) || 0.40;
      const sensitivity = parseFloat(document.getElementById("diagSens")?.value) || 0.65;
      const specificity = parseFloat(document.getElementById("diagSpec")?.value) || 0.88;

      const res = DiagnosticNomogram.analyze({ preTestProb, sensitivity, specificity });
      this.lastResults.diagnostic = res;

      const elPLR = document.getElementById("diagKpiPLR");
      if (elPLR) elPLR.textContent = res.plr.toFixed(2);

      const elPLRLabel = document.getElementById("diagKpiPLRLabel");
      if (elPLRLabel) {
        elPLRLabel.textContent = res.plr >= 10 ? "Large Shift (LR+ ≥ 10)" : res.plr >= 5 ? "Moderate Shift" : "Small Shift";
      }

      const elNLR = document.getElementById("diagKpiNLR");
      if (elNLR) elNLR.textContent = res.nlr.toFixed(3);

      const elNLRLabel = document.getElementById("diagKpiNLRLabel");
      if (elNLRLabel) {
        elNLRLabel.textContent = res.nlr <= 0.1 ? "Large Shift (LR- ≤ 0.1)" : res.nlr <= 0.2 ? "Moderate Shift" : "Small Shift";
      }

      const elDOR = document.getElementById("diagKpiDOR");
      if (elDOR) elDOR.textContent = res.dor.toFixed(1);

      const elPostPos = document.getElementById("diagKpiPostPos");
      if (elPostPos) elPostPos.textContent = `${(res.postTestPositive.prob * 100).toFixed(1)}% (95% HDI ${(res.postTestPositive.hdi95.low * 100).toFixed(1)}–${(res.postTestPositive.hdi95.high * 100).toFixed(1)}%)`;

      const elPostNeg = document.getElementById("diagKpiPostNeg");
      if (elPostNeg) elPostNeg.textContent = `${(res.postTestNegative.prob * 100).toFixed(1)}% (95% HDI ${(res.postTestNegative.hdi95.low * 100).toFixed(1)}–${(res.postTestNegative.hdi95.high * 100).toFixed(1)}%)`;

      const container = document.getElementById("nomogramSvgContainer");
      if (container) {
        CanvasCharts.renderFaganNomogram(container, preTestProb, res.plr, res.nlr, res.postTestPositive.prob, res.postTestNegative.prob);
      }
    },

    setupExportActions() {
      const copyPropBtn = document.getElementById("btnCopyReportTwoProp");
      if (copyPropBtn) {
        copyPropBtn.addEventListener("click", () => {
          const text = document.getElementById("reportTwoProp")?.textContent;
          if (text) {
            navigator.clipboard.writeText(text).then(() => {
              copyPropBtn.textContent = "Copied to Clipboard!";
              setTimeout(() => (copyPropBtn.textContent = "Copy APA/Clinical Report"), 2000);
            });
          }
        });
      }

      const copyMeansBtn = document.getElementById("btnCopyReportTwoMeans");
      if (copyMeansBtn) {
        copyMeansBtn.addEventListener("click", () => {
          const text = document.getElementById("reportTwoMeans")?.textContent;
          if (text) {
            navigator.clipboard.writeText(text).then(() => {
              copyMeansBtn.textContent = "Copied to Clipboard!";
              setTimeout(() => (copyMeansBtn.textContent = "Copy APA/Clinical Report"), 2000);
            });
          }
        });
      }

      document.querySelectorAll(".btn-export-png").forEach(btn => {
        btn.addEventListener("click", () => {
          const targetCanvasId = btn.getAttribute("data-target-canvas");
          const canvas = document.getElementById(targetCanvasId);
          if (canvas) {
            CanvasCharts.exportPNG(canvas, `${targetCanvasId}_bayes.png`);
          }
        });
      });
    },

    redrawActiveCharts() {
      if (this.activeTab === "twoProportions") {
        this.runTwoProportions();
      } else if (this.activeTab === "twoMeans") {
        this.runTwoMeans();
      } else if (this.activeTab === "singleProportion") {
        this.runSingleProportion();
      } else if (this.activeTab === "singleMean") {
        this.runSingleMean();
      } else if (this.activeTab === "diagnostic") {
        this.runDiagnosticNomogram();
      }
    }
  };

  if (typeof module !== "undefined" && module.exports) {
    module.exports = {
      SpecialFunctions, Distributions, MCMCSampler, BootstrapEngine,
      TwoProportions, TwoMeans, SingleProportion, SingleMean, DiagnosticNomogram
    };
  }

  if (typeof document !== "undefined") {
    if (document.readyState === "loading") {
      document.addEventListener("DOMContentLoaded", () => App.init());
    } else {
      App.init();
    }
  }

})();
