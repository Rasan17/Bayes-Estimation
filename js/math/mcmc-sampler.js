/**
 * Bayes-Estimation: Fast Monte Carlo Sampler & Posterior Analysis Engine
 * High-performance random generators, HDI (Highest Density Interval), ROPE, and Kernel Density Estimation.
 */

export const MCMCSampler = {
  /**
   * Generates a standard normal random variate N(0, 1) via Box-Muller transform
   */
  randomNormal(mean = 0, sd = 1) {
    let u1 = 0, u2 = 0;
    while (u1 === 0) u1 = Math.random();
    while (u2 === 0) u2 = Math.random();
    const z0 = Math.sqrt(-2.0 * Math.log(u1)) * Math.cos(2.0 * Math.PI * u2);
    return mean + z0 * sd;
  },

  /**
   * Generates a Gamma random variate Gamma(shape, scale)
   * Uses Marsaglia and Tsang method (2000) for shape >= 1
   */
  randomGamma(shape, scale = 1) {
    if (shape < 1) {
      // Johnk's or multiplication property: Gamma(a) = Gamma(a+1) * U^(1/a)
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

      if (u < 1 - 0.0331 * x * x * x * x) {
        return d * v * scale;
      }
      if (Math.log(u) < 0.5 * x * x + d * (1 - v + Math.log(v))) {
        return d * v * scale;
      }
    }
  },

  /**
   * Generates a Beta random variate Beta(a, b)
   * via ratio of two independent Gamma random variates: X = G_a / (G_a + G_b)
   */
  randomBeta(a, b) {
    if (a <= 0 || b <= 0) return 0.5;
    const ga = this.randomGamma(a, 1);
    const gb = this.randomGamma(b, 1);
    const sum = ga + gb;
    if (sum === 0) return 0.5;
    return ga / sum;
  },

  /**
   * Generates a Student-t random variate t(df, loc, scale)
   */
  randomStudentT(df, loc = 0, scale = 1) {
    const z = this.randomNormal(0, 1);
    const chi2 = this.randomGamma(df / 2, 2); // Chi-squared(df) = Gamma(df/2, 2)
    return loc + scale * (z / Math.sqrt(chi2 / df));
  },

  /**
   * Computes summary statistics (mean, median, sd, mode) of a numeric sample array
   */
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

    // Sort a copy for median, quantiles, and HDI
    const sorted = new Float64Array(samples).sort();
    const mid = Math.floor(n / 2);
    const median = n % 2 !== 0 ? sorted[mid] : (sorted[mid - 1] + sorted[mid]) / 2;

    // Estimate Mode via KDE peak
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

  /**
   * Computes Highest Density Interval (HDI) from sorted samples
   * Kruschke method: shortest interval enclosing credibleMass (e.g. 0.95)
   */
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

  /**
   * Computes Equal-Tailed Credible Interval (CrI) from sorted samples
   */
  computeCrI(sortedSamples, credibleMass = 0.95) {
    const n = sortedSamples.length;
    const alpha = 1 - credibleMass;
    const lowIdx = Math.max(0, Math.floor((alpha / 2) * n));
    const highIdx = Math.min(n - 1, Math.floor((1 - alpha / 2) * n));
    return {
      low: sortedSamples[lowIdx],
      high: sortedSamples[highIdx],
      mass: credibleMass
    };
  },

  /**
   * Evaluates Region of Practical Equivalence (ROPE)
   * @param {Float64Array|Array<number>} sortedSamples
   * @param {number} ropeLow - Lower ROPE bound (e.g., -0.05)
   * @param {number} ropeHigh - Upper ROPE bound (e.g., +0.05)
   * @param {Object} hdi95 - The 95% HDI { low, high }
   */
  evaluateROPE(sortedSamples, ropeLow, ropeHigh, hdi95) {
    const n = sortedSamples.length;
    let countBelow = 0;
    let countIn = 0;
    let countAbove = 0;

    for (let i = 0; i < n; i++) {
      const val = sortedSamples[i];
      if (val < ropeLow) countBelow++;
      else if (val > ropeHigh) countAbove++;
      else countIn++;
    }

    const pctBelow = (countBelow / n) * 100;
    const pctIn = (countIn / n) * 100;
    const pctAbove = (countAbove / n) * 100;

    // Decision rule (Kruschke Bayesian decision framework)
    let decision = "Undecided / Data Inconclusive";
    let decisionDetail = "The 95% HDI overlaps the ROPE boundary. More samples or narrower clinical equivalence margins required.";
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

    return {
      ropeLow,
      ropeHigh,
      pctBelow,
      pctIn,
      pctAbove,
      decision,
      decisionDetail,
      statusClass
    };
  },

  /**
   * Gaussian Kernel Density Estimation (KDE) with Silverman's rule of thumb bandwidth
   */
  estimateKDE(sortedSamples, numGridPoints = 120, xMin = null, xMax = null) {
    const n = sortedSamples.length;
    if (n === 0) return [];

    let sum = 0;
    for (let i = 0; i < n; i++) sum += sortedSamples[i];
    const mean = sum / n;

    let sqDiff = 0;
    for (let i = 0; i < n; i++) sqDiff += Math.pow(sortedSamples[i] - mean, 2);
    const sd = Math.sqrt(sqDiff / (n - 1)) || 1e-4;

    // Interquartile range (IQR)
    const q25 = sortedSamples[Math.floor(0.25 * n)];
    const q75 = sortedSamples[Math.floor(0.75 * n)];
    const iqr = (q75 - q25) || sd;

    // Silverman's rule: h = 0.9 * min(SD, IQR/1.34) * n^(-1/5)
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

  /**
   * Savage-Dickey Density Ratio at point null value (e.g. delta = 0)
   * BF_10 = Prior density at null / Posterior density at null
   */
  savageDickeyBF(priorDensityAtNull, posteriorSortedSamples, nullVal = 0) {
    const n = posteriorSortedSamples.length;
    if (n === 0) return { bf10: 1, bf01: 1, interpretation: "Insufficient data" };

    // Density of posterior at null value via Gaussian KDE
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

    return {
      bf10,
      bf01,
      priorDensityAtNull,
      posteriorDensityAtNull,
      interpretation
    };
  }
};
