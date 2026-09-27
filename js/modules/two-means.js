/**
 * Bayes-Estimation: Bayesian Two-Group Continuous Means Module (BEST / Bayesian t-Test)
 * Unequal variances Bayesian estimation for continuous clinical endpoints.
 * Calculates difference in means, Bayesian Cohen's d, SD ratio, HDI, ROPE, and Bayes Factors.
 */

import { Distributions } from "../math/distributions.js";
import { MCMCSampler } from "../math/mcmc-sampler.js";

export const TwoMeans = {
  /**
   * Parses raw numeric text (comma, space, or newline delimited) into a clean array of numbers
   */
  parseData(text) {
    if (!text || typeof text !== "string") return [];
    return text
      .split(/[\s,;\t\n]+/)
      .map(v => parseFloat(v.trim()))
      .filter(v => !isNaN(v) && isFinite(v));
  },

  /**
   * Computes sample statistics from raw array
   */
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

  /**
   * Analyzes two continuous groups
   */
  analyze(params) {
    let {
      dataA, dataB,
      nA, meanA, sdA,
      nB, meanB, sdB,
      ropeLow = -0.5, ropeHigh = 0.5,
      priorScale = 0.707,
      numSamples = 50000
    } = params;

    // Use raw data if provided
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

    const dfA = nA - 1;
    const dfB = nB - 1;
    const semA = sdA / Math.sqrt(nA);
    const semB = sdB / Math.sqrt(nB);

    const tCritA = Distributions.studentTQuantile(0.975, dfA, 0, 1);
    const tCritB = Distributions.studentTQuantile(0.975, dfB, 0, 1);

    const cri95A = { low: meanA - tCritA * semA, high: meanA + tCritA * semA };
    const cri95B = { low: meanB - tCritB * semB, high: meanB + tCritB * semB };

    // Joint Monte Carlo sampling of muA, muB, sigmaA, sigmaB
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

      const pooledSd = Math.sqrt((varA_sample + varB_sample) / 2);
      cohenSamples[i] = diff / pooledSd;
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
    const pooledVarClassic = ((nA - 1) * sdA * sdA + (nB - 1) * sdB * sdB) / (nA + nB - 2);
    const tClassic = (meanB - meanA) / Math.sqrt(pooledVarClassic * (1 / nA + 1 / nB));
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
    else if (bf10 > 30) bfInterpretation = "Very strong evidence for difference (H₁) ";
    else if (bf10 > 10) bfInterpretation = "Strong evidence for difference (H₁)";
    else if (bf10 > 3) bfInterpretation = "Moderate evidence for difference (H₁)";
    else if (bf10 > 1) bfInterpretation = "Anecdotal evidence for difference (H₁)";
    else if (bf10 === 1) bfInterpretation = "No evidence either way";
    else if (bf10 > 1 / 3) bfInterpretation = "Anecdotal evidence for equality (H₀)";
    else if (bf10 > 1 / 10) bfInterpretation = "Moderate evidence for equality (H₀)";
    else if (bf10 > 1 / 30) bfInterpretation = "Strong evidence for equality (H₀)";
    else bfInterpretation = "Extreme evidence for equality (H₀)";

    return {
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
      sdRatio: {
        mean: sdRatioSummary.mean,
        median: sdRatioSummary.median,
        hdi95: sdRatioHdi95
      },
      bayesFactor: {
        bf10,
        bf01,
        priorScale,
        tClassic,
        interpretation: bfInterpretation
      },
      rawData: {
        dataA: dataA || null,
        dataB: dataB || null
      }
    };
  }
};
