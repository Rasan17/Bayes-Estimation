/**
 * Bayes-Estimation: Bayesian Two-Group Proportions & A/B Testing Module
 * Exact Beta-Binomial conjugate updating, Monte Carlo joint posterior simulation,
 * Difference in proportions, Relative Risk, Odds Ratio, NNT, HDI, ROPE, and Bayes Factors.
 */

import { SpecialFunctions } from "../math/special-functions.js";
import { Distributions } from "../math/distributions.js";
import { MCMCSampler } from "../math/mcmc-sampler.js";

export const TwoProportions = {
  /**
   * Evaluates two-group Bayesian proportion comparison
   * @param {Object} params
   * @param {number} params.kA - Successes in Group A (Control)
   * @param {number} params.nA - Total trials in Group A
   * @param {number} params.aA - Prior alpha for Group A (default 1)
   * @param {number} params.bA - Prior beta for Group A (default 1)
   * @param {number} params.kB - Successes in Group B (Treatment)
   * @param {number} params.nB - Total trials in Group B
   * @param {number} params.aB - Prior alpha for Group B (default 1)
   * @param {number} params.bB - Prior beta for Group B (default 1)
   * @param {number} params.ropeLow - Lower ROPE threshold for difference (default -0.05)
   * @param {number} params.ropeHigh - Upper ROPE threshold for difference (default +0.05)
   * @param {number} params.numSamples - Monte Carlo iterations (default 50000)
   */
  analyze({
    kA, nA, aA = 1, bA = 1,
    kB, nB, aB = 1, bB = 1,
    ropeLow = -0.05, ropeHigh = 0.05,
    numSamples = 50000
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

    // Posterior Beta parameters
    const postAA = aA + kA;
    const postBA = bA + (nA - kA);
    const postAB = aB + kB;
    const postBB = bB + (nB - kB);

    // Analytical posterior moments for Group A
    const meanA = postAA / (postAA + postBA);
    const varA = (postAA * postBA) / (Math.pow(postAA + postBA, 2) * (postAA + postBA + 1));
    const sdA = Math.sqrt(varA);
    const modeA = (postAA > 1 && postBA > 1) ? (postAA - 1) / (postAA + postBA - 2) : meanA;
    const cri95A = {
      low: Distributions.betaQuantile(0.025, postAA, postBA),
      high: Distributions.betaQuantile(0.975, postAA, postBA)
    };

    // Analytical posterior moments for Group B
    const meanB = postAB / (postAB + postBB);
    const varB = (postAB * postBB) / (Math.pow(postAB + postBB, 2) * (postAB + postBB + 1));
    const sdB = Math.sqrt(varB);
    const modeB = (postAB > 1 && postBB > 1) ? (postAB - 1) / (postAB + postBB - 2) : meanB;
    const cri95B = {
      low: Distributions.betaQuantile(0.025, postAB, postBB),
      high: Distributions.betaQuantile(0.975, postAB, postBB)
    };

    // Monte Carlo simulation for differences, RR, and OR
    const diffSamples = new Float64Array(numSamples);
    const rrSamples = new Float64Array(numSamples);
    const orSamples = new Float64Array(numSamples);

    let superiorCount = 0; // pB > pA
    let countAboveZero = 0;

    for (let i = 0; i < numSamples; i++) {
      const pA = MCMCSampler.randomBeta(postAA, postBA);
      const pB = MCMCSampler.randomBeta(postAB, postBB);

      const diff = pB - pA;
      diffSamples[i] = diff;

      if (diff > 0) {
        superiorCount++;
        countAboveZero++;
      }

      // Avoid division by zero
      const safePA = Math.max(1e-7, Math.min(1 - 1e-7, pA));
      const safePB = Math.max(1e-7, Math.min(1 - 1e-7, pB));

      rrSamples[i] = safePB / safePA;
      orSamples[i] = (safePB / (1 - safePB)) / (safePA / (1 - safePA));
    }

    // Probability of Superiority
    const probSuperiority = superiorCount / numSamples;
    const probInferiority = 1 - probSuperiority;

    // Difference posterior analysis
    const diffSummary = MCMCSampler.summarize(diffSamples);
    const diffHdi95 = MCMCSampler.computeHDI(diffSummary.sorted, 0.95);
    const diffHdi90 = MCMCSampler.computeHDI(diffSummary.sorted, 0.90);
    const diffCri95 = MCMCSampler.computeCrI(diffSummary.sorted, 0.95);

    // ROPE evaluation for difference
    const ropeAnalysis = MCMCSampler.evaluateROPE(diffSummary.sorted, ropeLow, ropeHigh, diffHdi95);

    // Relative Risk analysis
    const rrSummary = MCMCSampler.summarize(rrSamples);
    const rrHdi95 = MCMCSampler.computeHDI(rrSummary.sorted, 0.95);
    const rrCri95 = MCMCSampler.computeCrI(rrSummary.sorted, 0.95);

    // Odds Ratio analysis
    const orSummary = MCMCSampler.summarize(orSamples);
    const orHdi95 = MCMCSampler.computeHDI(orSummary.sorted, 0.95);
    const orCri95 = MCMCSampler.computeCrI(orSummary.sorted, 0.95);

    // Number Needed to Treat (NNT = 1 / |mean diff|)
    const absDiff = Math.abs(diffSummary.mean);
    const nnt = absDiff > 1e-5 ? (1 / absDiff) : Infinity;

    // Prior density at null delta = 0 for Savage-Dickey Bayes Factor
    // Under independent Beta(aA, bA) and Beta(aB, bB), prior density of delta=0 can be sampled or computed
    let priorCountAroundZero = 0;
    const priorSimN = 25000;
    const priorDeltaSamples = new Float64Array(priorSimN);
    for (let i = 0; i < priorSimN; i++) {
      const priorA = MCMCSampler.randomBeta(aA, bA);
      const priorB = MCMCSampler.randomBeta(aB, bB);
      priorDeltaSamples[i] = priorB - priorA;
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
      groupA: {
        kA, nA,
        prior: { a: aA, b: bA },
        posterior: { a: postAA, b: postBA },
        mean: meanA,
        median: Distributions.betaQuantile(0.5, postAA, postBA),
        mode: modeA,
        sd: sdA,
        cri95: cri95A
      },
      groupB: {
        kB, nB,
        prior: { a: aB, b: bB },
        posterior: { a: postAB, b: postBB },
        mean: meanB,
        median: Distributions.betaQuantile(0.5, postAB, postBB),
        mode: modeB,
        sd: sdB,
        cri95: cri95B
      },
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
      relativeRisk: {
        mean: rrSummary.mean,
        median: rrSummary.median,
        sd: rrSummary.sd,
        hdi95: rrHdi95,
        cri95: rrCri95
      },
      oddsRatio: {
        mean: orSummary.mean,
        median: orSummary.median,
        sd: orSummary.sd,
        hdi95: orHdi95,
        cri95: orCri95
      },
      nnt: isFinite(nnt) ? nnt : null,
      bayesFactor,
      samplesCount: numSamples
    };
  }
};
