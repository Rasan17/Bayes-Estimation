/**
 * Bayes-Estimation: Clinical Diagnostic Bayesian Reasoning & Fagan Nomogram
 * Computes pre/post-test odds, positive/negative likelihood ratios, and generates
 * the Fagan Nomogram alignment coordinates and posterior uncertainty intervals.
 */

import { Distributions } from "../math/distributions.js";
import { MCMCSampler } from "../math/mcmc-sampler.js";

export const DiagnosticNomogram = {
  /**
   * Analyzes clinical diagnostic test probabilities
   * @param {Object} params
   * @param {number} params.preTestProb - Pre-test disease probability (0 to 1)
   * @param {number} params.sensitivity - Sensitivity (0 to 1)
   * @param {number} params.specificity - Specificity (0 to 1)
   * @param {number} [params.sampleN_disease] - Sample size for disease cohort (for uncertainty)
   * @param {number} [params.sampleN_healthy] - Sample size for healthy cohort
   */
  analyze({
    preTestProb = 0.20,
    sensitivity = 0.85,
    specificity = 0.90,
    sampleN_disease = 100,
    sampleN_healthy = 100
  }) {
    preTestProb = Math.max(0.001, Math.min(0.999, parseFloat(preTestProb) || 0.20));
    sensitivity = Math.max(0.001, Math.min(0.999, parseFloat(sensitivity) || 0.85));
    specificity = Math.max(0.001, Math.min(0.999, parseFloat(specificity) || 0.90));

    // Likelihood Ratios
    const plr = sensitivity / (1 - specificity);
    const nlr = (1 - sensitivity) / specificity;

    // Bayes' Rule in Odds Form: Post-test Odds = Pre-test Odds * LR
    const preTestOdds = preTestProb / (1 - preTestProb);

    const postOddsPos = preTestOdds * plr;
    const postProbPos = postOddsPos / (1 + postOddsPos);

    const postOddsNeg = preTestOdds * nlr;
    const postProbNeg = postOddsNeg / (1 + postOddsNeg);

    // Monte Carlo uncertainty estimation for post-test probabilities
    const simN = 10000;
    const postPosSamples = new Float64Array(simN);
    const postNegSamples = new Float64Array(simN);

    // Beta distributions representing uncertainty in sensitivity and specificity
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
      preTestProb,
      preTestOdds,
      sensitivity,
      specificity,
      plr,
      nlr,
      postTestPositive: {
        prob: postProbPos,
        hdi95: posHdi95,
        kde: MCMCSampler.estimateKDE(posSummary.sorted, 60)
      },
      postTestNegative: {
        prob: postProbNeg,
        hdi95: negHdi95,
        kde: MCMCSampler.estimateKDE(negSummary.sorted, 60)
      }
    };
  }
};
