/**
 * Bayes-Estimation: Single Continuous Mean Estimation (Normal-Normal Conjugate Model)
 * Prior, Likelihood, and Posterior distributions with Credible Intervals and Benchmark probabilities.
 */

import { Distributions } from "../math/distributions.js";
import { MCMCSampler } from "../math/mcmc-sampler.js";

export const SingleMean = {
  analyze({
    sampleMean, sampleSd, sampleN,
    priorMean = 0, priorSd = 10,
    benchmark = 0
  }) {
    sampleMean = parseFloat(sampleMean) || 0;
    sampleSd = Math.max(1e-6, parseFloat(sampleSd) || 1);
    sampleN = Math.max(1, parseInt(sampleN, 10) || 1);
    priorMean = parseFloat(priorMean) || 0;
    priorSd = Math.max(1e-6, parseFloat(priorSd) || 10);
    benchmark = parseFloat(benchmark) || 0;

    const sampleSem = sampleSd / Math.sqrt(sampleN);

    // Normal-Normal conjugate precision weighting:
    // tau_prior = 1 / priorSd^2
    // tau_data  = n / sampleSd^2 = 1 / sampleSem^2
    // tau_post  = tau_prior + tau_data
    // mu_post   = (tau_prior * priorMean + tau_data * sampleMean) / tau_post
    // sd_post   = 1 / sqrt(tau_post)
    const tauPrior = 1 / (priorSd * priorSd);
    const tauData = 1 / (sampleSem * sampleSem);
    const tauPost = tauPrior + tauData;

    const postMean = (tauPrior * priorMean + tauData * sampleMean) / tauPost;
    const postSd = 1 / Math.sqrt(tauPost);

    // 95% Credible Interval and HDI (for Normal, symmetric CrI = HDI)
    const zCrit = 1.959963984540054;
    const cri95 = {
      low: postMean - zCrit * postSd,
      high: postMean + zCrit * postSd
    };

    // Probability exceeding benchmark
    const probExceed = 1 - Distributions.normalCDF(benchmark, postMean, postSd);

    // Generate Tri-Plot data grid
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
      const priorDens = Distributions.normalPDF(x, priorMean, priorSd);
      const likDens = Distributions.normalPDF(x, sampleMean, sampleSem);
      const postDens = Distributions.normalPDF(x, postMean, postSd);
      triPlotData.push({ x, prior: priorDens, likelihood: likDens, posterior: postDens });
    }

    return {
      sample: { mean: sampleMean, sd: sampleSd, n: sampleN, sem: sampleSem },
      prior: { mean: priorMean, sd: priorSd },
      posterior: { mean: postMean, sd: postSd, cri95, hdi95: cri95 },
      benchmark: { value: benchmark, probExceed },
      triPlotData
    };
  }
};
