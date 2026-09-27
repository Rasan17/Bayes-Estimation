/**
 * Bayes-Estimation: Single Proportion Estimation (Beta-Binomial Conjugate Model)
 * Live Tri-plot (Prior, Likelihood, Posterior), Credible Interval, HDI, and Benchmark testing.
 */

import { Distributions } from "../math/distributions.js";
import { MCMCSampler } from "../math/mcmc-sampler.js";

export const SingleProportion = {
  analyze({ k, n, priorA = 1, priorB = 1, benchmark = 0.50, numSamples = 30000 }) {
    k = Math.max(0, parseInt(k, 10) || 0);
    n = Math.max(1, parseInt(n, 10) || 1);
    if (k > n) k = n;

    priorA = Math.max(0.001, parseFloat(priorA) || 1);
    priorB = Math.max(0.001, parseFloat(priorB) || 1);
    benchmark = Math.max(0, Math.min(1, parseFloat(benchmark) || 0.50));

    // Conjugate Posterior
    const postA = priorA + k;
    const postB = priorB + (n - k);

    // Prior moments
    const priorMean = priorA / (priorA + priorB);
    const priorVar = (priorA * priorB) / (Math.pow(priorA + priorB, 2) * (priorA + priorB + 1));
    const priorSd = Math.sqrt(priorVar);

    // Posterior moments
    const postMean = postA / (postA + postB);
    const postVar = (postA * postB) / (Math.pow(postA + postB, 2) * (postA + postB + 1));
    const postSd = Math.sqrt(postVar);
    const postMode = (postA > 1 && postB > 1) ? (postA - 1) / (postA + postB - 2) : postMean;
    const postMedian = Distributions.betaQuantile(0.5, postA, postB);

    // Credible Intervals
    const cri95 = {
      low: Distributions.betaQuantile(0.025, postA, postB),
      high: Distributions.betaQuantile(0.975, postA, postB)
    };
    const cri90 = {
      low: Distributions.betaQuantile(0.05, postA, postB),
      high: Distributions.betaQuantile(0.95, postA, postB)
    };

    // HDI via sample simulation
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

    // Grid curves for Tri-Plot (Prior, Normalized Likelihood, Posterior)
    const gridPoints = 200;
    const triPlotData = [];
    let maxLikelihood = 0;

    // Find max binomial likelihood for scaling: L(p) = p^k * (1-p)^(n-k)
    // Mode of likelihood is k/n
    const mP = k / n;
    const maxLikVal = Math.pow(mP, k) * Math.pow(1 - mP, n - k) || 1e-10;

    // Peak of posterior for visual normalization
    const maxPostDensity = Math.max(1, Distributions.betaPDF(postMode, postA, postB));

    for (let i = 0; i <= gridPoints; i++) {
      const p = i / gridPoints;
      const priorDens = Distributions.betaPDF(p, priorA, priorB);
      const postDens = Distributions.betaPDF(p, postA, postB);

      // Normalized likelihood to match posterior scale
      const rawLik = (p === 0 && k === 0) || (p === 1 && k === n) ? 1 :
                    (p > 0 && p < 1) ? Math.pow(p, k) * Math.pow(1 - p, n - k) : 0;
      const normLik = (rawLik / maxLikVal) * maxPostDensity;

      triPlotData.push({
        p,
        prior: priorDens,
        likelihood: normLik,
        posterior: postDens
      });
    }

    return {
      k, n,
      prior: { a: priorA, b: priorB, mean: priorMean, sd: priorSd },
      posterior: {
        a: postA, b: postB,
        mean: postMean,
        median: postMedian,
        mode: postMode,
        sd: postSd,
        cri95,
        cri90,
        hdi95
      },
      benchmark: {
        value: benchmark,
        probExceed: probExceedBenchmark
      },
      triPlotData
    };
  }
};
