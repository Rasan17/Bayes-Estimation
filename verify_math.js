/**
 * Bayes-Estimation: Mathematical Verification & Benchmark Test Suite
 * Validates special functions, distributions, MCMC sampler, HDI, ROPE, and statistical modules.
 */

import { SpecialFunctions } from "./js/math/special-functions.js";
import { Distributions } from "./js/math/distributions.js";
import { MCMCSampler } from "./js/math/mcmc-sampler.js";
import { TwoProportions } from "./js/modules/two-proportions.js";
import { TwoMeans } from "./js/modules/two-means.js";
import { SingleProportion } from "./js/modules/single-proportion.js";
import { SingleMean } from "./js/modules/single-mean.js";
import { DiagnosticNomogram } from "./js/modules/diagnostic-nomogram.js";

let testsRun = 0;
let testsPassed = 0;

function assert(condition, message) {
  testsRun++;
  if (condition) {
    testsPassed++;
    console.log(`  ✓ PASS: ${message}`);
  } else {
    console.error(`  ✗ FAIL: ${message}`);
    process.exitCode = 1;
  }
}

function assertClose(actual, expected, tol, message) {
  const diff = Math.abs(actual - expected);
  assert(diff <= tol, `${message} (Expected ~${expected}, Got ${actual}, Diff: ${diff.toExponential(3)})`);
}

console.log("================================================================================");
console.log("               BAYES-ESTIMATION MATHEMATICAL VERIFICATION SUITE                  ");
console.log("================================================================================\n");

// 1. Special Functions
console.log("1. Special Functions & Approximations:");
// lnGamma(5) = ln(24) = 3.1780538303479458
assertClose(SpecialFunctions.lnGamma(5), Math.log(24), 1e-10, "lnGamma(5) equals ln(24)");
// Gamma(0.5) = sqrt(pi) = 1.772453850905516
assertClose(SpecialFunctions.gamma(0.5), Math.sqrt(Math.PI), 1e-10, "gamma(0.5) equals sqrt(pi)");
// Incomplete Beta incBeta(0.5, 2, 2) for symmetric Beta(2,2) at median 0.5 must be 0.5
assertClose(SpecialFunctions.incBeta(0.5, 2, 2), 0.5, 1e-10, "incBeta(0.5, 2, 2) equals 0.5");
// Error function erf(1) ~ 0.84270079
assertClose(SpecialFunctions.erf(1), 0.84270079, 1e-6, "erf(1) matches standard value");
// Inverse normal CDF at 0.975 must be ~ 1.95996
assertClose(SpecialFunctions.invNormalCDF(0.975), 1.959963984540054, 1e-5, "invNormalCDF(0.975) is ~ 1.960");

// 2. Distributions
console.log("\n2. Probability Distributions:");
// Beta(2, 2) PDF at 0.5 = 6 * 0.5 * 0.5 = 1.5
assertClose(Distributions.betaPDF(0.5, 2, 2), 1.5, 1e-10, "Beta(2,2) PDF at mode 0.5 is 1.5");
// Beta(2, 2) CDF at 0.5 is 0.5
assertClose(Distributions.betaCDF(0.5, 2, 2), 0.5, 1e-10, "Beta(2,2) CDF at 0.5 is 0.5");
// Normal PDF at mean 0 with sd 1 is 1/sqrt(2*pi) ~ 0.39894228
assertClose(Distributions.normalPDF(0, 0, 1), 1 / Math.sqrt(2 * Math.PI), 1e-10, "Standard normal PDF peak is 0.3989");
// Student-t CDF at t=0 is 0.5
assertClose(Distributions.studentTCDF(0, 10), 0.5, 1e-10, "Student-t CDF at t=0 is 0.5");

// 3. MCMC Sampler & Highest Density Interval (HDI)
console.log("\n3. MCMC Sampler & Highest Density Interval:");
const nSamples = 10000;
const normalSamples = new Float64Array(nSamples);
for (let i = 0; i < nSamples; i++) {
  normalSamples[i] = MCMCSampler.randomNormal(10, 2);
}
const normalSorted = normalSamples.sort();
const hdi95 = MCMCSampler.computeHDI(normalSorted, 0.95);
// For N(10, 2), 95% interval should be approx 10 - 1.96*2 = 6.08 to 10 + 1.96*2 = 13.92
assertClose(hdi95.low, 10 - 1.96 * 2, 0.3, "95% HDI lower bound for N(10, 2) ~ 6.08");
assertClose(hdi95.high, 10 + 1.96 * 2, 0.3, "95% HDI upper bound for N(10, 2) ~ 13.92");
assert(hdi95.width < 8.2 && hdi95.width > 7.4, "HDI width spans expected interval");

// 4. Two Proportions Module
console.log("\n4. Two-Group Proportions (A/B Testing) Module:");
const propResult = TwoProportions.analyze({
  kA: 20, nA: 100, // 20%
  kB: 40, nB: 100, // 40%
  aA: 1, bA: 1,
  aB: 1, bB: 1,
  ropeLow: -0.05, ropeHigh: 0.05,
  numSamples: 20000
});
assertClose(propResult.groupA.mean, 21 / 102, 0.01, "Group A posterior mean matches Beta(21, 81)");
assertClose(propResult.groupB.mean, 41 / 102, 0.01, "Group B posterior mean matches Beta(41, 61)");
assertClose(propResult.difference.mean, (41 - 21) / 102, 0.02, "Posterior difference mean is ~ 0.196");
assert(propResult.difference.probSuperiority > 0.99, "Treatment is superior with > 99% probability");
assertClose(propResult.relativeRisk.median, 2.0, 0.25, "Relative Risk is ~ 2.0");
assertClose(propResult.nnt, 5.1, 1.0, "NNT is ~ 5.1");
assert(propResult.difference.rope.decision === "Superiority Established", "ROPE decision identifies superiority");

// 5. Two Means Module
console.log("\n5. Two Continuous Means (BEST) Module:");
const meansResult = TwoMeans.analyze({
  meanA: 100, sdA: 15, nA: 50,
  meanB: 110, sdB: 15, nB: 50,
  ropeLow: -2.0, ropeHigh: 2.0,
  numSamples: 20000
});
assertClose(meansResult.difference.mean, 10.0, 0.5, "Mean difference is ~ 10.0");
assertClose(meansResult.cohensD.median, 10.0 / 15.0, 0.15, "Bayesian Cohen's d is ~ 0.67");
assert(meansResult.difference.probSuperiority > 0.99, "Treatment superiority probability is > 99%");
assert(meansResult.difference.rope.decision === "Superiority Established", "ROPE confirms superiority outside [-2, 2]");

// 6. Clinical Diagnostic Nomogram
console.log("\n6. Clinical Diagnostic Bayes & Fagan Nomogram:");
const diagResult = DiagnosticNomogram.analyze({
  preTestProb: 0.20,
  sensitivity: 0.80,
  specificity: 0.90
});
// LR+ = 0.80 / 0.10 = 8.0
assertClose(diagResult.plr, 8.0, 1e-8, "Positive Likelihood Ratio (LR+) is 8.0");
// LR- = 0.20 / 0.90 = 0.222
assertClose(diagResult.nlr, 0.222222, 1e-4, "Negative Likelihood Ratio (LR-) is 0.222");
// Pre-test odds = 0.20 / 0.80 = 0.25. Post odds(+) = 0.25 * 8 = 2.0 -> Post prob(+) = 2/3 = 0.6667
assertClose(diagResult.postTestPositive.prob, 2 / 3, 1e-4, "Post-test probability given positive test is ~ 66.7%");
// Post odds(-) = 0.25 * (2/9) = 1/18 -> Post prob(-) = (1/18) / (19/18) = 1/19 ~ 0.0526
assertClose(diagResult.postTestNegative.prob, 1 / 19, 1e-4, "Post-test probability given negative test is ~ 5.26%");

console.log("\n================================================================================");
console.log(`Results: ${testsPassed} / ${testsRun} tests passed (${Math.round(testsPassed/testsRun*100)}%)`);
console.log("================================================================================\n");
