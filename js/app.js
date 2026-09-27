/**
 * Bayes-Estimation: Main Application Controller
 * Manages tab switching, module synchronization, reactive UI updates, clinical reporting, and preset loading.
 */

import { ThemeManager } from "./ui/theme.js";
import { Presets } from "./ui/presets.js";
import { CanvasCharts } from "./ui/canvas-charts.js";
import { TwoProportions } from "./modules/two-proportions.js";
import { TwoMeans } from "./modules/two-means.js";
import { SingleProportion } from "./modules/single-proportion.js";
import { SingleMean } from "./modules/single-mean.js";
import { DiagnosticNomogram } from "./modules/diagnostic-nomogram.js";

export const App = {
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

    // Run initial calculations
    this.runTwoProportions();
    this.runTwoMeans();
    this.runSingleProportion();
    this.runSingleMean();
    this.runDiagnosticNomogram();

    // Re-draw charts when window resizes
    window.addEventListener("resize", () => {
      this.redrawActiveCharts();
    });

    // Re-draw when theme changes
    const themeBtn = document.getElementById("themeToggleBtn");
    if (themeBtn) {
      themeBtn.addEventListener("click", () => {
        setTimeout(() => this.redrawActiveCharts(), 50);
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
    // Two Proportions Presets
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

    // Two Means Presets
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

    // Single Prop Presets
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

    // Diagnostic Presets
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

  /* ---------------- Two Proportions Module ---------------- */
  setupTwoProportionsModule() {
    const btn = document.getElementById("btnCalcTwoProp");
    if (btn) btn.addEventListener("click", () => this.runTwoProportions());

    const autoInputs = ["propKA", "propNA", "propPriorA", "propKB", "propNB", "propPriorB", "propRopeLow", "propRopeHigh"];
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

    // Update KPI metrics
    const pctSup = (res.difference.probSuperiority * 100).toFixed(1);
    const diffPct = (res.difference.mean * 100).toFixed(2);
    const hdiLow = (res.difference.hdi95.low * 100).toFixed(2);
    const hdiHigh = (res.difference.hdi95.high * 100).toFixed(2);

    document.getElementById("propKpiProbSup").textContent = `${pctSup}%`;
    document.getElementById("propKpiDiff").textContent = `${diffPct > 0 ? "+" : ""}${diffPct}%`;
    document.getElementById("propKpiDiff").className = `kpi-value ${res.difference.mean >= 0 ? "text-success" : "text-danger"}`;
    document.getElementById("propKpiHdi").textContent = `[${hdiLow}%, ${hdiHigh}%]`;
    document.getElementById("propKpiRR").textContent = `${res.relativeRisk.median.toFixed(2)} (95% CrI ${res.relativeRisk.cri95.low.toFixed(2)}–${res.relativeRisk.cri95.high.toFixed(2)})`;
    document.getElementById("propKpiOR").textContent = `${res.oddsRatio.median.toFixed(2)} (95% CrI ${res.oddsRatio.cri95.low.toFixed(2)}–${res.oddsRatio.cri95.high.toFixed(2)})`;
    document.getElementById("propKpiNNT").textContent = res.nnt ? Math.round(res.nnt).toString() : "N/A";

    // ROPE decision badge
    const badge = document.getElementById("propRopeBadge");
    if (badge) {
      badge.textContent = res.difference.rope.decision;
      badge.className = `status-badge ${res.difference.rope.statusClass}`;
    }
    const ropeDetail = document.getElementById("propRopeDetail");
    if (ropeDetail) {
      ropeDetail.textContent = `${res.difference.rope.decisionDetail} (${res.difference.rope.pctIn.toFixed(1)}% in ROPE, ${res.difference.rope.pctAbove.toFixed(1)}% above, ${res.difference.rope.pctBelow.toFixed(1)}% below).`;
    }

    // Bayes Factor
    const bfEl = document.getElementById("propBayesFactor");
    if (bfEl) {
      bfEl.textContent = `BF₁₀ = ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
    }

    // Render Canvas Charts
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

    // Generate Clinical Report
    const reportBox = document.getElementById("reportTwoProp");
    if (reportBox) {
      reportBox.textContent = `CLINICAL & STATISTICAL BAYESIAN REPORT: TWO PROPORTIONS (A/B)
================================================================================
Sample Data:
  • Group A (Control): ${res.groupA.kA} / ${res.groupA.nA} (${(res.groupA.kA/res.groupA.nA*100).toFixed(1)}%)
    Posterior Beta(${res.groupA.posterior.a}, ${res.groupA.posterior.b}): Mean = ${(res.groupA.mean*100).toFixed(1)}%, 95% CrI [${(res.groupA.cri95.low*100).toFixed(1)}%, ${(res.groupA.cri95.high*100).toFixed(1)}%]
  • Group B (Treatment): ${res.groupB.kB} / ${res.groupB.nB} (${(res.groupB.kB/res.groupB.nB*100).toFixed(1)}%)
    Posterior Beta(${res.groupB.posterior.a}, ${res.groupB.posterior.b}): Mean = ${(res.groupB.mean*100).toFixed(1)}%, 95% CrI [${(res.groupB.cri95.low*100).toFixed(1)}%, ${(res.groupB.cri95.high*100).toFixed(1)}%]

Posterior Difference Analysis (Treatment - Control):
  • Posterior Mean Difference (ARR): ${(res.difference.mean*100).toFixed(2)}% (SD = ${(res.difference.sd*100).toFixed(2)}%)
  • 95% Highest Density Interval (HDI): [${hdiLow}%, ${hdiHigh}%]
  • Probability of Superiority P(p_B > p_A | data): ${pctSup}%
  • Relative Risk (RR): ${res.relativeRisk.median.toFixed(3)} [95% HDI: ${res.relativeRisk.hdi95.low.toFixed(3)} to ${res.relativeRisk.hdi95.high.toFixed(3)}]
  • Odds Ratio (OR): ${res.oddsRatio.median.toFixed(3)} [95% HDI: ${res.oddsRatio.hdi95.low.toFixed(3)} to ${res.oddsRatio.hdi95.high.toFixed(3)}]
  • Number Needed to Treat (NNT): ${res.nnt ? res.nnt.toFixed(1) : "Undefined"}

Region of Practical Equivalence (ROPE) [${(res.difference.rope.ropeLow*100).toFixed(1)}%, ${(res.difference.rope.ropeHigh*100).toFixed(1)}%]:
  • ${res.difference.rope.decision}
  • ${res.difference.rope.pctIn.toFixed(1)}% in ROPE | ${res.difference.rope.pctAbove.toFixed(1)}% Superior | ${res.difference.rope.pctBelow.toFixed(1)}% Inferior
  • Bayes Factor BF₁₀: ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
    }
  },

  /* ---------------- Two Means Module ---------------- */
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

    // Update KPIs
    const pctSup = (res.difference.probSuperiority * 100).toFixed(1);
    document.getElementById("meansKpiProbSup").textContent = `${pctSup}%`;
    document.getElementById("meansKpiDiff").textContent = `${res.difference.mean > 0 ? "+" : ""}${res.difference.mean.toFixed(2)}`;
    document.getElementById("meansKpiHdi").textContent = `[${res.difference.hdi95.low.toFixed(2)}, ${res.difference.hdi95.high.toFixed(2)}]`;
    document.getElementById("meansKpiCohensD").textContent = `${res.cohensD.median.toFixed(2)} (95% HDI ${res.cohensD.hdi95.low.toFixed(2)}–${res.cohensD.hdi95.high.toFixed(2)})`;
    document.getElementById("meansKpiSdRatio").textContent = `${res.sdRatio.median.toFixed(2)} (95% HDI ${res.sdRatio.hdi95.low.toFixed(2)}–${res.sdRatio.hdi95.high.toFixed(2)})`;
    document.getElementById("meansKpiBF").textContent = `BF₁₀ = ${res.bayesFactor.bf10.toFixed(2)}`;

    // ROPE badge
    const badge = document.getElementById("meansRopeBadge");
    if (badge) {
      badge.textContent = res.difference.rope.decision;
      badge.className = `status-badge ${res.difference.rope.statusClass}`;
    }
    const ropeDetail = document.getElementById("meansRopeDetail");
    if (ropeDetail) {
      ropeDetail.textContent = `${res.difference.rope.decisionDetail} (${res.difference.rope.pctIn.toFixed(1)}% in ROPE, ${res.difference.rope.pctAbove.toFixed(1)}% above, ${res.difference.rope.pctBelow.toFixed(1)}% below).`;
    }

    // Canvas charts
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

    // Report
    const reportBox = document.getElementById("reportTwoMeans");
    if (reportBox) {
      reportBox.textContent = `CLINICAL & STATISTICAL BAYESIAN REPORT: TWO CONTINUOUS MEANS (BEST)
================================================================================
Sample Statistics:
  • Group A (Control): N = ${res.groupA.n}, Mean = ${res.groupA.mean.toFixed(2)}, SD = ${res.groupA.sd.toFixed(2)} (95% CrI [${res.groupA.cri95.low.toFixed(2)}, ${res.groupA.cri95.high.toFixed(2)}])
  • Group B (Treatment): N = ${res.groupB.n}, Mean = ${res.groupB.mean.toFixed(2)}, SD = ${res.groupB.sd.toFixed(2)} (95% CrI [${res.groupB.cri95.low.toFixed(2)}, ${res.groupB.cri95.high.toFixed(2)}])

Bayesian Mean Difference (μ_B - μ_A):
  • Posterior Mean Difference: ${res.difference.mean.toFixed(3)} (SD = ${res.difference.sd.toFixed(3)})
  • 95% Highest Density Interval (HDI): [${res.difference.hdi95.low.toFixed(3)}, ${res.difference.hdi95.high.toFixed(3)}]
  • Probability of Superiority P(μ_B > μ_A | data): ${pctSup}%
  • Bayesian Cohen's d: ${res.cohensD.median.toFixed(3)} [95% HDI: ${res.cohensD.hdi95.low.toFixed(3)} to ${res.cohensD.hdi95.high.toFixed(3)}]
  • Standard Deviation Ratio (σ_B / σ_A): ${res.sdRatio.median.toFixed(3)} [95% HDI: ${res.sdRatio.hdi95.low.toFixed(3)} to ${res.sdRatio.hdi95.high.toFixed(3)}]

Bayesian Model Comparison & ROPE [${res.difference.rope.ropeLow}, ${res.difference.rope.ropeHigh}]:
  • Decision: ${res.difference.rope.decision}
  • JZS Bayes Factor BF₁₀: ${res.bayesFactor.bf10.toFixed(2)} (${res.bayesFactor.interpretation})`;
    }
  },

  /* ---------------- Single Proportion Module ---------------- */
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

    document.getElementById("spKpiMean").textContent = `${(res.posterior.mean * 100).toFixed(2)}%`;
    document.getElementById("spKpiHdi").textContent = `[${(res.posterior.hdi95.low * 100).toFixed(2)}%, ${(res.posterior.hdi95.high * 100).toFixed(2)}%]`;
    document.getElementById("spKpiProbBench").textContent = `${(res.benchmark.probExceed * 100).toFixed(1)}%`;
    document.getElementById("spKpiPostA").textContent = `Beta(${res.posterior.a}, ${res.posterior.b})`;

    const canvas = document.getElementById("chartSingleProp");
    if (canvas) {
      CanvasCharts.drawTriPlot(canvas, res.triPlotData);
    }
  },

  /* ---------------- Single Mean Module ---------------- */
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

    document.getElementById("smKpiMean").textContent = res.posterior.mean.toFixed(2);
    document.getElementById("smKpiHdi").textContent = `[${res.posterior.cri95.low.toFixed(2)}, ${res.posterior.cri95.high.toFixed(2)}]`;
    document.getElementById("smKpiProbBench").textContent = `${(res.benchmark.probExceed * 100).toFixed(1)}%`;
    document.getElementById("smKpiPostSd").textContent = `SD = ${res.posterior.sd.toFixed(2)}`;

    const canvas = document.getElementById("chartSingleMean");
    if (canvas) {
      CanvasCharts.drawTriPlot(canvas, res.triPlotData);
    }
  },

  /* ---------------- Diagnostic Module ---------------- */
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

    document.getElementById("diagKpiPLR").textContent = res.plr.toFixed(2);
    document.getElementById("diagKpiNLR").textContent = res.nlr.toFixed(3);
    document.getElementById("diagKpiPostPos").textContent = `${(res.postTestPositive.prob * 100).toFixed(1)}% (95% HDI ${(res.postTestPositive.hdi95.low * 100).toFixed(1)}–${(res.postTestPositive.hdi95.high * 100).toFixed(1)}%)`;
    document.getElementById("diagKpiPostNeg").textContent = `${(res.postTestNegative.prob * 100).toFixed(1)}% (95% HDI ${(res.postTestNegative.hdi95.low * 100).toFixed(1)}–${(res.postTestNegative.hdi95.high * 100).toFixed(1)}%)`;

    const container = document.getElementById("nomogramSvgContainer");
    if (container) {
      CanvasCharts.renderFaganNomogram(container, preTestProb, res.plr, res.nlr, res.postTestPositive.prob, res.postTestNegative.prob);
    }
  },

  /* ---------------- Export & Copy Utilities ---------------- */
  setupExportActions() {
    // Copy Two Prop Report
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

    // Copy Two Means Report
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

    // Export PNG buttons
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

// Auto-boot upon DOM ready
document.addEventListener("DOMContentLoaded", () => {
  App.init();
});
