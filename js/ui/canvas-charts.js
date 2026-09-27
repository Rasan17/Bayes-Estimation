/**
 * Bayes-Estimation: High-DPI Interactive Canvas & SVG Charting Engine
 * Custom rendering for Bayesian posterior distributions, HDI regions, ROPE bands, Tri-plots, and Fagan Nomograms.
 */

export const CanvasCharts = {
  /**
   * Sets up high-DPI canvas resolution
   */
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

  /**
   * Formats a number to clean decimal places
   */
  formatNum(n, d = 3) {
    if (typeof n !== "number" || isNaN(n)) return "--";
    return Number(n.toFixed(d)).toString();
  },

  /**
   * Draws Posterior Difference Distribution with 95% HDI and ROPE
   * @param {HTMLCanvasElement} canvas
   * @param {Array<{x: number, density: number}>} kde - Grid points
   * @param {Object} hdi95 - { low, high }
   * @param {Object} [rope] - { ropeLow, ropeHigh, pctBelow, pctIn, pctAbove }
   * @param {Object} [pointEstimates] - { mean, mode }
   * @param {Object} [options]
   */
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

    // Determine domain bounds
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
    yMax = (yMax || 1) * 1.15; // 15% top margin

    // Coordinate transforms
    const getX = x => padding.left + ((x - xMin) / (xMax - xMin)) * plotW;
    const getY = y => padding.top + plotH - (y / yMax) * plotH;

    ctx.clearRect(0, 0, width, height);

    // 1. Draw ROPE band if specified
    if (rope && rope.ropeLow !== undefined && rope.ropeHigh !== undefined) {
      const rx1 = Math.max(padding.left, getX(rope.ropeLow));
      const rx2 = Math.min(padding.left + plotW, getX(rope.ropeHigh));
      const rw = rx2 - rx1;

      if (rw > 0) {
        ctx.fillStyle = ropeColor;
        ctx.fillRect(rx1, padding.top, rw, plotH);

        // ROPE vertical boundary lines
        ctx.beginPath();
        ctx.strokeStyle = ropeBorder;
        ctx.lineWidth = 1.5;
        ctx.setLineDash([4, 3]);
        ctx.moveTo(rx1, padding.top); ctx.lineTo(rx1, padding.top + plotH);
        ctx.moveTo(rx2, padding.top); ctx.lineTo(rx2, padding.top + plotH);
        ctx.stroke();
        ctx.setLineDash([]);

        // ROPE label
        ctx.fillStyle = ropeBorder;
        ctx.font = "600 11px system-ui, -apple-system, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`ROPE [${this.formatNum(rope.ropeLow, 2)}, ${this.formatNum(rope.ropeHigh, 2)}]`, (rx1 + rx2) / 2, padding.top + 16);
      }
    }

    // 2. Draw zero reference line (Null Difference = 0)
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

    // 3. Shaded 95% HDI region under curve
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

        // HDI base line with end markers
        const hdiY = getY(0) - 2;
        const hx1 = getX(hdi95.low);
        const hx2 = getX(hdi95.high);

        ctx.beginPath();
        ctx.strokeStyle = textBright;
        ctx.lineWidth = 3;
        ctx.moveTo(hx1, hdiY);
        ctx.lineTo(hx2, hdiY);
        // Vertical ticks
        ctx.moveTo(hx1, hdiY - 6); ctx.lineTo(hx1, hdiY + 6);
        ctx.moveTo(hx2, hdiY - 6); ctx.lineTo(hx2, hdiY + 6);
        ctx.stroke();

        // HDI text labels
        ctx.fillStyle = textBright;
        ctx.font = "bold 11px system-ui, -apple-system, sans-serif";
        ctx.textAlign = "center";
        ctx.fillText(`95% HDI: [${this.formatNum(hdi95.low, 3)}, ${this.formatNum(hdi95.high, 3)}]`, (hx1 + hx2) / 2, hdiY - 12);
      }
    }

    // 4. Fill full density curve
    ctx.beginPath();
    ctx.moveTo(getX(kde[0].x), getY(0));
    for (const pt of kde) {
      ctx.lineTo(getX(pt.x), getY(pt.density));
    }
    ctx.lineTo(getX(kde[kde.length - 1].x), getY(0));
    ctx.closePath();
    ctx.fillStyle = fillColor;
    ctx.fill();

    // 5. Outline density curve
    ctx.beginPath();
    ctx.moveTo(getX(kde[0].x), getY(kde[0].density));
    for (let i = 1; i < kde.length; i++) {
      ctx.lineTo(getX(kde[i].x), getY(kde[i].density));
    }
    ctx.strokeStyle = lineColor;
    ctx.lineWidth = 2.5;
    ctx.stroke();

    // 6. Mean and Mode markers
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

    // 7. Axes and Ticks
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 1;

    // X Axis
    ctx.beginPath();
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

    // Axis label
    ctx.fillStyle = textBright;
    ctx.font = "600 12px system-ui, -apple-system, sans-serif";
    ctx.fillText(options.xLabel || "Parameter Difference (Group B − Group A)", padding.left + plotW / 2, padding.top + plotH + 46);
  },

  /**
   * Draws Dual Group Density Distributions (Group A vs Group B)
   */
  drawGroupDensities(canvas, seriesList, options = {}) {
    if (!canvas || !seriesList || seriesList.length === 0) return;
    const { ctx, width, height } = this.setupCanvas(canvas);

    const isDark = document.documentElement.getAttribute("data-theme") === "dark";
    const textColor = isDark ? "#94a3b8" : "#475569";
    const textBright = isDark ? "#f8fafc" : "#0f172a";
    const gridColor = isDark ? "rgba(148, 163, 184, 0.12)" : "rgba(100, 116, 139, 0.15)";

    const padding = { top: 35, right: 30, bottom: 45, left: 55 };
    const plotW = width - padding.left - padding.right;
    const plotH = height - padding.top - padding.bottom;

    // Find global X and Y bounds
    let xMin = Infinity, xMax = -Infinity, yMax = 0;
    for (const s of seriesList) {
      for (const pt of s.data) {
        if (pt.x < xMin) xMin = pt.x;
        if (pt.x > xMax) xMax = pt.x;
        if (pt.y > yMax) yMax = pt.y;
      }
    }
    const xPad = (xMax - xMin) * 0.08 || 0.05;
    xMin -= xPad;
    xMax += xPad;
    yMax *= 1.15;

    const getX = x => padding.left + ((x - xMin) / (xMax - xMin)) * plotW;
    const getY = y => padding.top + plotH - (y / yMax) * plotH;

    ctx.clearRect(0, 0, width, height);

    // Draw grid lines
    ctx.strokeStyle = gridColor;
    ctx.lineWidth = 1;
    for (let i = 0; i <= 4; i++) {
      const yVal = (i / 4) * yMax;
      const py = getY(yVal);
      ctx.beginPath();
      ctx.moveTo(padding.left, py);
      ctx.lineTo(padding.left + plotW, py);
      ctx.stroke();
    }

    // Draw each series
    for (const s of seriesList) {
      if (!s.data || s.data.length === 0) continue;

      // Filled area
      ctx.beginPath();
      ctx.moveTo(getX(s.data[0].x), getY(0));
      for (const pt of s.data) {
        ctx.lineTo(getX(pt.x), getY(pt.y));
      }
      ctx.lineTo(getX(s.data[s.data.length - 1].x), getY(0));
      ctx.closePath();
      ctx.fillStyle = s.fillColor || "rgba(14, 165, 233, 0.15)";
      ctx.fill();

      // Stroke line
      ctx.beginPath();
      ctx.moveTo(getX(s.data[0].x), getY(s.data[0].y));
      for (let i = 1; i < s.data.length; i++) {
        ctx.lineTo(getX(s.data[i].x), getY(s.data[i].y));
      }
      ctx.strokeStyle = s.strokeColor || "#0ea5e9";
      ctx.lineWidth = 2.5;
      ctx.stroke();
    }

    // Axes
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
      const px = getX(val);
      ctx.fillText(this.formatNum(val, 2), px, padding.top + plotH + 18);
    }

    // Legend
    let legX = padding.left + 10;
    ctx.font = "600 12px system-ui, -apple-system, sans-serif";
    for (const s of seriesList) {
      ctx.fillStyle = s.strokeColor;
      ctx.fillRect(legX, padding.top - 20, 14, 10);
      ctx.fillStyle = textBright;
      ctx.textAlign = "left";
      ctx.fillText(s.name, legX + 20, padding.top - 11);
      legX += ctx.measureText(s.name).width + 45;
    }
  },

  /**
   * Draws Tri-Plot (Prior, Normalized Likelihood, and Posterior)
   */
  drawTriPlot(canvas, triPlotData, options = {}) {
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

    // Posterior Filled Area
    ctx.beginPath();
    ctx.moveTo(getX(triPlotData[0][xProp]), getY(0));
    for (const pt of triPlotData) {
      ctx.lineTo(getX(pt[xProp]), getY(pt.posterior));
    }
    ctx.lineTo(getX(triPlotData[triPlotData.length - 1][xProp]), getY(0));
    ctx.closePath();
    ctx.fillStyle = isDark ? "rgba(56, 189, 248, 0.20)" : "rgba(2, 132, 199, 0.15)";
    ctx.fill();

    // 1. Prior curve (dashed amber)
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

    // 2. Likelihood curve (dotted green)
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

    // 3. Posterior curve (solid blue)
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

  /**
   * Renders Interactive SVG Fagan Nomogram
   */
  renderFaganNomogram(container, preTestProb, plr, nlr, postProbPos, postProbNeg) {
    if (!container) return;
    const isDark = document.documentElement.getAttribute("data-theme") === "dark";
    const strokeColor = isDark ? "#334155" : "#cbd5e1";
    const textColor = isDark ? "#94a3b8" : "#475569";
    const textBright = isDark ? "#f8fafc" : "#0f172a";
    const linePosColor = isDark ? "#38bdf8" : "#0284c7";
    const lineNegColor = isDark ? "#f43f5e" : "#e11d48";

    const width = 640;
    const height = 460;
    const margin = { top: 40, bottom: 40, left: 70, right: 70 };
    const axisH = height - margin.top - margin.bottom;

    // Fagan scale mapping: log-odds y = log(p / (1-p)) for 0.001 to 0.999
    const minProb = 0.001;
    const maxProb = 0.999;
    const minOdds = Math.log(minProb / (1 - minProb));
    const maxOdds = Math.log(maxProb / (1 - maxProb));

    const probToY = p => {
      p = Math.max(minProb, Math.min(maxProb, p));
      const odds = Math.log(p / (1 - p));
      return margin.top + (1 - (odds - minOdds) / (maxOdds - minOdds)) * axisH;
    };

    // LR axis in center: log(LR) from -3 to +3 (0.001 to 1000)
    const minLR = -3;
    const maxLR = 3;
    const lrToY = lr => {
      const logLR = Math.log10(Math.max(0.001, Math.min(1000, lr)));
      return margin.top + (1 - (logLR - minLR) / (maxLR - minLR)) * axisH;
    };

    const xPre = margin.left;
    const xLR = width / 2;
    const xPost = width - margin.right;

    // Y coordinates of trajectories
    const yPre = probToY(preTestProb);
    const yPosLR = lrToY(plr);
    const yPostPos = probToY(postProbPos);

    const yNegLR = lrToY(nlr);
    const yPostNeg = probToY(postProbNeg);

    const probTicks = [0.001, 0.005, 0.01, 0.02, 0.05, 0.10, 0.20, 0.30, 0.50, 0.70, 0.80, 0.90, 0.95, 0.98, 0.99, 0.999];
    const lrTicks = [1000, 500, 200, 100, 50, 20, 10, 5, 2, 1, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005, 0.002, 0.001];

    let svg = `<svg viewBox="0 0 ${width} ${height}" class="fagan-nomogram-svg" style="width: 100%; height: auto; font-family: system-ui, -apple-system, sans-serif;">
      <!-- Title -->
      <text x="${width/2}" y="22" text-anchor="middle" font-size="14" font-weight="700" fill="${textBright}">Bayesian Clinical Fagan Nomogram</text>

      <!-- Axes Vertical Lines -->
      <line x1="${xPre}" y1="${margin.top}" x2="${xPre}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />
      <line x1="${xLR}" y1="${margin.top}" x2="${xLR}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />
      <line x1="${xPost}" y1="${margin.top}" x2="${xPost}" y2="${margin.top + axisH}" stroke="${strokeColor}" stroke-width="2" />

      <!-- Axis Headers -->
      <text x="${xPre}" y="${margin.top - 12}" text-anchor="middle" font-size="12" font-weight="600" fill="${textBright}">Pre-Test Prob</text>
      <text x="${xLR}" y="${margin.top - 12}" text-anchor="middle" font-size="12" font-weight="600" fill="${textBright}">Likelihood Ratio</text>
      <text x="${xPost}" y="${margin.top - 12}" text-anchor="middle" font-size="12" font-weight="600" fill="${textBright}">Post-Test Prob</text>
    `;

    // Draw Pre & Post probability tick marks
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

    // Draw LR tick marks
    for (const lr of lrTicks) {
      const y = lrToY(lr);
      svg += `
        <line x1="${xLR - 5}" y1="${y}" x2="${xLR + 5}" y2="${y}" stroke="${strokeColor}" stroke-width="1.5" />
        <text x="${xLR + 8}" y="${y + 3}" text-anchor="start" font-size="10" fill="${textColor}">${lr}</text>
      `;
    }

    // Trajectory for Positive Test (Pre-test -> LR+ -> Post-test+)
    svg += `
      <!-- Positive test line -->
      <line x1="${xPre}" y1="${yPre}" x2="${xPost}" y2="${yPostPos}" stroke="${linePosColor}" stroke-width="3" stroke-linecap="round" opacity="0.9" />
      <circle cx="${xPre}" cy="${yPre}" r="5" fill="${linePosColor}" />
      <circle cx="${xLR}" cy="${yPosLR}" r="5" fill="${linePosColor}" />
      <circle cx="${xPost}" cy="${yPostPos}" r="6" fill="${linePosColor}" />
      <text x="${xPost + 42}" y="${yPostPos + 4}" font-size="11" font-weight="700" fill="${linePosColor}">+Test: ${(postProbPos * 100).toFixed(1)}%</text>

      <!-- Negative test line -->
      <line x1="${xPre}" y1="${yPre}" x2="${xPost}" y2="${yPostNeg}" stroke="${lineNegColor}" stroke-width="2.5" stroke-dasharray="6,4" stroke-linecap="round" opacity="0.9" />
      <circle cx="${xLR}" cy="${yNegLR}" r="4" fill="${lineNegColor}" />
      <circle cx="${xPost}" cy="${yPostNeg}" r="5" fill="${lineNegColor}" />
      <text x="${xPost + 42}" y="${yPostNeg + 4}" font-size="11" font-weight="700" fill="${lineNegColor}">−Test: ${(postProbNeg * 100).toFixed(1)}%</text>
    </svg>`;

    container.innerHTML = svg;
  },

  /**
   * Exports canvas as high-resolution PNG
   */
  exportPNG(canvas, filename = "bayes_chart.png") {
    if (!canvas) return;
    const link = document.createElement("a");
    link.download = filename;
    link.href = canvas.toDataURL("image/png");
    link.click();
  }
};
