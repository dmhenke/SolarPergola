/* =====================================================================
 * Solar Pergola — Annual PV Generation Model (Salt Lake City, 84108)
 * ---------------------------------------------------------------------
 * Engine: solar position (Spencer 1971 / NOAA), sunrise/sunset, clear-sky
 * shape (Haurwitz) scaled to NREL long-term monthly GHI and DHI for
 * Salt Lake City (WBAN 24127), Hay–Davies transposition onto a south-
 * facing (azimuth 180°) monoslope roof, Martin–Ruiz incidence-angle
 * modifier, NOCT cell temperature, linear power temperature coefficient,
 * PVWatts-style lumped system losses and inverter efficiency.
 *
 * The engine (SolarEngine) has no DOM dependencies, so it can be tested
 * with Node: `node -e "require('./solar_yield.js')"`.
 * ===================================================================== */
(function (root) {
  "use strict";

  // ------------------------------------------------------------------
  // Site & climate inputs
  // ------------------------------------------------------------------
  var SITE = { lat: 40.755, lon: -111.83, tz: -7, elevFt: 5500 };

  // NREL Solar Radiation Data Manual (1961–1990), Salt Lake City WBAN 24127,
  // as tabulated by NOAA NCEI Engineering Weather Data (Btu/ft²/day).
  var GHI_BTU = [610, 920, 1300, 1710, 2060, 2330, 2320, 2060, 1660, 1170, 710, 530];
  var DHI_BTU = [350, 440, 560, 680, 720, 670, 610, 550, 450, 360, 330, 300];
  var BTU2KWH = 0.0031546; // 1 Btu/ft² = 3.1546 Wh/m²

  // NOAA 1991–2020 normals, Salt Lake City Intl (°F)
  var TMAX_F = [38.6, 44.7, 55.3, 61.9, 72.6, 84.1, 94.0, 91.7, 80.6, 65.5, 50.7, 39.0];
  var TMIN_F = [24.2, 28.6, 36.3, 41.8, 50.4, 59.1, 68.2, 66.6, 56.3, 43.6, 32.8, 25.3];

  // Ground albedo (higher with intermittent winter snow cover)
  var ALBEDO = [0.35, 0.35, 0.25, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.25, 0.35];

  var MDAYS = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];
  var MNAME = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"];
  var ND = 365;
  var STEP = 10; // minutes
  var NS = (24 * 60) / STEP;
  var DEG = Math.PI / 180;

  // 2026 US DST: Sun Mar 8 (day index 66) through Sat Oct 31 (index 303)
  var DST_START = 66, DST_END = 304;

  var MSTART = [], MMID = [], DMONTH = new Int8Array(ND);
  (function () {
    var s = 0;
    for (var m = 0; m < 12; m++) {
      MSTART.push(s);
      MMID.push(s + MDAYS[m] / 2 - 0.5);
      for (var d = 0; d < MDAYS[m]; d++) DMONTH[s + d] = m;
      s += MDAYS[m];
    }
  })();

  // Smooth (circular, piecewise-linear between mid-month points) daily series
  function interpMonthly(arr) {
    var out = new Float64Array(ND);
    for (var d = 0; d < ND; d++) {
      var m = 11, a, b;
      for (var k = 11; k >= 0; k--) { if (MMID[k] <= d) { m = k; break; } }
      if (d < MMID[0]) { a = MMID[11] - ND; b = MMID[0]; m = 11; }
      else { a = MMID[m]; b = m === 11 ? MMID[0] + ND : MMID[m + 1]; }
      var t = (d - a) / (b - a);
      out[d] = arr[m] * (1 - t) + arr[(m + 1) % 12] * t;
    }
    return out;
  }

  function fToC(f) { return (f - 32) * 5 / 9; }

  // Spencer (1971) declination [rad] and equation of time [min]
  function solarTerms(dayIdx, hourStd) {
    var g = 2 * Math.PI / 365 * (dayIdx + (hourStd - 12) / 24);
    var decl = 0.006918 - 0.399912 * Math.cos(g) + 0.070257 * Math.sin(g)
      - 0.006758 * Math.cos(2 * g) + 0.000907 * Math.sin(2 * g)
      - 0.002697 * Math.cos(3 * g) + 0.00148 * Math.sin(3 * g);
    var eot = 229.18 * (0.000075 + 0.001868 * Math.cos(g) - 0.032077 * Math.sin(g)
      - 0.014615 * Math.cos(2 * g) - 0.040849 * Math.sin(2 * g));
    var e0 = 1.00011 + 0.034221 * Math.cos(g) + 0.00128 * Math.sin(g)
      + 0.000719 * Math.cos(2 * g) + 0.000077 * Math.sin(2 * g);
    return { decl: decl, eot: eot, i0: 1361 * e0 };
  }

  // Solar position for a day index and standard-time hour
  function sunPos(dayIdx, hourStd) {
    var st = solarTerms(dayIdx, hourStd);
    var phi = SITE.lat * DEG;
    var tst = hourStd * 60 + st.eot + 4 * SITE.lon - 60 * SITE.tz; // true solar time, min
    var w = (tst / 4 - 180) * DEG; // hour angle
    var cosZ = Math.sin(phi) * Math.sin(st.decl) + Math.cos(phi) * Math.cos(st.decl) * Math.cos(w);
    cosZ = Math.max(-1, Math.min(1, cosZ));
    var az = Math.atan2(Math.sin(w), Math.cos(w) * Math.sin(phi) - Math.tan(st.decl) * Math.cos(phi)) + Math.PI;
    return {
      cosZ: cosZ, sinZ: Math.sqrt(1 - cosZ * cosZ), alt: 90 - Math.acos(cosZ) / DEG,
      az: az / DEG, i0: st.i0, solarHour: tst / 60, decl: st.decl, eot: st.eot
    };
  }

  // Sunrise / sunset / solar noon in standard time hours (apparent, -0.833°)
  function sunTimes(dayIdx) {
    var st = solarTerms(dayIdx, 12);
    var phi = SITE.lat * DEG;
    var c = (Math.sin(-0.833 * DEG) - Math.sin(phi) * Math.sin(st.decl)) / (Math.cos(phi) * Math.cos(st.decl));
    c = Math.max(-1, Math.min(1, c));
    var w0 = Math.acos(c) / DEG; // deg
    var corr = (st.eot + 4 * SITE.lon - 60 * SITE.tz) / 60; // h
    var noon = 12 - corr;
    return { rise: noon - w0 / 15, set: noon + w0 / 15, noon: noon, dayLen: 2 * w0 / 15, decl: st.decl / DEG };
  }

  function isDST(d) { return d >= DST_START && d < DST_END; }

  // ------------------------------------------------------------------
  // Precompute the sky & sun on a 365 × 144 grid (independent of pitch)
  // ------------------------------------------------------------------
  var N = ND * NS;
  var gCosZ = new Float32Array(N), gSinZ = new Float32Array(N), gCosAzS = new Float32Array(N);
  var gAlt = new Float32Array(N), gAz = new Float32Array(N);
  var gGHI = new Float32Array(N), gDHI = new Float32Array(N), gBh = new Float32Array(N);
  var gAi = new Float32Array(N), gTa = new Float32Array(N);
  var dayGHI = new Float64Array(ND), dayDHI = new Float64Array(ND), dayAlb = interpMonthly(ALBEDO);
  var TIMES = [];

  (function precompute() {
    var cs = new Float32Array(N), csDay = new Float64Array(ND);
    var i, d, s;
    for (d = 0; d < ND; d++) {
      TIMES.push(sunTimes(d));
      for (s = 0; s < NS; s++) {
        i = d * NS + s;
        var p = sunPos(d, (s + 0.5) * STEP / 60);
        gCosZ[i] = p.cosZ; gSinZ[i] = p.sinZ; gAlt[i] = p.alt; gAz[i] = p.az;
        gCosAzS[i] = Math.cos((p.az - 180) * DEG);
        gAi[i] = p.i0; // temporarily hold extraterrestrial normal irradiance
        gTa[i] = p.solarHour;
        if (p.cosZ > 0) {
          // Haurwitz clear-sky GHI (W/m²) — used only for intra-day/inter-day shape
          cs[i] = 1098 * p.cosZ * Math.exp(-0.057 / p.cosZ);
          csDay[d] += cs[i] * STEP / 60 / 1000; // kWh/m²
        }
      }
    }
    // Daily clear-sky index, iterated so each month reproduces the NREL mean exactly
    var target = GHI_BTU.map(function (v) { return v * BTU2KWH; });
    var kM = [];
    for (var m = 0; m < 12; m++) {
      var sum = 0; for (d = MSTART[m]; d < MSTART[m] + MDAYS[m]; d++) sum += csDay[d];
      kM.push(target[m] / (sum / MDAYS[m]));
    }
    var kd = interpMonthly(kM);
    for (var it = 0; it < 8; it++) {
      var ratio = [];
      for (m = 0; m < 12; m++) {
        var acc = 0; for (d = MSTART[m]; d < MSTART[m] + MDAYS[m]; d++) acc += kd[d] * csDay[d];
        ratio.push(target[m] / (acc / MDAYS[m]));
      }
      var rd = interpMonthly(ratio);
      for (d = 0; d < ND; d++) kd[d] *= rd[d];
    }
    var fd = interpMonthly(DHI_BTU.map(function (v, k) { return v / GHI_BTU[k]; }));
    var tMean = interpMonthly(TMAX_F.map(function (v, k) { return fToC((v + TMIN_F[k]) / 2); }));
    var tAmp = interpMonthly(TMAX_F.map(function (v, k) { return (v - TMIN_F[k]) * 5 / 9 / 2; }));

    for (d = 0; d < ND; d++) {
      for (s = 0; s < NS; s++) {
        i = d * NS + s;
        var solarHour = gTa[i], i0 = gAi[i];
        // Diurnal temperature: min near dawn, max ~15:00 solar time
        gTa[i] = tMean[d] + tAmp[d] * Math.cos(2 * Math.PI * (solarHour - 15) / 24);
        if (gCosZ[i] <= 0) { gAi[i] = 0; continue; }
        var G = kd[d] * cs[i];
        var D = fd[d] * G;
        var Bh = G - D;
        var dni = Bh / Math.max(gCosZ[i], 0.0872);
        gGHI[i] = G; gDHI[i] = D; gBh[i] = Bh;
        gAi[i] = Math.min(1, dni / i0); // Hay–Davies anisotropy index
        dayGHI[d] += G * STEP / 60 / 1000;
        dayDHI[d] += D * STEP / 60 / 1000;
      }
    }
  })();

  // Martin & Ruiz (2001) incidence angle modifier, a_r = 0.16
  var IAM_DEN = 1 - Math.exp(-1 / 0.16);
  function iam(cosT) { return cosT <= 0 ? 0 : (1 - Math.exp(-cosT / 0.16)) / IAM_DEN; }

  // ------------------------------------------------------------------
  // Simulate one pitch
  //   opts: { pdcW, lossFrac, gamma, invEff, horizonE, snowFrac, keepGrid }
  // ------------------------------------------------------------------
  function simulate(pitchDeg, opts) {
    var b = pitchDeg * DEG, cb = Math.cos(b), sb = Math.sin(b);
    var iso = (1 + cb) / 2, gnd = (1 - cb) / 2;
    var pdc = opts.pdcW, gamma = opts.gamma, derate = (1 - opts.lossFrac) * opts.invEff;
    var hE = opts.horizonE || 0, snow = opts.snowFrac || 0;
    var keep = !!opts.keepGrid;
    var acGrid = keep ? new Float32Array(N) : null;
    var poaGrid = keep ? new Float32Array(N) : null;
    var cosTGrid = keep ? new Float32Array(N) : null;
    var dayKWh = new Float64Array(ND), dayPOA = new Float64Array(ND);
    var dt = STEP / 60;
    for (var d = 0; d < ND; d++) {
      var m = DMONTH[d];
      var snowF = (m === 11 || m === 0 || m === 1) ? 1 - snow : 1;
      var alb = dayAlb[d];
      var eDay = 0, pDay = 0;
      for (var s = 0; s < NS; s++) {
        var i = d * NS + s;
        var cosZ = gCosZ[i];
        if (cosZ <= 0) continue;
        var cosT = cosZ * cb + gSinZ[i] * sb * gCosAzS[i];
        if (cosT < 0) cosT = 0;
        var shaded = hE > 0 && gAz[i] > 45 && gAz[i] < 150 && gAlt[i] < hE;
        var rb = cosT / Math.max(cosZ, 0.0872);
        var G = gGHI[i], D = gDHI[i], Ai = gAi[i];
        var beam = shaded ? 0 : gBh[i] * rb;
        var circ = shaded ? 0 : D * Ai * rb;
        var isoD = D * (1 - Ai) * iso;
        var grd = G * alb * gnd;
        var poa = beam + circ + isoD + grd;
        var k = iam(cosT);
        var eff = (beam + circ) * k + (isoD + grd) * 0.95;
        var tc = gTa[i] + poa * 25 / 800; // NOCT 45 °C, open-rack
        var ac = pdc * eff / 1000 * (1 + gamma * (tc - 25)) * derate * snowF;
        if (ac < 0) ac = 0;
        eDay += ac * dt; pDay += poa * dt;
        if (keep) { acGrid[i] = ac; poaGrid[i] = poa; cosTGrid[i] = cosT; }
      }
      dayKWh[d] = eDay / 1000; dayPOA[d] = pDay / 1000;
    }
    var month = [], monthPOA = [], annual = 0, annualPOA = 0;
    for (var mm = 0; mm < 12; mm++) {
      var e = 0, p = 0;
      for (var dd = MSTART[mm]; dd < MSTART[mm] + MDAYS[mm]; dd++) { e += dayKWh[dd]; p += dayPOA[dd]; }
      month.push(e); monthPOA.push(p); annual += e; annualPOA += p;
    }
    return {
      pitch: pitchDeg, dayKWh: dayKWh, dayPOA: dayPOA, month: month, monthPOA: monthPOA,
      annual: annual, annualPOA: annualPOA, acGrid: acGrid, poaGrid: poaGrid, cosTGrid: cosTGrid
    };
  }

  var engine = {
    SITE: SITE, ND: ND, NS: NS, STEP: STEP, MDAYS: MDAYS, MSTART: MSTART, MNAME: MNAME, DMONTH: DMONTH,
    GHI_KWH: GHI_BTU.map(function (v) { return v * BTU2KWH; }),
    DHI_KWH: DHI_BTU.map(function (v) { return v * BTU2KWH; }),
    TIMES: TIMES, dayGHI: dayGHI, dayDHI: dayDHI,
    gAlt: gAlt, gAz: gAz, gTa: gTa, gGHI: gGHI,
    isDST: isDST, sunPos: sunPos, sunTimes: sunTimes, simulate: simulate,
    // Combined relative uncertainty (approx. 95%): NREL data ±9%,
    // transposition/temperature model ±4%, loss assumptions ±3% (root-sum-square)
    REL_UNC: Math.sqrt(0.09 * 0.09 + 0.04 * 0.04 + 0.03 * 0.03)
  };

  if (typeof module !== "undefined" && module.exports) module.exports = engine;
  root.SolarEngine = engine;
})(typeof window !== "undefined" ? window : this);

/* =====================================================================
 * UI layer — only runs in the browser when the #solar tab is present
 * ===================================================================== */
(function () {
  "use strict";
  if (typeof document === "undefined") return;
  var E = window.SolarEngine;
  var MOD_LEN_FT = 93.8 / 12, MOD_WID_FT = 44.6 / 12;

  function $(id) { return document.getElementById(id); }
  function cssVar(n) { return getComputedStyle(document.documentElement).getPropertyValue(n).trim(); }
  function isDark() { return document.documentElement.getAttribute("data-theme") === "dark"; }
  function fmt(n, dp) { return n.toLocaleString(undefined, { minimumFractionDigits: dp || 0, maximumFractionDigits: dp || 0 }); }
  function hhmm(h) {
    h = ((h % 24) + 24) % 24;
    var hr = Math.floor(h), mn = Math.round((h - hr) * 60);
    if (mn === 60) { hr += 1; mn = 0; }
    var ap = hr >= 12 ? "pm" : "am", h12 = hr % 12 === 0 ? 12 : hr % 12;
    return h12 + ":" + (mn < 10 ? "0" : "") + mn + " " + ap;
  }
  function dateLabel(d) {
    var m = E.DMONTH[d];
    return E.MNAME[m] + " " + (d - E.MSTART[m] + 1);
  }
  function palette() {
    var dark = isDark();
    return {
      text: cssVar("--global-text-color"), muted: cssVar("--global-text-color-light"),
      grid: cssVar("--global-divider-color"), surface: cssVar("--global-surface-color"),
      bg: cssVar("--global-bg-color"),
      s1: cssVar("--global-theme-color"),          // selected pitch
      s2: dark ? "#f59e0b" : "#b45309",              // optimum pitch (amber)
      s3: dark ? "#9ca3af" : "#6b7280",              // flat reference (neutral)
      band: dark ? "rgba(45,212,191,0.16)" : "rgba(13,148,136,0.12)",
      zone: dark ? "rgba(16,185,129,0.12)" : "rgba(5,150,105,0.09)",
      // Sequential amber ramp for power (light → dark in light mode; inverse in dark)
      ramp: dark ? [[42, 31, 14], [120, 72, 14], [217, 119, 6], [251, 191, 36], [254, 243, 199]]
                 : [[255, 247, 230], [253, 211, 120], [245, 158, 11], [180, 83, 9], [120, 45, 8]]
    };
  }
  function rampColor(ramp, t) {
    t = Math.max(0, Math.min(1, t));
    var x = t * (ramp.length - 1), i = Math.min(ramp.length - 2, Math.floor(x)), f = x - i;
    var a = ramp[i], b = ramp[i + 1];
    return [Math.round(a[0] + (b[0] - a[0]) * f), Math.round(a[1] + (b[1] - a[1]) * f), Math.round(a[2] + (b[2] - a[2]) * f)];
  }

  // ---------------- SVG helpers ----------------
  var NS_SVG = "http://www.w3.org/2000/svg";
  function el(tag, attrs, parent) {
    var e = document.createElementNS(NS_SVG, tag);
    for (var k in attrs) e.setAttribute(k, attrs[k]);
    if (parent) parent.appendChild(e);
    return e;
  }
  function txt(parent, x, y, s, attrs) {
    var a = { x: x, y: y, "font-size": 11, "font-family": "Inter, sans-serif" };
    for (var k in attrs) a[k] = attrs[k];
    var t = el("text", a, parent); t.textContent = s; return t;
  }
  function niceTicks(lo, hi, n) {
    var span = hi - lo, step = Math.pow(10, Math.floor(Math.log10(span / n)));
    var err = span / n / step;
    if (err >= 7.5) step *= 10; else if (err >= 3.5) step *= 5; else if (err >= 1.5) step *= 2;
    var out = [];
    for (var v = Math.ceil(lo / step) * step; v <= hi + 1e-9; v += step) out.push(+v.toFixed(10));
    return out;
  }
  function newSvg(host, h) {
    host.innerHTML = "";
    var w = Math.max(300, host.clientWidth);
    var svg = el("svg", { viewBox: "0 0 " + w + " " + h, width: w, height: h, role: "img" }, host);
    return { svg: svg, w: w, h: h };
  }

  // ---------------- Tooltip ----------------
  var tip;
  function showTip(html, cx, cy) {
    if (!tip) { tip = document.createElement("div"); tip.className = "sy-tip"; document.body.appendChild(tip); }
    tip.innerHTML = html; tip.style.display = "block";
    var r = tip.getBoundingClientRect();
    var x = cx + 14, y = cy + 14;
    if (x + r.width > window.innerWidth - 8) x = cx - r.width - 14;
    if (y + r.height > window.innerHeight - 8) y = cy - r.height - 14;
    tip.style.left = Math.max(8, x) + "px"; tip.style.top = Math.max(8, y) + "px";
  }
  function hideTip() { if (tip) tip.style.display = "none"; }

  // ---------------- State ----------------
  var state = { res: null, flat: null, best: null, sweep: null, opts: null, nMod: 0, pitch: 7 };

  function readInputs() {
    var pitch = parseFloat($("sy-pitch").value);
    var W = parseFloat(($("calc-depth") || { value: 16 }).value);
    var L = parseFloat(($("calc-span") || { value: 20 }).value);
    var ovFt = parseFloat(($("calc-mod-overhang") || { value: 6 }).value) / 12;
    var auto = $("sy-auto").checked;
    var rows = Math.floor((W / Math.cos(pitch * Math.PI / 180) + 2 * ovFt) / MOD_LEN_FT + 1e-6);
    var cols = Math.floor((L + 2 * ovFt) / MOD_WID_FT + 1e-6);
    var nAuto = Math.max(1, rows * cols);
    var nMod = auto ? nAuto : Math.max(1, parseInt($("sy-nmod").value, 10) || 1);
    if (auto) $("sy-nmod").value = nAuto;
    $("sy-nmod").disabled = auto;
    var watt = parseFloat($("sy-watt").value);
    return {
      pitch: pitch, nMod: nMod, rows: rows, cols: cols, W: W, L: L, watt: watt,
      opts: {
        pdcW: nMod * watt, lossFrac: parseFloat($("sy-loss").value) / 100, gamma: -0.0030, invEff: 0.96,
        horizonE: parseFloat($("sy-horizon").value), snowFrac: parseFloat($("sy-snow").value) / 100
      }
    };
  }

  function compute() {
    var inp = readInputs();
    $("disp-sy-pitch").textContent = inp.pitch + "°";
    $("disp-sy-watt").textContent = inp.watt + " W";
    $("disp-sy-loss").textContent = $("sy-loss").value + "%";
    $("disp-sy-horizon").textContent = $("sy-horizon").value + "°";
    $("disp-sy-snow").textContent = $("sy-snow").value + "%";
    $("sy-array-note").textContent = $("sy-auto").checked
      ? inp.rows + " rows × " + inp.cols + " columns of 93.8 × 44.6 in modules (portrait) on the " +
        fmt(inp.W, 1) + " × " + fmt(inp.L, 1) + " ft designer footprint"
      : "Manual module count";
    var o = inp.opts;
    state.opts = o; state.nMod = inp.nMod; state.pitch = inp.pitch;
    var keep = {}; for (var k in o) keep[k] = o[k]; keep.keepGrid = true;
    state.res = E.simulate(inp.pitch, keep);
    state.flat = E.simulate(0, o);
    state.sweep = [];
    var best = null;
    for (var p = 0; p <= 45; p++) {
      var r = E.simulate(p, o);
      state.sweep.push({ p: p, e: r.annual });
      if (!best || r.annual > best.e) best = { p: p, e: r.annual };
    }
    state.best = E.simulate(best.p, o);
    renderStats(inp);
    renderAll();
  }

  function renderStats(inp) {
    var r = state.res, u = E.REL_UNC, kwp = inp.opts.pdcW / 1000;
    $("sy-annual").textContent = fmt(r.annual) + " kWh";
    $("sy-annual-sub").textContent = "95% UI " + fmt(r.annual * (1 - u)) + "–" + fmt(r.annual * (1 + u)) +
      " kWh · " + inp.nMod + " × " + inp.watt + " W = " + fmt(kwp, 2) + " kWp";
    $("sy-spec").textContent = fmt(r.annual / kwp) + " kWh/kWp";
    $("sy-spec-sub").textContent = "Capacity factor " + fmt(100 * r.annual / (kwp * 8760), 1) + "% (AC energy / DC nameplate)";
    var loss = 100 * (1 - r.annual / state.best.annual);
    $("sy-tilt").textContent = (loss < 0.05 ? "0.0" : "−" + fmt(loss, 1)) + "%";
    $("sy-tilt-sub").textContent = "vs. irradiance-optimal " + state.best.pitch + "° (" + fmt(state.best.annual) +
      " kWh); flat 0° roof yields " + fmt(state.flat.annual) + " kWh";
    var jun = r.month[5] / 30, dec = r.month[11] / 31;
    $("sy-season").textContent = fmt(jun, 1) + " / " + fmt(dec, 1);
    $("sy-season-sub").textContent = "kWh/day, June vs. December (" + fmt(jun / dec, 1) + "× seasonal swing)";
    var tJ = E.TIMES[171], tD = E.TIMES[354];
    $("sy-day").textContent = fmt(tJ.dayLen, 1) + " / " + fmt(tD.dayLen, 1) + " h";
    $("sy-day-sub").textContent = "Daylength Jun 21 vs. Dec 21; noon sun " + fmt(90 - 40.755 + 23.44, 1) + "° vs. " +
      fmt(90 - 40.755 - 23.44, 1) + "° above horizon";
  }

  // ---------------- Figure 1: Heatmap ----------------
  var HM = { y0: 4.5, y1: 21.5 }; // clock-hour window
  function renderHeatmap() {
    var host = $("sy-heatmap"), P = palette(), r = state.res;
    host.innerHTML = "";
    var w = Math.max(320, host.clientWidth), h = w < 560 ? 300 : 360;
    var m = { l: 52, r: 14, t: 12, b: 34 };
    var pw = w - m.l - m.r, ph = h - m.t - m.b;
    var dpr = window.devicePixelRatio || 1;
    var cv = document.createElement("canvas");
    cv.width = w * dpr; cv.height = h * dpr; cv.style.width = w + "px"; cv.style.height = h + "px";
    host.appendChild(cv);
    var ctx = cv.getContext("2d"); ctx.scale(dpr, dpr);
    var maxAC = 0; for (var i = 0; i < r.acGrid.length; i++) if (r.acGrid[i] > maxAC) maxAC = r.acGrid[i];
    var hrs = HM.y1 - HM.y0, stepH = E.STEP / 60;
    var cw = pw / E.ND, ch = ph / (hrs / stepH);
    ctx.fillStyle = P.surface; ctx.fillRect(m.l, m.t, pw, ph);
    ctx.strokeStyle = P.grid; ctx.lineWidth = 1;
    for (var gh = 6; gh <= 21; gh += 3) {
      var gy = m.t + (gh - HM.y0) / hrs * ph;
      ctx.beginPath(); ctx.moveTo(m.l, gy); ctx.lineTo(m.l + pw, gy); ctx.stroke();
    }
    for (var d = 0; d < E.ND; d++) {
      var shift = E.isDST(d) ? 1 : 0;
      for (var s = 0; s < E.NS; s++) {
        var v = r.acGrid[d * E.NS + s];
        if (v <= 0) continue;
        var clock = s * stepH + shift;
        if (clock < HM.y0 || clock >= HM.y1) continue;
        var c = rampColor(P.ramp, Math.pow(v / maxAC, 0.85));
        ctx.fillStyle = "rgb(" + c[0] + "," + c[1] + "," + c[2] + ")";
        ctx.fillRect(m.l + d * cw, m.t + (clock - HM.y0) / stepH * ch, cw + 0.6, ch + 0.6);
      }
    }
    function yOf(clock) { return m.t + (clock - HM.y0) / hrs * ph; }
    // Sunrise, sunset, solar noon
    [["rise", [] ], ["set", []], ["noon", []]].forEach(function (pair) {
      ctx.beginPath();
      ctx.strokeStyle = pair[0] === "noon" ? P.muted : P.text;
      ctx.lineWidth = pair[0] === "noon" ? 1 : 1.5;
      ctx.setLineDash(pair[0] === "noon" ? [4, 4] : []);
      var prevDst = null;
      for (var d = 0; d < E.ND; d++) {
        var dst = E.isDST(d), t = E.TIMES[d][pair[0]] + (dst ? 1 : 0);
        var x = m.l + (d + 0.5) * cw, y = yOf(t);
        if (d === 0 || dst !== prevDst) ctx.moveTo(x, y); else ctx.lineTo(x, y);
        prevDst = dst;
      }
      ctx.stroke();
    });
    ctx.setLineDash([]);
    // Axes
    ctx.fillStyle = P.muted; ctx.font = "11px Inter, sans-serif"; ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (var hh = 6; hh <= 21; hh += 3) {
      var yy = yOf(hh);
      ctx.fillText(hhmm(hh).replace(":00", ""), m.l - 6, yy);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (var mo = 0; mo < 12; mo++) {
      if (w < 520 && mo % 2) continue;
      ctx.fillText(E.MNAME[mo], m.l + (E.MSTART[mo] + E.MDAYS[mo] / 2) * cw, m.t + ph + 8);
    }
    // Direct labels on the overlay curves
    ctx.fillStyle = P.text; ctx.textAlign = "left"; ctx.font = "600 11px Inter, sans-serif";
    ctx.textBaseline = "bottom"; ctx.fillText("sunset", m.l + 140 * cw, yOf(E.TIMES[140].set + 1) - 3);
    ctx.textBaseline = "top"; ctx.fillText("sunrise", m.l + 140 * cw, yOf(E.TIMES[140].rise + 1) + 3);
    ctx.fillStyle = P.muted; ctx.textBaseline = "bottom"; ctx.font = "11px Inter, sans-serif";
    ctx.fillText("solar noon", m.l + 4, yOf(E.TIMES[2].noon) - 3);
    // Legend
    var lg = $("sy-heat-legend");
    var stops = [0, 0.25, 0.5, 0.75, 1].map(function (t) { var c = rampColor(P.ramp, t); return "rgb(" + c.join(",") + ") " + (t * 100) + "%"; });
    lg.innerHTML = '<span>0 kW</span><span class="sy-ramp" style="background:linear-gradient(90deg,' + stops.join(",") +
      ')"></span><span>' + fmt(maxAC / 1000, 1) + " kW AC</span>";
    // Hover
    cv.addEventListener("mousemove", function (ev) {
      var b = cv.getBoundingClientRect(), x = ev.clientX - b.left, y = ev.clientY - b.top;
      if (x < m.l || x > m.l + pw || y < m.t || y > m.t + ph) { hideTip(); return; }
      var d = Math.min(E.ND - 1, Math.floor((x - m.l) / cw));
      var clock = HM.y0 + (y - m.t) / ph * hrs;
      var s = Math.floor((clock - (E.isDST(d) ? 1 : 0)) / stepH);
      if (s < 0 || s >= E.NS) { hideTip(); return; }
      var i = d * E.NS + s, T = E.TIMES[d], dst = E.isDST(d) ? 1 : 0;
      var up = E.gAlt[i] > 0;
      var html = "<strong>" + dateLabel(d) + ", " + hhmm(s * stepH + dst) + (dst ? " MDT" : " MST") + "</strong><br>" +
        (up ? "Sun altitude " + fmt(E.gAlt[i], 1) + "°, azimuth " + fmt(E.gAz[i], 0) + "°<br>" +
          "Incidence on roof " + fmt(Math.acos(Math.min(1, r.cosTGrid[i])) * 180 / Math.PI, 1) + "°<br>" +
          "Plane-of-array " + fmt(r.poaGrid[i]) + " W/m²<br>Output <strong>" + fmt(r.acGrid[i] / 1000, 2) + " kW</strong>"
          : "Sun below horizon") +
        "<br><span class='sy-tip-muted'>Sunrise " + hhmm(T.rise + dst) + " · sunset " + hhmm(T.set + dst) +
        " · " + fmt(T.dayLen, 1) + " h</span>";
      showTip(html, ev.clientX, ev.clientY);
    });
    cv.addEventListener("mouseleave", hideTip);
  }

  // ---------------- Generic line chart (single y-axis) ----------------
  function lineChart(host, cfg) {
    var P = palette(), o = newSvg(host, cfg.h || 300), svg = o.svg, w = o.w, h = o.h;
    var m = { l: 52, r: 16, t: 14, b: 34 }, pw = w - m.l - m.r, ph = h - m.t - m.b;
    var X = function (v) { return m.l + (v - cfg.x0) / (cfg.x1 - cfg.x0) * pw; };
    var Y = function (v) { return m.t + ph - (v - cfg.y0) / (cfg.y1 - cfg.y0) * ph; };
    if (cfg.zones) cfg.zones.forEach(function (z) {
      el("rect", { x: X(z.a), y: m.t, width: X(z.b) - X(z.a), height: ph, fill: P.zone }, svg);
      txt(svg, (X(z.a) + X(z.b)) / 2, m.t + 12, z.label, { fill: P.muted, "text-anchor": "middle", "font-size": 10 });
    });
    niceTicks(cfg.y0, cfg.y1, 5).forEach(function (t) {
      el("line", { x1: m.l, x2: m.l + pw, y1: Y(t), y2: Y(t), stroke: P.grid, "stroke-width": 1 }, svg);
      txt(svg, m.l - 6, Y(t) + 4, cfg.yFmt ? cfg.yFmt(t) : fmt(t), { fill: P.muted, "text-anchor": "end" });
    });
    cfg.xTicks.forEach(function (t) { txt(svg, X(t.v), h - 12, t.s, { fill: P.muted, "text-anchor": "middle" }); });
    if (cfg.yLabel) txt(svg, 12, m.t + ph / 2, cfg.yLabel, { fill: P.muted, "text-anchor": "middle", transform: "rotate(-90 12 " + (m.t + ph / 2) + ")" });
    if (cfg.band) {
      var d = "";
      cfg.band.hi.forEach(function (v, k) { d += (k ? "L" : "M") + X(cfg.xs[k]).toFixed(1) + " " + Y(v).toFixed(1); });
      for (var k = cfg.band.lo.length - 1; k >= 0; k--) d += "L" + X(cfg.xs[k]).toFixed(1) + " " + Y(cfg.band.lo[k]).toFixed(1);
      el("path", { d: d + "Z", fill: P.band, stroke: "none" }, svg);
    }
    cfg.series.forEach(function (s) {
      var d = "";
      s.ys.forEach(function (v, k) { d += (k ? "L" : "M") + X(cfg.xs[k]).toFixed(1) + " " + Y(v).toFixed(1); });
      el("path", { d: d, fill: "none", stroke: P[s.color], "stroke-width": 2, "stroke-dasharray": s.dash || "none", "stroke-linejoin": "round" }, svg);
    });
    if (cfg.marks) cfg.marks.forEach(function (mk) {
      el("circle", { cx: X(mk.x), cy: Y(mk.y), r: 5, fill: P[mk.color], stroke: P.bg, "stroke-width": 2 }, svg);
      txt(svg, X(mk.x) + (mk.anchor === "end" ? -9 : 9), Y(mk.y) + (mk.below ? 18 : -8), mk.label, { fill: P.text, "text-anchor": mk.anchor || "start", "font-weight": 600 });
    });
    // Crosshair + tooltip
    var cross = el("line", { y1: m.t, y2: m.t + ph, stroke: P.muted, "stroke-width": 1, visibility: "hidden" }, svg);
    var dots = cfg.series.map(function (s) { return el("circle", { r: 4, fill: P[s.color], stroke: P.bg, "stroke-width": 2, visibility: "hidden" }, svg); });
    var hit = el("rect", { x: m.l, y: m.t, width: pw, height: ph, fill: "transparent" }, svg);
    hit.addEventListener("mousemove", function (ev) {
      var b = svg.getBoundingClientRect(), x = (ev.clientX - b.left) * (w / b.width);
      var v = cfg.x0 + (x - m.l) / pw * (cfg.x1 - cfg.x0);
      var k = 0, best = Infinity;
      cfg.xs.forEach(function (xx, j) { var dd = Math.abs(xx - v); if (dd < best) { best = dd; k = j; } });
      cross.setAttribute("x1", X(cfg.xs[k])); cross.setAttribute("x2", X(cfg.xs[k])); cross.setAttribute("visibility", "visible");
      dots.forEach(function (dt, j) { dt.setAttribute("cx", X(cfg.xs[k])); dt.setAttribute("cy", Y(cfg.series[j].ys[k])); dt.setAttribute("visibility", "visible"); });
      showTip(cfg.tip(k), ev.clientX, ev.clientY);
    });
    hit.addEventListener("mouseleave", function () {
      cross.setAttribute("visibility", "hidden"); dots.forEach(function (dt) { dt.setAttribute("visibility", "hidden"); }); hideTip();
    });
  }

  function legend(id, items) {
    var P = palette();
    $(id).innerHTML = items.map(function (it) {
      var sw = it.band ? '<span class="sy-sw sy-sw-band" style="background:' + P.band + '"></span>'
        : '<span class="sy-sw" style="border-top:2px ' + (it.dash ? "dashed " : "solid ") + P[it.color] + '"></span>';
      return '<span class="sy-leg-item">' + sw + it.label + "</span>";
    }).join("");
  }

  // ---------------- Figure 2: Daily energy ----------------
  function renderDaily() {
    var r = state.res, f = state.flat, b = state.best, u = E.REL_UNC;
    var xs = [], ys = [], yf = [], yb = [], lo = [], hi = [];
    // 7-day centered running mean for readability (daily values are climatological, not weather)
    function smooth(a, d) { var s = 0, n = 0; for (var k = -3; k <= 3; k++) { var j = (d + k + E.ND) % E.ND; s += a[j]; n++; } return s / n; }
    for (var d = 0; d < E.ND; d++) {
      xs.push(d); var v = smooth(r.dayKWh, d);
      ys.push(v); lo.push(v * (1 - u)); hi.push(v * (1 + u));
      yf.push(smooth(f.dayKWh, d)); yb.push(smooth(b.dayKWh, d));
    }
    var ymax = Math.max.apply(null, hi.concat(yb)) * 1.05;
    var series = [{ ys: ys, color: "s1" }, { ys: yf, color: "s3", dash: "5 4" }];
    if (b.pitch !== r.pitch) series.push({ ys: yb, color: "s2" });
    lineChart($("sy-daily"), {
      h: 300, xs: xs, x0: 0, x1: E.ND - 1, y0: 0, y1: ymax, yLabel: "kWh / day",
      xTicks: E.MNAME.map(function (n, k) { return { v: E.MSTART[k] + E.MDAYS[k] / 2, s: n }; })
        .filter(function (t, k) { return $("sy-daily").clientWidth >= 520 || k % 2 === 0; }),
      series: series, band: { lo: lo, hi: hi },
      tip: function (k) {
        var t = E.TIMES[k];
        return "<strong>" + dateLabel(k) + "</strong><br>" + state.pitch + "° roof: <strong>" + fmt(ys[k], 1) +
          " kWh</strong> (" + fmt(lo[k], 1) + "–" + fmt(hi[k], 1) + ")<br>Flat 0°: " + fmt(yf[k], 1) + " kWh" +
          (b.pitch !== r.pitch ? "<br>Optimal " + b.pitch + "°: " + fmt(yb[k], 1) + " kWh" : "") +
          "<br><span class='sy-tip-muted'>Daylength " + fmt(t.dayLen, 1) + " h · noon altitude " +
          fmt(90 - E.SITE.lat + t.decl, 1) + "°</span>";
      }
    });
    var items = [{ color: "s1", label: state.pitch + "° pitch (this design)" }, { band: true, label: "95% uncertainty interval" },
      { color: "s3", dash: true, label: "Flat 0° reference" }];
    if (b.pitch !== r.pitch) items.push({ color: "s2", label: "Irradiance-optimal " + b.pitch + "°" });
    legend("sy-daily-legend", items);
  }

  // ---------------- Figure 3: Monthly bars ----------------
  function renderMonthly() {
    var P = palette(), r = state.res, u = E.REL_UNC;
    var o = newSvg($("sy-monthly"), 280), svg = o.svg, w = o.w, h = o.h;
    var m = { l: 52, r: 12, t: 14, b: 30 }, pw = w - m.l - m.r, ph = h - m.t - m.b;
    var ymax = Math.max.apply(null, r.month) * (1 + u) * 1.08;
    var Y = function (v) { return m.t + ph - v / ymax * ph; };
    niceTicks(0, ymax, 5).forEach(function (t) {
      el("line", { x1: m.l, x2: m.l + pw, y1: Y(t), y2: Y(t), stroke: P.grid }, svg);
      txt(svg, m.l - 6, Y(t) + 4, fmt(t), { fill: P.muted, "text-anchor": "end" });
    });
    txt(svg, 12, m.t + ph / 2, "kWh / month", { fill: P.muted, "text-anchor": "middle", transform: "rotate(-90 12 " + (m.t + ph / 2) + ")" });
    var slot = pw / 12, bw = Math.min(34, slot * 0.62);
    r.month.forEach(function (v, k) {
      var cx = m.l + slot * (k + 0.5), y = Y(v), rr = Math.min(4, bw / 2);
      el("path", { d: "M" + (cx - bw / 2) + " " + Y(0) + "V" + (y + rr) + "Q" + (cx - bw / 2) + " " + y + " " + (cx - bw / 2 + rr) + " " + y +
        "H" + (cx + bw / 2 - rr) + "Q" + (cx + bw / 2) + " " + y + " " + (cx + bw / 2) + " " + (y + rr) + "V" + Y(0) + "Z", fill: P.s1 }, svg);
      el("line", { x1: cx, x2: cx, y1: Y(v * (1 - u)), y2: Y(v * (1 + u)), stroke: P.text, "stroke-width": 1.25 }, svg);
      el("line", { x1: cx - 4, x2: cx + 4, y1: Y(v * (1 + u)), y2: Y(v * (1 + u)), stroke: P.text, "stroke-width": 1.25 }, svg);
      el("line", { x1: cx - 4, x2: cx + 4, y1: Y(v * (1 - u)), y2: Y(v * (1 - u)), stroke: P.text, "stroke-width": 1.25 }, svg);
      if (w >= 520 || k % 2 === 0) txt(svg, cx, h - 10, E.MNAME[k], { fill: P.muted, "text-anchor": "middle" });
      var hit = el("rect", { x: cx - slot / 2, y: m.t, width: slot, height: ph, fill: "transparent" }, svg);
      hit.addEventListener("mousemove", function (ev) {
        showTip("<strong>" + E.MNAME[k] + "</strong><br>" + fmt(v) + " kWh (95% UI " + fmt(v * (1 - u)) + "–" + fmt(v * (1 + u)) +
          ")<br>" + fmt(v / E.MDAYS[k], 1) + " kWh/day · " + fmt(100 * v / r.annual, 1) + "% of annual<br>" +
          "<span class='sy-tip-muted'>GHI " + fmt(E.GHI_KWH[k], 2) + " → POA " + fmt(r.monthPOA[k] / E.MDAYS[k], 2) + " kWh/m²/day</span>",
          ev.clientX, ev.clientY);
      });
      hit.addEventListener("mouseleave", hideTip);
    });
    // Table view
    var rows = r.month.map(function (v, k) {
      return "<tr><td>" + E.MNAME[k] + "</td><td>" + fmt(E.GHI_KWH[k], 2) + "</td><td>" + fmt(r.monthPOA[k] / E.MDAYS[k], 2) +
        "</td><td>" + fmt(v) + "</td><td>" + fmt(v * (1 - u)) + "–" + fmt(v * (1 + u)) + "</td><td>" + fmt(v / E.MDAYS[k], 1) + "</td></tr>";
    }).join("");
    $("sy-table").innerHTML = "<thead><tr><th>Month</th><th>GHI (kWh/m²/d)</th><th>POA (kWh/m²/d)</th><th>AC energy (kWh)</th><th>95% UI (kWh)</th><th>kWh/day</th></tr></thead><tbody>" +
      rows + "<tr><td><strong>Year</strong></td><td>" + fmt(E.GHI_KWH.reduce(function (a, v, k) { return a + v * E.MDAYS[k]; }, 0) / 365, 2) +
      "</td><td>" + fmt(r.annualPOA / 365, 2) + "</td><td><strong>" + fmt(r.annual) + "</strong></td><td>" + fmt(r.annual * (1 - u)) + "–" +
      fmt(r.annual * (1 + u)) + "</td><td>" + fmt(r.annual / 365, 1) + "</td></tr></tbody>";
  }

  // ---------------- Figure 4: Sun path ----------------
  var PATH_DAYS = [
    { d: 171, s: "Jun 21" }, { d: 140, s: "May 21 / Jul 21" }, { d: 109, s: "Apr 20 / Aug 22" },
    { d: 78, s: "Mar 20 / Sep 22" }, { d: 50, s: "Feb 20 / Oct 22" }, { d: 20, s: "Jan 21 / Nov 21" }, { d: 354, s: "Dec 21" }
  ];
  function renderSunPath() {
    var P = palette(), o = newSvg($("sy-sunpath"), 320), svg = o.svg, w = o.w, h = o.h;
    var m = { l: 44, r: 16, t: 14, b: 34 }, pw = w - m.l - m.r, ph = h - m.t - m.b;
    var X = function (az) { return m.l + (az - 50) / 260 * pw; }, Y = function (a) { return m.t + ph - a / 80 * ph; };
    [0, 20, 40, 60, 80].forEach(function (a) {
      el("line", { x1: m.l, x2: m.l + pw, y1: Y(a), y2: Y(a), stroke: P.grid }, svg);
      txt(svg, m.l - 6, Y(a) + 4, a + "°", { fill: P.muted, "text-anchor": "end" });
    });
    [[60, "ENE"], [90, "E"], [135, "SE"], [180, "S"], [225, "SW"], [270, "W"], [300, "WNW"]].forEach(function (t) {
      if (w < 520 && (t[0] === 60 || t[0] === 300)) return;
      el("line", { x1: X(t[0]), x2: X(t[0]), y1: m.t, y2: m.t + ph, stroke: P.grid, "stroke-dasharray": t[0] === 180 ? "none" : "2 4" }, svg);
      txt(svg, X(t[0]), h - 12, t[1] + " " + t[0] + "°", { fill: P.muted, "text-anchor": "middle" });
    });
    var hE = state.opts.horizonE;
    if (hE > 0) {
      el("rect", { x: X(50), y: Y(hE), width: X(150) - X(50), height: Y(0) - Y(hE), fill: P.muted, opacity: 0.18 }, svg);
      txt(svg, X(100), Y(hE) - 4, "eastern horizon mask " + hE + "°", { fill: P.muted, "text-anchor": "middle", "font-size": 10 });
    }
    var ramp = P.ramp, n = PATH_DAYS.length;
    PATH_DAYS.forEach(function (pd, k) {
      var c = rampColor(ramp, isDark() ? 0.95 - k / (n - 1) * 0.6 : 0.35 + k / (n - 1) * 0.65);
      var col = "rgb(" + c.join(",") + ")", d = "", started = false;
      for (var t = 3; t <= 21; t += 0.05) {
        var p = E.sunPos(pd.d, t);
        if (p.alt < 0 || p.az < 50 || p.az > 310) { started = false; continue; }
        d += (started ? "L" : "M") + X(p.az).toFixed(1) + " " + Y(p.alt).toFixed(1); started = true;
      }
      el("path", { d: d, fill: "none", stroke: col, "stroke-width": 2 }, svg);
      // hourly dots in clock time
      var dst = E.isDST(pd.d) ? 1 : 0;
      for (var hc = 5; hc <= 21; hc++) {
        var p2 = E.sunPos(pd.d, hc - dst);
        if (p2.alt <= 0 || p2.az < 50 || p2.az > 310) continue;
        var dot = el("circle", { cx: X(p2.az), cy: Y(p2.alt), r: 3.2, fill: col, stroke: P.bg, "stroke-width": 1.5 }, svg);
        (function (p2, hc, pd) {
          dot.addEventListener("mousemove", function (ev) {
            showTip("<strong>" + pd.s + ", " + hhmm(hc) + (dst ? " MDT" : " MST") + "</strong><br>Altitude " + fmt(p2.alt, 1) +
              "°, azimuth " + fmt(p2.az, 0) + "°", ev.clientX, ev.clientY);
          });
          dot.addEventListener("mouseleave", hideTip);
        })(p2, hc, pd);
      }
      var pn = E.sunPos(pd.d, E.TIMES[pd.d].noon);
      pd.col = col;
      if (k === 0) txt(svg, X(180) - 8, Y(pn.alt) - 6, pd.s, { fill: P.text, "font-size": 10.5, "font-weight": 600, "text-anchor": "end" });
      if (k === 3) txt(svg, X(180) + 8, Y(pn.alt) - 6, pd.s, { fill: P.text, "font-size": 10.5, "font-weight": 600 });
      if (k === n - 1) txt(svg, X(180), Y(pn.alt) + 16, pd.s, { fill: P.text, "font-size": 10.5, "font-weight": 600, "text-anchor": "middle" });
    });
    $("sy-sunpath-legend").innerHTML = PATH_DAYS.map(function (pd) {
      return '<span class="sy-leg-item"><span class="sy-sw" style="border-top:2px solid ' + pd.col + '"></span>' + pd.s + "</span>";
    }).join("") + '<span class="sy-leg-item"><span class="sy-sw-ring" style="border-color:' + P.s1 + '"></span>Roof surface normal</span>';
    // Roof normal
    var nAlt = 90 - state.pitch;
    el("circle", { cx: X(180), cy: Y(Math.min(80, nAlt)), r: 6, fill: "none", stroke: P.s1, "stroke-width": 2.5 }, svg);
    txt(svg, X(180) + 10, Y(Math.min(80, nAlt)) + 4, nAlt > 80 ? "roof normal ↑ " + nAlt + "°" : "roof normal (" + nAlt + "°)",
      { fill: P.s1, "text-anchor": "start", "font-weight": 600, "font-size": 10.5 });
  }

  // ---------------- Figure 5: Annual yield vs pitch ----------------
  function renderSweep() {
    var sw = state.sweep, xs = sw.map(function (s) { return s.p; }), ys = sw.map(function (s) { return s.e; });
    var lo = Math.min.apply(null, ys), hi = Math.max.apply(null, ys), pad = (hi - lo) * 0.25 + 1;
    var cur = sw[Math.round(state.pitch)];
    var marks = [{ x: state.best.pitch, y: state.best.annual, color: "s2", label: "optimum " + state.best.pitch + "°", anchor: "end" }];
    if (Math.abs(state.pitch - state.best.pitch) >= 2) marks.push({ x: cur.p, y: cur.e, color: "s1", label: state.pitch + "° · " + fmt(cur.e) + " kWh", below: true });
    lineChart($("sy-sweep"), {
      h: 260, xs: xs, x0: 0, x1: 45, y0: Math.max(0, lo - pad), y1: hi + pad * 0.6, yLabel: "kWh / year",
      xTicks: [0, 5, 10, 15, 20, 25, 30, 35, 40, 45].map(function (v) { return { v: v, s: v + "°" }; }),
      zones: [{ a: 7, b: 10, label: "7–10° structural" }],
      series: [{ ys: ys, color: "s1" }], marks: marks,
      tip: function (k) {
        return "<strong>" + xs[k] + "° pitch</strong><br>" + fmt(ys[k]) + " kWh/yr (95% UI " + fmt(ys[k] * (1 - E.REL_UNC)) + "–" +
          fmt(ys[k] * (1 + E.REL_UNC)) + ")<br>" + fmt(100 * ys[k] / state.best.annual, 1) + "% of optimum" +
          (xs[k] > 10 ? "<br><span class='sy-tip-muted'>Outside the 7–10° structural/zoning window</span>" : "");
      }
    });
  }

  function renderAll() {
    if (!$("solar").classList.contains("active") || !state.res) return;
    renderHeatmap(); renderDaily(); renderMonthly(); renderSunPath(); renderSweep();
  }

  function init() {
    if (!$("solar")) return;
    var calcPitch = $("calc-pitch");
    if (calcPitch) $("sy-pitch").value = calcPitch.value;
    var deb;
    function queue() { clearTimeout(deb); deb = setTimeout(compute, 60); }
    ["sy-pitch", "sy-watt", "sy-loss", "sy-horizon", "sy-snow", "sy-nmod", "sy-auto"].forEach(function (id) {
      $(id).addEventListener("input", function () {
        if (id === "sy-pitch" && calcPitch) { calcPitch.value = $("sy-pitch").value; calcPitch.dispatchEvent(new Event("input")); }
        queue();
      });
      $(id).addEventListener("change", queue);
    });
    // Keep in sync with the Interactive Sizer
    ["calc-depth", "calc-span", "calc-pitch", "calc-mod-overhang"].forEach(function (id) {
      var e = $(id); if (!e) return;
      e.addEventListener("input", function () { if (id === "calc-pitch") $("sy-pitch").value = e.value; queue(); });
    });
    document.querySelectorAll(".preset-chip").forEach(function (b) {
      b.addEventListener("click", function () { if (calcPitch) $("sy-pitch").value = calcPitch.value; queue(); });
    });
    document.querySelectorAll('.tab-button[data-tab="solar"]').forEach(function (b) {
      b.addEventListener("click", function () { setTimeout(renderAll, 20); });
    });
    var rz; window.addEventListener("resize", function () { clearTimeout(rz); rz = setTimeout(renderAll, 120); });
    new MutationObserver(function () { setTimeout(renderAll, 30); })
      .observe(document.documentElement, { attributes: true, attributeFilter: ["data-theme"] });
    compute();
  }
  if (document.readyState === "loading") document.addEventListener("DOMContentLoaded", init); else init();
})();
