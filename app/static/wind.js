/* HK Wind Sketch: a qualitative 2D wind sketch at pedestrian level for Hong Kong.
 *
 * A D2Q9 lattice Boltzmann solver (BGK + Smagorinsky) on a square slice about
 * 2 m above ground, through LandsD building footprints. WebGL2 (float render
 * targets) when available, plain JavaScript otherwise. Everything runs in the
 * browser: no server, no upload.
 *
 * Not an Air Ventilation Assessment (HPLB/ETWB TC No. 1/06). No VR numbers.
 * Open source, MIT licence. Data: Lands Department and Planning Department,
 * HKSAR Government, via DATA.GOV.HK / CSDI.
 */
(function () {
  'use strict';
  const root = document.getElementById('hkw');
  if (!root) return;
  const $ = (s) => root.querySelector(s);
  const params = new URLSearchParams(location.search);

  // ---------------------------------------------------------------- constants
  const DOMAIN_M = 600;            // simulated square, metres
  const U0 = 0.1;                  // incoming speed in lattice units
  const NU = 0.004;                // lattice viscosity (turbulence added by Smagorinsky)
  const TAU0 = 0.5 + 3 * NU;
  const CS = 0.18;                 // Smagorinsky constant
  const BLOCK_MIN_H = 2.5;         // a block at least this tall (m, from its own base) stands in the 2 m slice
  const FPS_CAP = 30;
  const IDLE_PAUSE_MS = 90000;     // pause when nobody has touched the page for this long
  const QUALITY = { light: 240, normal: 300, fine: 400 };
  const CPU_N = 150;
  const DIRS = [
    ['N', '北'], ['NNE', '北北東'], ['NE', '東北'], ['ENE', '東北東'],
    ['E', '東'], ['ESE', '東南東'], ['SE', '東南'], ['SSE', '南南東'],
    ['S', '南'], ['SSW', '南南西'], ['SW', '西南'], ['WSW', '西南西'],
    ['W', '西'], ['WNW', '西北西'], ['NW', '西北'], ['NNW', '北北西'],
  ];
  const SEASONS = { annual: '全年 Annual', summer: '夏季 Summer', winter: '冬季 Winter' };
  const EX = [0, 1, 0, -1, 0, 1, -1, -1, 1];
  const EY = [0, 0, 1, 0, -1, 1, 1, -1, -1];      // +y = next row = south on screen
  const WT = [4 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 9, 1 / 36, 1 / 36, 1 / 36, 1 / 36];
  const OPP = [0, 3, 4, 1, 2, 7, 8, 5, 6];

  // live open data (all answer cross-origin requests, no key)
  const WFS = 'https://portal.csdi.gov.hk/server/services/common/landsd_rcd_1637211194312_35158/MapServer/WFSServer';
  const PLACE_SEARCH = 'https://www.map.gov.hk/gs/api/v1.0.0/locationSearch?q=';
  const ALS = 'https://www.als.gov.hk/lookup?n=5&q=';
  const STATIC = window.HKW_STATIC || 'app/static/';
  const DATA_URL = window.HKW_DATA_URL || (STATIC + 'data/');
  const HK_BOUNDS = [22.13, 113.82, 22.58, 114.45];     // lat0, lon0, lat1, lon1
  const FETCH_PAD_M = 20;
  const M_PER_DEG_LAT = 110574;
  const mPerDegLon = (lat) => 111320 * Math.cos(lat * Math.PI / 180);

  // HK1980 Grid (EPSG:2326) <-> WGS84. Transverse Mercator on the International 1924
  // ellipsoid plus LandsD's published approximate datum shift (lat -5.5", lon +8.8");
  // good to about 1-2 m, plenty for 500 m wind cells and a 600 m square.
  const HK80 = (() => {
    const a = 6378388, f = 1 / 297, e2 = 2 * f - f * f, ep2 = e2 / (1 - e2), e4 = e2 * e2, e6 = e4 * e2, Dg = Math.PI / 180;
    const lat0 = (22 + 18 / 60 + 43.68 / 3600) * Dg, lon0 = (114 + 10 / 60 + 42.80 / 3600) * Dg, FE = 836694.05, FN = 819069.80;
    const M = (p) => a * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * p - (3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * Math.sin(2 * p)
      + (15 * e4 / 256 + 45 * e6 / 1024) * Math.sin(4 * p) - (35 * e6 / 3072) * Math.sin(6 * p));
    const M0 = M(lat0), SLAT = 5.5 / 3600, SLON = 8.8 / 3600;
    return {
      toGrid(latW, lonW) {
        const p = (latW + SLAT) * Dg, l = (lonW - SLON) * Dg, s = Math.sin(p), c = Math.cos(p), t = Math.tan(p);
        const N = a / Math.sqrt(1 - e2 * s * s), T = t * t, C = ep2 * c * c, A = (l - lon0) * c;
        const x = N * (A + (1 - T + C) * A ** 3 / 6 + (5 - 18 * T + T * T + 72 * C - 58 * ep2) * A ** 5 / 120);
        const y = M(p) - M0 + N * t * (A * A / 2 + (5 - T + 9 * C + 4 * C * C) * A ** 4 / 24 + (61 - 58 * T + T * T + 600 * C - 330 * ep2) * A ** 6 / 720);
        return [FE + x, FN + y];
      },
      toWGS(E, Nn) {
        const mu = (M0 + Nn - FN) / (a * (1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256));
        const e1 = (1 - Math.sqrt(1 - e2)) / (1 + Math.sqrt(1 - e2));
        const p1 = mu + (3 * e1 / 2 - 27 * e1 ** 3 / 32) * Math.sin(2 * mu) + (21 * e1 * e1 / 16 - 55 * e1 ** 4 / 32) * Math.sin(4 * mu)
          + (151 * e1 ** 3 / 96) * Math.sin(6 * mu) + (1097 * e1 ** 4 / 512) * Math.sin(8 * mu);
        const s = Math.sin(p1), c = Math.cos(p1), t = Math.tan(p1);
        const C1 = ep2 * c * c, T1 = t * t, N1 = a / Math.sqrt(1 - e2 * s * s), R1 = a * (1 - e2) / Math.pow(1 - e2 * s * s, 1.5), d = (E - FE) / N1;
        const p = p1 - (N1 * t / R1) * (d * d / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * ep2) * d ** 4 / 24 + (61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * ep2 - 3 * C1 * C1) * d ** 6 / 720);
        const l = lon0 + (d - (1 + 2 * T1 + C1) * d ** 3 / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * ep2 + 24 * T1 * T1) * d ** 5 / 120) / c;
        return [p / Dg - SLAT, l / Dg + SLON];
      },
    };
  })();

  // Wind FROM bearing b (deg clockwise from north) blows towards b+180.
  function inflow(dirIdx) {
    const b = dirIdx * 22.5 * Math.PI / 180;
    return [-Math.sin(b) * U0, Math.cos(b) * U0];
  }

  // ---------------------------------------------------------------- GPU solver
  const VS = `#version 300 es
in vec2 p; void main(){ gl_Position = vec4(p, 0.0, 1.0); }`;
  const FS = `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D T0; uniform sampler2D T1; uniform sampler2D T2; uniform sampler2D M;
uniform vec2 uIn; uniform float tau0; uniform float cs2; uniform int N; uniform float sponge; uniform int uInit;
layout(location=0) out vec4 o0; layout(location=1) out vec4 o1; layout(location=2) out vec4 o2;
const ivec2 E[9] = ivec2[9](ivec2(0,0),ivec2(1,0),ivec2(0,1),ivec2(-1,0),ivec2(0,-1),ivec2(1,1),ivec2(-1,1),ivec2(-1,-1),ivec2(1,-1));
const float W[9] = float[9](4./9.,1./9.,1./9.,1./9.,1./9.,1./36.,1./36.,1./36.,1./36.);
const int OPP[9] = int[9](0,3,4,1,2,7,8,5,6);
float getF(ivec2 q, int i){
  if (i < 4) return texelFetch(T0, q, 0)[i];
  if (i < 8) return texelFetch(T1, q, 0)[i-4];
  return texelFetch(T2, q, 0).x;
}
float feq(int i, float rho, vec2 u){ float eu = dot(vec2(E[i]), u); return W[i]*rho*(1.+3.*eu+4.5*eu*eu-1.5*dot(u,u)); }
bool solidAt(ivec2 q){ return texelFetch(M, q, 0).r > 0.5; }
void writeEq(float rho, vec2 u){
  o0 = vec4(feq(0,rho,u), feq(1,rho,u), feq(2,rho,u), feq(3,rho,u));
  o1 = vec4(feq(4,rho,u), feq(5,rho,u), feq(6,rho,u), feq(7,rho,u));
  o2 = vec4(feq(8,rho,u), u.x, u.y, rho);
}
void main(){
  ivec2 p = ivec2(gl_FragCoord.xy);
  if (solidAt(p)) { writeEq(1., vec2(0.)); o2.w = -1.; return; }
  if (uInit == 1) { writeEq(1., uIn); return; }
  float f[9];
  for (int i = 0; i < 9; i++) {
    ivec2 q = p - E[i];
    if (q.x < 0 || q.y < 0 || q.x >= N || q.y >= N) f[i] = feq(i, 1., uIn);
    else if (solidAt(q)) f[i] = getF(p, OPP[i]);
    else f[i] = getF(q, i);
  }
  float rho = 0.; vec2 m = vec2(0.);
  for (int i = 0; i < 9; i++) { rho += f[i]; m += vec2(E[i]) * f[i]; }
  if (!(rho > 0.2 && rho < 5.0)) { writeEq(1., uIn); return; }
  vec2 u = m / rho; float sp = length(u); if (sp > 0.3) u *= 0.3 / sp;
  float fe[9]; float pxx = 0., pyy = 0., pxy = 0.;
  for (int i = 0; i < 9; i++) {
    fe[i] = feq(i, rho, u); float d = f[i] - fe[i]; vec2 e = vec2(E[i]);
    pxx += e.x*e.x*d; pyy += e.y*e.y*d; pxy += e.x*e.y*d;
  }
  float Q = sqrt(pxx*pxx + pyy*pyy + 2.*pxy*pxy);
  float tau = 0.5*(tau0 + sqrt(tau0*tau0 + 18.*1.41421356*cs2*Q/rho));
  float om = 1./tau;
  float dd = float(min(min(p.x, p.y), min(N-1-p.x, N-1-p.y)));
  float s = 0.; if (dd < sponge) { float t = 1. - dd/sponge; s = 0.25*t*t; }
  float g[9];
  for (int i = 0; i < 9; i++) { float post = f[i] - om*(f[i]-fe[i]); g[i] = mix(post, feq(i, 1., uIn), s); }
  o0 = vec4(g[0], g[1], g[2], g[3]); o1 = vec4(g[4], g[5], g[6], g[7]); o2 = vec4(g[8], u.x, u.y, rho);
}`;

  function makeGPU(N) {
    const cv = document.createElement('canvas');
    cv.width = cv.height = 1;
    let gl;
    try {
      // low-power: on a laptop with two GPUs this asks for the integrated one.
      gl = cv.getContext('webgl2', { antialias: false, depth: false, stencil: false, powerPreference: 'low-power' });
    } catch (e) { gl = null; }
    if (!gl || !gl.getExtension('EXT_color_buffer_float')) return null;
    const sh = (type, src) => {
      const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      return s;
    };
    const prog = gl.createProgram();
    try {
      gl.attachShader(prog, sh(gl.VERTEX_SHADER, VS));
      gl.attachShader(prog, sh(gl.FRAGMENT_SHADER, FS));
      gl.linkProgram(prog);
      if (!gl.getProgramParameter(prog, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(prog));
    } catch (e) { console.warn('[wind] GPU shader failed, using CPU:', e.message); return null; }
    gl.useProgram(prog);
    const buf = gl.createBuffer();
    gl.bindBuffer(gl.ARRAY_BUFFER, buf);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    const loc = gl.getAttribLocation(prog, 'p');
    gl.enableVertexAttribArray(loc);
    gl.vertexAttribPointer(loc, 2, gl.FLOAT, false, 0, 0);
    const U = (n) => gl.getUniformLocation(prog, n);
    const uni = { uIn: U('uIn'), tau0: U('tau0'), cs2: U('cs2'), N: U('N'), sponge: U('sponge'), uInit: U('uInit') };
    gl.uniform1i(U('T0'), 0); gl.uniform1i(U('T1'), 1); gl.uniform1i(U('T2'), 2); gl.uniform1i(U('M'), 3);
    gl.uniform1f(uni.tau0, TAU0); gl.uniform1f(uni.cs2, CS * CS); gl.uniform1i(uni.N, N);
    gl.uniform1f(uni.sponge, Math.round(N * 0.06));
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    const mk = (internal, fmt, type, data) => {
      const t = gl.createTexture(); gl.bindTexture(gl.TEXTURE_2D, t);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MIN_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_MAG_FILTER, gl.NEAREST);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE);
      gl.texParameteri(gl.TEXTURE_2D, gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE);
      gl.texImage2D(gl.TEXTURE_2D, 0, internal, N, N, 0, fmt, type, data);
      return t;
    };
    const sets = [0, 1].map(() => {
      const t = [0, 1, 2].map(() => mk(gl.RGBA32F, gl.RGBA, gl.FLOAT, null));
      const fb = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
      t.forEach((tx, k) => gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + k, gl.TEXTURE_2D, tx, 0));
      gl.drawBuffers([gl.COLOR_ATTACHMENT0, gl.COLOR_ATTACHMENT1, gl.COLOR_ATTACHMENT2]);
      return { t, fb, ok: gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE };
    });
    if (!sets[0].ok || !sets[1].ok) return null;
    const mask = mk(gl.R8, gl.RED, gl.UNSIGNED_BYTE, new Uint8Array(N * N));
    let cur = 0, lost = false;
    cv.addEventListener('webglcontextlost', (e) => { e.preventDefault(); lost = true; });
    const out = new Float32Array(N * N * 4);
    function pass(init) {
      const src = sets[cur], dst = sets[1 - cur];
      gl.bindFramebuffer(gl.FRAMEBUFFER, dst.fb);
      for (let k = 0; k < 3; k++) { gl.activeTexture(gl.TEXTURE0 + k); gl.bindTexture(gl.TEXTURE_2D, src.t[k]); }
      gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, mask);
      gl.uniform1i(uni.uInit, init ? 1 : 0);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      cur = 1 - cur;
    }
    gl.viewport(0, 0, N, N);
    return {
      kind: 'gpu', N,
      get lost() { return lost || gl.isContextLost(); },
      init(u) { gl.uniform2f(uni.uIn, u[0], u[1]); pass(true); },
      setInflow(u) { gl.uniform2f(uni.uIn, u[0], u[1]); },
      setMask(m) {
        const b = new Uint8Array(N * N); for (let i = 0; i < b.length; i++) b[i] = m[i] ? 255 : 0;
        gl.activeTexture(gl.TEXTURE3); gl.bindTexture(gl.TEXTURE_2D, mask);
        gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, N, N, gl.RED, gl.UNSIGNED_BYTE, b);
      },
      step(n) { for (let k = 0; k < n; k++) pass(false); },
      read() {
        gl.bindFramebuffer(gl.FRAMEBUFFER, sets[cur].fb);
        gl.readBuffer(gl.COLOR_ATTACHMENT2);
        gl.readPixels(0, 0, N, N, gl.RGBA, gl.FLOAT, out);
        return out;   // per cell: f8, ux, uy, rho (rho = -1 inside buildings)
      },
      dispose() { const x = gl.getExtension('WEBGL_lose_context'); if (x) x.loseContext(); },
    };
  }

  // ---------------------------------------------------------------- CPU solver (same scheme)
  function makeCPU(N) {
    const NN = N * N;
    let f = new Float32Array(9 * NN), g = new Float32Array(9 * NN);
    let mask = new Uint8Array(NN);
    const out = new Float32Array(NN * 4);
    let ux0 = 0, uy0 = 0;
    const sponge = Math.round(N * 0.06);
    const feq = (i, rho, ux, uy) => {
      const eu = EX[i] * ux + EY[i] * uy;
      return WT[i] * rho * (1 + 3 * eu + 4.5 * eu * eu - 1.5 * (ux * ux + uy * uy));
    };
    const fin = new Float64Array(9), fe = new Float64Array(9), fInf = new Float64Array(9);
    function setEq(arr, c, rho, ux, uy) { for (let i = 0; i < 9; i++) arr[i * NN + c] = feq(i, rho, ux, uy); }
    function step1() {
      for (let i = 0; i < 9; i++) fInf[i] = feq(i, 1, ux0, uy0);
      for (let y = 0; y < N; y++) {
        for (let x = 0; x < N; x++) {
          const c = y * N + x;
          if (mask[c]) { setEq(g, c, 1, 0, 0); out[c * 4 + 1] = 0; out[c * 4 + 2] = 0; out[c * 4 + 3] = -1; continue; }
          let rho = 0, mx = 0, my = 0;
          for (let i = 0; i < 9; i++) {
            const qx = x - EX[i], qy = y - EY[i];
            let v;
            if (qx < 0 || qy < 0 || qx >= N || qy >= N) v = fInf[i];
            else { const q = qy * N + qx; v = mask[q] ? f[OPP[i] * NN + c] : f[i * NN + q]; }
            fin[i] = v; rho += v; mx += EX[i] * v; my += EY[i] * v;
          }
          if (!(rho > 0.2 && rho < 5)) { setEq(g, c, 1, ux0, uy0); continue; }
          let ux = mx / rho, uy = my / rho;
          const sp = Math.hypot(ux, uy); if (sp > 0.3) { ux *= 0.3 / sp; uy *= 0.3 / sp; }
          let pxx = 0, pyy = 0, pxy = 0;
          for (let i = 0; i < 9; i++) {
            fe[i] = feq(i, rho, ux, uy); const d = fin[i] - fe[i];
            pxx += EX[i] * EX[i] * d; pyy += EY[i] * EY[i] * d; pxy += EX[i] * EY[i] * d;
          }
          const Q = Math.sqrt(pxx * pxx + pyy * pyy + 2 * pxy * pxy);
          const tau = 0.5 * (TAU0 + Math.sqrt(TAU0 * TAU0 + 18 * Math.SQRT2 * CS * CS * Q / rho));
          const om = 1 / tau;
          const dd = Math.min(x, y, N - 1 - x, N - 1 - y);
          let s = 0; if (dd < sponge) { const t = 1 - dd / sponge; s = 0.25 * t * t; }
          for (let i = 0; i < 9; i++) {
            const post = fin[i] - om * (fin[i] - fe[i]);
            g[i * NN + c] = post + s * (fInf[i] - post);
          }
          out[c * 4 + 1] = ux; out[c * 4 + 2] = uy; out[c * 4 + 3] = rho;
        }
      }
      const t = f; f = g; g = t;
    }
    return {
      kind: 'cpu', N, lost: false,
      init(u) {
        ux0 = u[0]; uy0 = u[1];
        for (let c = 0; c < NN; c++) {
          if (mask[c]) { setEq(f, c, 1, 0, 0); out[c * 4 + 3] = -1; out[c * 4 + 1] = out[c * 4 + 2] = 0; }
          else { setEq(f, c, 1, ux0, uy0); out[c * 4 + 1] = ux0; out[c * 4 + 2] = uy0; out[c * 4 + 3] = 1; }
        }
      },
      setInflow(u) { ux0 = u[0]; uy0 = u[1]; },
      setMask(m) { mask = Uint8Array.from(m, (v) => (v ? 1 : 0)); },
      step(n) { for (let k = 0; k < n; k++) step1(); },
      read() { return out; },
      dispose() {},
    };
  }

  // ---------------------------------------------------------------- colour map
  const STOPS = [   // speed relative to incoming wind -> colour
    [0.00, [24, 44, 84]], [0.35, [38, 104, 160]], [0.70, [62, 165, 150]],
    [1.00, [228, 198, 78]], [1.30, [236, 122, 44]], [1.70, [176, 36, 30]],
  ];
  const LUT_MAX = 1.7, LUT = new Uint8ClampedArray(256 * 3);
  for (let k = 0; k < 256; k++) {
    const v = (k / 255) * LUT_MAX;
    let a = STOPS[0], b = STOPS[STOPS.length - 1];
    for (let s = 0; s < STOPS.length - 1; s++) if (v >= STOPS[s][0] && v <= STOPS[s + 1][0]) { a = STOPS[s]; b = STOPS[s + 1]; break; }
    const t = b[0] === a[0] ? 0 : (v - a[0]) / (b[0] - a[0]);
    for (let c = 0; c < 3; c++) LUT[k * 3 + c] = a[1][c] + (b[1][c] - a[1][c]) * t;
  }
  const lg = $('#hkwLegendBar');
  if (lg) lg.style.background = 'linear-gradient(90deg,' + STOPS.map((s) => `rgb(${s[1].join(',')}) ${(s[0] / LUT_MAX * 100).toFixed(0)}%`).join(',') + ')';

  // ---------------------------------------------------------------- DOM
  const stage = $('#hkwStage');
  const cvField = $('#hkwField'), cvTrail = $('#hkwTrail'), cvBld = $('#hkwBld');
  const ctxF = cvField.getContext('2d'), ctxT = cvTrail.getContext('2d'), ctxB = cvBld.getContext('2d');
  const tip = $('#hkwTip');
  const statusEl = $('#hkwStatus'), engineEl = $('#hkwEngine'), settleEl = $('#hkwSettle');
  const runBtn = $('#hkwRun'), pauseNote = $('#hkwPauseNote');

  const S = {
    data: null, areaId: params.get('lat') ? null : (params.get('area') || 'kwun_tong'), dirIdx: null, season: 'annual',
    view: '2d', loading: false, siteLabel: '',
    quality: params.get('quality') || 'normal', running: true, idlePaused: false,
    showField: true, showParticles: true, mode: 'look', user: [], removed: new Set(), sel: null,
    steps: 0, lastTouch: performance.now(), spf: 20, size: 600, dpr: 1, field: null, fps: 0,
  };
  let solver = null, maskArr = null, particles = [], buildings = [], nextUserId = 1;
  let view3d = null, view3dLoading = null, need3d = { site: true, user: true };
  window.HKW = { state: S, get solver() { return solver; }, get view3d() { return view3d; }, HK80 };   // for tests and the curious

  function forceCPU() { return params.get('engine') === 'cpu'; }

  function newSolver() {
    if (solver) solver.dispose();
    const N = QUALITY[S.quality] || 300;
    solver = forceCPU() ? null : makeGPU(N);
    if (!solver) solver = makeCPU(CPU_N);
    S.spf = solver.kind === 'gpu' ? 20 : 3;
    engineEl.textContent = solver.kind === 'gpu'
      ? `Graphics card (WebGL2) · ${solver.N} x ${solver.N} cells, ${(DOMAIN_M / solver.N).toFixed(1)} m each`
      : `Processor only (no WebGL2 float support) · ${solver.N} x ${solver.N} cells, ${(DOMAIN_M / solver.N).toFixed(1)} m each`;
    $('#hkwQuality').disabled = solver.kind !== 'gpu';
  }

  // ---------------------------------------------------------------- geometry
  function userPoly(u) {
    const a = u.rot * Math.PI / 180, c = Math.cos(a), s = Math.sin(a);
    return [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([px, py]) => {
      const x = px * u.w / 2, y = py * u.d / 2;
      return [u.cx + x * c - y * s, u.cy + x * s + y * c];
    });
  }
  const blocks = (b) => b.h >= BLOCK_MIN_H && b.type !== 'Open-sided Structure';
  function pathRing(ctx, ring, k) {
    ctx.moveTo(ring[0][0] * k, ring[0][1] * k);
    for (let i = 1; i < ring.length; i++) ctx.lineTo(ring[i][0] * k, ring[i][1] * k);
    ctx.closePath();
  }
  function rebuildMask() {
    const N = solver.N;
    const oc = rebuildMask.c || (rebuildMask.c = document.createElement('canvas'));
    oc.width = oc.height = N;
    const c = oc.getContext('2d', { willReadFrequently: true });
    c.clearRect(0, 0, N, N); c.fillStyle = '#000';
    const k = N / DOMAIN_M;
    c.beginPath();
    for (const b of buildings) if (blocks(b) && !S.removed.has(b.id)) for (const r of b.rings) pathRing(c, r, k);
    for (const u of S.user) pathRing(c, userPoly(u), k);
    c.fill('nonzero');
    need3d.user = true;
    const px = c.getImageData(0, 0, N, N).data;
    maskArr = new Uint8Array(N * N);
    for (let i = 0; i < N * N; i++) maskArr[i] = px[i * 4 + 3] >= 128 ? 1 : 0;
    solver.setMask(maskArr);
  }
  function inPoly(x, y, ring) {
    let ins = false;
    for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
      const xi = ring[i][0], yi = ring[i][1], xj = ring[j][0], yj = ring[j][1];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) ins = !ins;
    }
    return ins;
  }
  function hitUser(x, y) { for (let i = S.user.length - 1; i >= 0; i--) if (inPoly(x, y, userPoly(S.user[i]))) return S.user[i]; return null; }
  function hitBuilding(x, y) {
    let best = null;
    for (const b of buildings) {
      if (x < b.bb[0] || x > b.bb[2] || y < b.bb[1] || y > b.bb[3]) continue;
      if (b.rings.some((r) => inPoly(x, y, r)) && (!best || b.h > best.h)) best = b;
    }
    return best;
  }

  // ---------------------------------------------------------------- drawing
  function resize() {
    const w = Math.min(stage.parentElement.clientWidth, 760);
    S.size = Math.max(260, Math.floor(w));
    S.dpr = Math.min(window.devicePixelRatio || 1, 2);
    stage.style.width = stage.style.height = S.size + 'px';
    for (const cv of [cvTrail, cvBld]) { cv.width = cv.height = Math.round(S.size * S.dpr); }
    if (view3d) view3d.resize(S.size, S.size);
    drawBuildings();
  }
  function heightTone(h, dark) {
    const t = Math.min(1, Math.max(0, (h - 5) / 150));
    const a = dark ? [120, 117, 110] : [214, 209, 199], b = dark ? [215, 210, 200] : [92, 89, 84];
    return `rgb(${a.map((v, i) => Math.round(v + (b[i] - v) * t)).join(',')})`;
  }
  function isDark() {
    const t = document.documentElement.getAttribute('data-theme');
    if (t) return t === 'dark';
    return window.matchMedia && matchMedia('(prefers-color-scheme: dark)').matches;
  }
  function drawBuildings() {
    need3d.user = true;
    const c = ctxB, k = S.size * S.dpr / DOMAIN_M, dark = isDark();
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.clearRect(0, 0, cvBld.width, cvBld.height);
    c.lineJoin = 'round';
    const sorted = buildings.slice().sort((a, b) => a.h - b.h);
    for (const b of sorted) {
      c.beginPath(); for (const r of b.rings) pathRing(c, r, k);
      if (S.removed.has(b.id)) {
        c.setLineDash([4 * S.dpr, 3 * S.dpr]); c.strokeStyle = '#c0392b'; c.lineWidth = 1.2 * S.dpr; c.stroke(); c.setLineDash([]);
      } else if (!blocks(b)) {
        c.setLineDash([2 * S.dpr, 2 * S.dpr]); c.strokeStyle = dark ? 'rgba(255,255,255,.5)' : 'rgba(0,0,0,.45)'; c.lineWidth = 0.8 * S.dpr; c.stroke(); c.setLineDash([]);
      } else {
        c.fillStyle = heightTone(b.h, dark); c.fill('nonzero');
        c.strokeStyle = dark ? 'rgba(0,0,0,.6)' : 'rgba(30,28,25,.75)'; c.lineWidth = 0.7 * S.dpr; c.stroke();
      }
    }
    for (const u of S.user) {
      c.beginPath(); pathRing(c, userPoly(u), k);
      c.fillStyle = 'rgba(217,119,87,.92)'; c.fill();
      c.strokeStyle = u === S.sel ? '#111' : '#7a3a22'; c.lineWidth = (u === S.sel ? 2.4 : 1.2) * S.dpr; c.stroke();
      c.fillStyle = '#fff'; c.font = `600 ${11 * S.dpr}px ui-sans-serif,system-ui,sans-serif`; c.textAlign = 'center'; c.textBaseline = 'middle';
      c.fillText(`${u.h} m`, u.cx * k, u.cy * k);
    }
    if (S.draft) {
      c.beginPath(); pathRing(c, userPoly(S.draft), k);
      c.setLineDash([5 * S.dpr, 4 * S.dpr]); c.strokeStyle = '#d97757'; c.lineWidth = 2 * S.dpr; c.stroke(); c.setLineDash([]);
    }
  }

  // ---------------------------------------------------------------- field + particles
  function paintField(buf) {
    const N = solver.N;
    if (!S.field || S.field.width !== N) { cvField.width = cvField.height = N; S.field = ctxF.createImageData(N, N); }
    const d = S.field.data, inv = 1 / U0, sc = 255 / LUT_MAX;
    for (let c = 0, n = N * N; c < n; c++) {
      const o = c * 4;
      if (buf[o + 3] < 0) { d[o + 3] = 0; continue; }
      const sp = Math.hypot(buf[o + 1], buf[o + 2]) * inv;
      const li = Math.min(255, (sp * sc) | 0) * 3;
      d[o] = LUT[li]; d[o + 1] = LUT[li + 1]; d[o + 2] = LUT[li + 2]; d[o + 3] = 255;
    }
    ctxF.putImageData(S.field, 0, 0);
  }
  function sampleU(buf, xm, ym) {
    const N = solver.N, dx = DOMAIN_M / N;
    const gx = xm / dx - 0.5, gy = ym / dx - 0.5;
    const x0 = Math.max(0, Math.min(N - 2, Math.floor(gx))), y0 = Math.max(0, Math.min(N - 2, Math.floor(gy)));
    const tx = Math.min(1, Math.max(0, gx - x0)), ty = Math.min(1, Math.max(0, gy - y0));
    let ux = 0, uy = 0, solid = false;
    const cs = [[x0, y0, (1 - tx) * (1 - ty)], [x0 + 1, y0, tx * (1 - ty)], [x0, y0 + 1, (1 - tx) * ty], [x0 + 1, y0 + 1, tx * ty]];
    for (const [cx, cy, w] of cs) {
      const o = (cy * N + cx) * 4;
      if (buf[o + 3] < 0) { if (w > 0.25) solid = true; continue; }
      ux += buf[o + 1] * w; uy += buf[o + 2] * w;
    }
    return solid ? null : [ux, uy];
  }
  function seedParticle(p) {
    for (let tries = 0; tries < 20; tries++) {
      const x = Math.random() * DOMAIN_M, y = Math.random() * DOMAIN_M;
      const N = solver.N, cx = Math.floor(x / DOMAIN_M * N), cy = Math.floor(y / DOMAIN_M * N);
      if (!maskArr || !maskArr[cy * N + cx]) { p.x = x; p.y = y; p.age = Math.floor(Math.random() * 80); p.fresh = true; return p; }
    }
    p.x = Math.random() * DOMAIN_M; p.y = Math.random() * DOMAIN_M; p.age = 0; p.fresh = true; return p;
  }
  function resetParticles() {
    const n = solver.kind === 'gpu' ? 2600 : 1400;
    particles = Array.from({ length: n }, () => seedParticle({}));
    ctxT.setTransform(1, 0, 0, 1, 0, 0); ctxT.clearRect(0, 0, cvTrail.width, cvTrail.height);
  }
  function moveParticles(buf, steps) {
    const c = ctxT, k = S.size * S.dpr / DOMAIN_M, dx = DOMAIN_M / solver.N;
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.globalCompositeOperation = 'destination-out';
    c.fillStyle = 'rgba(0,0,0,0.10)'; c.fillRect(0, 0, cvTrail.width, cvTrail.height);
    c.globalCompositeOperation = 'source-over';
    if (!S.showParticles && S.view !== '3d') return;
    const draw2d = S.view === '2d';
    c.strokeStyle = S.showField ? 'rgba(255,255,255,0.78)' : (isDark() ? 'rgba(235,232,225,.8)' : 'rgba(30,30,30,.7)');
    c.lineWidth = 1.1 * S.dpr; c.beginPath();
    const gain = steps * dx * 0.55;
    for (const p of particles) {
      const u = sampleU(buf, p.x, p.y);
      p.age++;
      if (!u || p.age > 140) { seedParticle(p); p.age = 0; continue; }
      const nx = p.x + u[0] * gain, ny = p.y + u[1] * gain;
      if (nx < 0 || ny < 0 || nx >= DOMAIN_M || ny >= DOMAIN_M) { seedParticle(p); p.age = 0; continue; }
      if (draw2d) { c.moveTo(p.x * k, p.y * k); c.lineTo(nx * k, ny * k); }
      p.x = nx; p.y = ny;
    }
    if (draw2d) c.stroke();
    if (view3d && S.view === '3d') view3d.pushParticles(particles, S.showParticles);
  }

  // ---------------------------------------------------------------- loop
  let lastFrame = 0, fpsAcc = 0, fpsN = 0, fpsT = performance.now();
  function frame(now) {
    requestAnimationFrame(frame);
    if (document.hidden) return;                       // pause when the tab is hidden
    if (now - lastFrame < 1000 / FPS_CAP - 2) return;  // frame cap
    lastFrame = now;
    if (S.running && now - S.lastTouch > IDLE_PAUSE_MS) { S.running = false; S.idlePaused = true; syncRun(); }
    if (solver && S.data && !S.loading && S.running) {
      if (solver.lost) { console.warn('[wind] GPU context lost, switching to CPU'); params.set('engine', 'cpu'); restart(); return; }
      const t0 = performance.now();
      solver.step(S.spf);
      S.steps += S.spf;
      const buf = solver.read();
      if (S.showField) { paintField(buf); if (view3d) view3d.markField(); }
      moveParticles(buf, S.spf);
      const work = performance.now() - t0;
      // keep the page light: aim for well under one frame of work
      if (work > 24 && S.spf > 2) S.spf = Math.max(2, Math.floor(S.spf * 0.8));
      else if (work < 10 && S.spf < (solver.kind === 'gpu' ? 36 : 8)) S.spf += 1;
      fpsN++;
    }
    if (view3d && S.view === '3d') { sync3d(); view3d.frame(); updateNorth(); }
    if (now - fpsT > 1000) { S.fps = fpsN * 1000 / (now - fpsT); fpsN = 0; fpsT = now; if (solver) updateStatus(); }
  }
  function updateStatus() {
    const flow = S.steps * U0 / solver.N;      // how many times the wind has crossed the square
    const pct = Math.min(100, Math.round(flow / 1.5 * 100));
    settleEl.textContent = pct < 100 ? `Settling 形成中 ${pct}%` : 'Settled 已穩定 (eddies keep moving, as real wind does)';
    settleEl.classList.toggle('ok', pct >= 100);
    statusEl.textContent = `${S.fps.toFixed(0)} fps · ${S.spf} steps per frame`;
  }
  function syncRun() {
    runBtn.textContent = S.running ? 'Pause 暫停' : 'Run 開始';
    runBtn.classList.toggle('primary', !S.running);
    pauseNote.hidden = !S.idlePaused;
  }
  function touch() { S.lastTouch = performance.now(); if (S.idlePaused) { S.idlePaused = false; S.running = true; syncRun(); } }
  root.addEventListener('pointerdown', touch, true);
  root.addEventListener('keydown', touch, true);
  document.addEventListener('visibilitychange', () => { if (!document.hidden) lastFrame = 0; });

  // ---------------------------------------------------------------- area + wind
  const cache = {};
  function prepBuildings(d) {
    for (const b of d.buildings) {
      let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
      for (const r of b.rings) for (const [x, y] of r) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      b.bb = [x0, y0, x1, y1];
    }
    return d;
  }
  async function getJSON(url, opt) {
    const r = await fetch(url, opt);
    if (!r.ok) throw new Error(`${r.status} from ${new URL(url, location.href).host}`);
    return r.json();
  }
  async function loadArea(id) {
    if (cache[id]) return cache[id];
    const d = (window.HKW_DATA && window.HKW_DATA[id]) || await getJSON(DATA_URL + id + '.json');
    return (cache[id] = prepBuildings(d));
  }

  // any site: LandsD buildings from the CSDI web feature service, by bounding box
  async function fetchBuildings(lat, lon, onCount) {
    const half = DOMAIN_M / 2 + FETCH_PAD_M, dlat = half / M_PER_DEG_LAT, dlon = half / mPerDegLon(lat);
    const bbox = [lat - dlat, lon - dlon, lat + dlat, lon + dlon].map((v) => v.toFixed(6)).join(',');
    let feats = [];
    for (let start = 0; start < 20000; start += 1000) {
      const url = `${WFS}?service=wfs&version=2.0.0&request=GetFeature&typeNames=csdi:Building&outputFormat=GEOJSON` +
        `&count=1000&startIndex=${start}&bbox=${bbox},urn:ogc:def:crs:EPSG::4326`;
      const d = await getJSON(url);
      const batch = d.features || [];
      feats = feats.concat(batch); onCount(feats.length);
      if (batch.length < 1000) break;
    }
    return feats;
  }
  // same rules as prep_data.py: active blocks with both levels, outer rings, local metres (x east, y south)
  function convertFeatures(feats, lat, lon) {
    const kx = mPerDegLon(lat), half = DOMAIN_M / 2, out = [];
    const seen = new Set();
    for (const f of feats) {
      const p = f.properties || {}, g = f.geometry || {};
      if (p.Status != null && p.Status !== 'Active') continue;
      if (p.BaseHeight == null || p.TopHeight == null) continue;
      if (seen.has(p.OBJECTID)) continue; seen.add(p.OBJECTID);
      const polys = g.type === 'Polygon' ? [g.coordinates] : (g.coordinates || []);
      const rings = [];
      for (const poly of polys) if (poly && poly[0]) rings.push(poly[0].map(([x, y]) => [+((x - lon) * kx + half).toFixed(2), +(-(y - lat) * M_PER_DEG_LAT + half).toFixed(2)]));
      if (!rings.length) continue;
      out.push({ id: p.OBJECTID, name: (p.BuildingNameEN || '').trim(), name_zh: (p.BuildingNameTC || '').trim(), type: p.BuildingBlockType || '',
        base: +(+p.BaseHeight).toFixed(1), top: +(+p.TopHeight).toFixed(1), h: +(p.TopHeight - p.BaseHeight).toFixed(1), storeys: p.Storeys, rings });
    }
    return out;
  }
  // PlanD Site Wind Availability Data: the 500 m cell under the centre, 200 m level, three seasons
  let gridP = null;
  async function fetchSiteWind(lat, lon) {
    gridP = gridP || getJSON(DATA_URL + 'pland_grid.json');
    const G = await gridP;
    const [E, N] = HK80.toGrid(lat, lon);
    const c = Math.floor((E - G.x0) / G.step), r = Math.floor((N - G.y0) / G.step);
    const id = c >= 0 && r >= 0 && c < G.cols && r < G.rows ? G.ids[r * G.cols + c] : '';
    const res = { height_m: 200, seasons: {} };
    if (!id) return res;
    res.grid_id = id; res.cell_hk1980 = [G.x0 + c * G.step, G.y0 + r * G.step, G.x0 + (c + 1) * G.step, G.y0 + (r + 1) * G.step];
    await Promise.all(['annual', 'summer', 'winter'].map(async (season) => {
      const url = G.url.replace('{season}', season).replace('{id}', id);
      try {
        const r2 = await fetch(url); if (!r2.ok) return;
        const lines = (await r2.text()).split(/\r?\n/);
        const dirs = lines[0].split('\t').slice(2, 18).map((v) => v.trim());
        const freq = lines[1].split('\t').slice(2, 18).map(Number);
        if (dirs.length === 16 && freq.every(isFinite)) res.seasons[season] = { dirs, freq, source: url };
      } catch (e) { console.info('[wind] PlanD table not loaded', season, e.message); }
    }));
    return res;
  }
  async function loadSite(lat, lon, label, onCount) {
    const key = lat.toFixed(5) + ',' + lon.toFixed(5);
    if (cache[key]) return cache[key];
    const [feats, wind] = await Promise.all([fetchBuildings(lat, lon, onCount), fetchSiteWind(lat, lon)]);
    const d = prepBuildings({ id: key, name_en: label || 'Picked on the map', name_zh: label ? '' : '地圖揀點', live: true,
      centre_wgs84: [lat, lon], size_m: DOMAIN_M, buildings: convertFeatures(feats, lat, lon), site_wind: wind,
      fetched: new Date().toISOString().slice(0, 10) });
    return (cache[key] = d);
  }

  function seasonFreq() {
    const s = S.data && S.data.site_wind && S.data.site_wind.seasons[S.season];
    if (!s) return null;
    const m = {}; s.dirs.forEach((d, i) => { m[d] = s.freq[i]; });
    return DIRS.map(([d]) => m[d] || 0);
  }
  function topDir(season) {
    const s = S.data && S.data.site_wind && S.data.site_wind.seasons[season]; if (!s) return 4;
    let bi = 0; s.freq.forEach((f, i) => { if (f > s.freq[bi]) bi = i; });
    return Math.max(0, DIRS.findIndex(([d]) => d === s.dirs[bi]));
  }

  const loadingEl = $('#hkwLoading');
  let loadToken = 0;
  function showLoading(html) { loadingEl.innerHTML = html; loadingEl.hidden = !html; }
  function setURL(q) {
    try {
      const u = new URL(location.href);
      ['area', 'lat', 'lon', 'q'].forEach((k) => u.searchParams.delete(k));
      for (const [k, v] of Object.entries(q)) if (v != null && v !== '') u.searchParams.set(k, v);
      history.replaceState(null, '', u);
    } catch (e) { /* file:// or sandboxed */ }
  }
  async function selectArea(id) {
    const tok = ++loadToken;
    S.areaId = id; S.loading = true;
    root.querySelectorAll('[data-area]').forEach((b) => b.classList.toggle('on', b.dataset.area === id));
    showLoading('Loading 載入中…');
    let d;
    try { d = await loadArea(id); } catch (e) { if (tok === loadToken) { S.loading = false; showLoading(''); $('#hkwAreaNote').textContent = 'Could not load this area: ' + e.message; } return; }
    if (tok !== loadToken) return;
    setURL({ area: id });
    if (d.centre_wgs84) showSquare(d.centre_wgs84[0], d.centre_wgs84[1], false);
    useData(d);
  }
  async function selectSite(lat, lon, label) {
    if (!(lat > HK_BOUNDS[0] && lat < HK_BOUNDS[2] && lon > HK_BOUNDS[1] && lon < HK_BOUNDS[3])) {
      $('#hkwAreaNote').textContent = 'That point is outside Hong Kong. Pick a place inside Hong Kong. 請揀香港境內嘅地點。';
      return;
    }
    const tok = ++loadToken;
    S.areaId = null; S.loading = true; touch();
    root.querySelectorAll('[data-area]').forEach((b) => b.classList.remove('on'));
    showSquare(lat, lon, true);
    showLoading('Fetching buildings from Lands Department… 攞緊屋宇資料…');
    let d;
    try {
      d = await loadSite(lat, lon, label, (n) => { if (tok === loadToken) showLoading(`Fetching buildings from Lands Department… ${n} so far 攞緊屋宇資料…`); });
    } catch (e) {
      if (tok !== loadToken) return;
      console.info('[wind] site load failed:', e.message);
      showLoading(`<div>The Lands Department data service did not answer (${escapeHTML(e.message)}). 政府數據服務冇回應。<br><button type="button" class="hkw-btn" id="hkwRetry">Try again 再試</button></div>`);
      $('#hkwRetry').addEventListener('click', () => selectSite(lat, lon, label));
      return;
    }
    if (tok !== loadToken) return;
    setURL({ lat: lat.toFixed(5), lon: lon.toFixed(5), q: label || '' });
    useData(d);
  }
  function escapeHTML(t) { return String(t).replace(/[&<>"']/g, (c) => ({ '&': '&amp;', '<': '&lt;', '>': '&gt;', '"': '&quot;', "'": '&#39;' }[c])); }
  function useData(d) {
    S.loading = false; showLoading('');
    S.data = d; buildings = d.buildings; S.user = []; S.removed = new Set(); S.sel = null; syncSel();
    need3d.site = true; need3d.user = true;
    const nb = buildings.filter(blocks).length;
    const sw = d.site_wind || {};
    const windNote = sw.grid_id ? `Wind table: PlanD grid cell ${escapeHTML(sw.grid_id)}.` : 'No PlanD wind table for this point (sea or outside the grid): all 16 directions are shown as equally likely. 呢點冇 PlanD 風資料。';
    if (d.live) {
      const [la, lo] = d.centre_wgs84;
      $('#hkwAreaNote').innerHTML = `<b>${escapeHTML(d.name_en)}${d.name_zh ? ' ' + escapeHTML(d.name_zh) : ''}</b> · ${la.toFixed(5)} N, ${lo.toFixed(5)} E · ${d.size_m} m x ${d.size_m} m · ` +
        `${nb} building blocks stand in the 2 m slice (of ${buildings.length} in LandsD data, fetched live ${d.fetched}). ${windNote}` +
        (buildings.length === 0 ? '<br>No buildings here: open ground or sea. 呢度冇建築物。' : '');
    } else {
      $('#hkwAreaNote').innerHTML = `<b>${d.name_en} ${d.name_zh}</b> · ${d.size_m} m x ${d.size_m} m · ${nb} building blocks stand in the 2 m slice (of ${buildings.length} in LandsD data).<br>${d.note_en}<br>${d.note_zh}`;
    }
    const want = params.get('dir');
    const wi = want && S.dirIdx === null ? DIRS.findIndex(([e]) => e === want.toUpperCase()) : -1;
    S.dirIdx = wi >= 0 ? wi : topDir('annual');       // each site starts from its own most common wind
    drawRose(); restartFlow();
  }

  // ---------------------------------------------------------------- base map (Leaflet + LandsD tiles)
  let map = null, square = null;
  function initMap() {
    const el = $('#hkwMap');
    if (!window.L) { el.innerHTML = '<p class="small" style="padding:12px">The map could not load (offline?). Search still works. 地圖載入唔到，搜尋仍然用得。</p>'; return; }
    map = L.map(el, { minZoom: 10, maxZoom: 19, maxBounds: [[21.9, 113.6], [22.8, 114.7]], scrollWheelZoom: false }).setView([22.33, 114.17], 11);
    L.tileLayer('https://mapapi.geodata.gov.hk/gs/api/v1.0.0/xyz/basemap/wgs84/{z}/{x}/{y}.png', {
      maxZoom: 19, attribution: '<a href="https://api.portal.hkmapservice.gov.hk/disclaimer" target="_blank" rel="noopener">&copy; Map from Lands Department</a>',
    }).addTo(map);
    L.tileLayer('https://mapapi.geodata.gov.hk/gs/api/v1.0.0/xyz/label/hk/en/wgs84/{z}/{x}/{y}.png', { maxZoom: 19 }).addTo(map);
    const Logo = L.Control.extend({ options: { position: 'bottomleft' }, onAdd() {
      const d = L.DomUtil.create('div', 'hkw-landsd');
      d.innerHTML = '<img src="https://api.hkmapservice.gov.hk/mapapi/landsdlogo.jpg" alt="Lands Department">';
      return d;
    } });
    new Logo().addTo(map);
    map.on('focus', () => map.scrollWheelZoom.enable());
    map.on('click', (e) => { qEl.value = ''; selectSite(e.latlng.lat, e.latlng.lng, null); });
  }
  function showSquare(lat, lon, fly) {
    if (!map) return;
    const half = DOMAIN_M / 2, dlat = half / M_PER_DEG_LAT, dlon = half / mPerDegLon(lat);
    const b = [[lat - dlat, lon - dlon], [lat + dlat, lon + dlon]];
    if (square) square.setBounds(b);
    else square = L.rectangle(b, { color: '#c46240', weight: 2, fillOpacity: 0.08, interactive: false }).addTo(map);
    if (fly) map.flyToBounds(b, { padding: [30, 30], maxZoom: 16, duration: 0.8 });
    else map.fitBounds(b, { padding: [30, 30], maxZoom: 15 });
  }

  // ---------------------------------------------------------------- search (LandsD place names + ALS addresses)
  const qEl = $('#hkwQ'), suggEl = $('#hkwSugg');
  let sTimer = 0, sAbort = null, sItems = [], sIdx = -1;
  const titleCase = (t) => (t || '').toLowerCase().replace(/\b([a-z])/g, (m) => m.toUpperCase());
  async function placeSearch(v, signal) {
    const arr = await getJSON(PLACE_SEARCH + encodeURIComponent(v), { signal });
    // district and estate names (no street address) before shops and kiosks
    const rows = (arr || []).map((r, i) => [r.addressEN ? 1 : 0, i, r]).sort((x, y) => x[0] - y[0] || x[1] - y[1]).map((x) => x[2]);
    return rows.slice(0, 7).map((r) => {
      const [lat, lon] = HK80.toWGS(+r.x, +r.y);
      return { kind: 'Place 地方', title: `${titleCase(r.nameEN)} ${r.nameZH || ''}`.trim(), sub: [titleCase(r.addressEN), r.districtEN].filter(Boolean).join(' · '), lat, lon };
    });
  }
  async function addrSearch(v, signal) {
    const d = await getJSON(ALS + encodeURIComponent(v), { signal, headers: { Accept: 'application/json' } });
    return (d.SuggestedAddress || []).slice(0, 4).map((s) => {
      const pa = s.Address.PremisesAddress, e = pa.EngPremisesAddress || {}, c = pa.ChiPremisesAddress || {}, g = pa.GeospatialInformation || {};
      const st = e.EngStreet ? [e.EngStreet.BuildingNoFrom, titleCase(e.EngStreet.StreetName)].filter(Boolean).join(' ') : '';
      const en = [titleCase(e.BuildingName), st, e.EngEstate && titleCase(e.EngEstate.EstateName)].filter(Boolean).join(', ');
      const zh = [c.ChiStreet && (c.ChiStreet.StreetName || '') + (c.ChiStreet.BuildingNoFrom ? c.ChiStreet.BuildingNoFrom + '號' : ''), c.BuildingName].filter(Boolean).join(' ');
      const dist = e.EngDistrict && e.EngDistrict.DcDistrict;
      return { kind: 'Address 地址', title: en || zh, sub: [zh, titleCase(dist)].filter(Boolean).join(' · '), lat: +g.Latitude, lon: +g.Longitude };
    }).filter((x) => isFinite(x.lat) && isFinite(x.lon));
  }
  function renderSugg(items, msg) {
    sItems = items; sIdx = -1;
    suggEl.innerHTML = '';
    if (msg) { const li = document.createElement('li'); li.className = 'msg'; li.textContent = msg; suggEl.appendChild(li); }
    items.forEach((it, i) => {
      const li = document.createElement('li'); li.setAttribute('role', 'option'); li.id = 'hkwS' + i;
      const k = document.createElement('span'); k.className = 'k'; k.textContent = it.kind;
      const s = document.createElement('span'); s.className = 's'; s.textContent = it.sub;
      li.append(k, document.createTextNode(it.title), s);
      li.addEventListener('mousedown', (e) => { e.preventDefault(); pickSugg(i); });
      suggEl.appendChild(li);
    });
    suggEl.hidden = !msg && !items.length;
    qEl.setAttribute('aria-expanded', String(!suggEl.hidden));
  }
  function hideSugg() { suggEl.hidden = true; qEl.setAttribute('aria-expanded', 'false'); }
  function pickSugg(i) {
    const it = sItems[i]; if (!it) return;
    qEl.value = it.title; hideSugg(); qEl.blur();
    selectSite(it.lat, it.lon, it.title);
  }
  async function runSearch(v) {
    if (sAbort) sAbort.abort();
    const ac = sAbort = new AbortController();
    renderSugg([], 'Searching… 搜尋中…');
    const [pl, ad] = await Promise.allSettled([placeSearch(v, ac.signal), addrSearch(v, ac.signal)]);
    if (ac.signal.aborted) return;
    const items = [];
    for (const it of [...(pl.value || []), ...(ad.value || [])]) {
      if (!items.some((o) => Math.abs(o.lat - it.lat) < 0.0004 && Math.abs(o.lon - it.lon) < 0.0004)) items.push(it);
    }
    if (items.length) renderSugg(items.slice(0, 9));
    else renderSugg([], pl.status === 'rejected' && ad.status === 'rejected' ? 'The search services did not answer. Try again, or click the map. 搜尋服務冇回應。' : 'No match. Try another name, or click the map. 搵唔到，試吓第二個名或者撳地圖。');
  }
  qEl.addEventListener('input', () => {
    clearTimeout(sTimer);
    const v = qEl.value.trim();
    if (v.length < 2) { hideSugg(); return; }
    sTimer = setTimeout(() => runSearch(v), 380);     // be gentle with the government APIs
  });
  qEl.addEventListener('keydown', (e) => {
    if (e.key === 'ArrowDown' || e.key === 'ArrowUp') {
      if (!sItems.length) return; e.preventDefault();
      sIdx = (sIdx + (e.key === 'ArrowDown' ? 1 : -1) + sItems.length) % sItems.length;
      suggEl.querySelectorAll('[role=option]').forEach((li, i) => li.setAttribute('aria-selected', String(i === sIdx)));
      qEl.setAttribute('aria-activedescendant', 'hkwS' + sIdx);
    } else if (e.key === 'Enter') {
      e.preventDefault();
      if (sItems.length) pickSugg(Math.max(0, sIdx));
      else if (qEl.value.trim().length >= 2) { clearTimeout(sTimer); runSearch(qEl.value.trim()); }
    } else if (e.key === 'Escape') hideSugg();
  });
  qEl.addEventListener('blur', () => setTimeout(hideSugg, 150));

  // ---------------------------------------------------------------- 3D view (three.js, loaded on first use)
  const stateFor3d = () => ({ buildings, removed: S.removed, user: S.user, sel: S.sel, blocks });
  let lastRemovedKey = '';
  function sync3d() {
    const key = [...S.removed].join(',');
    if (need3d.site || key !== lastRemovedKey) { view3d.rebuild(stateFor3d(), need3d.site); need3d.site = false; lastRemovedKey = key; }
    if (need3d.user) { view3d.rebuildUser(stateFor3d()); need3d.user = false; }
  }
  function ensure3d() {
    if (!view3dLoading) {
      view3dLoading = import(new URL(STATIC + 'wind3d.js' + (window.HKW_VER ? '?v=' + window.HKW_VER : ''), document.baseURI).href).then((m) => {
        view3d = m.create({ container: $('#hkw3d'), fieldCanvas: cvField, domain: DOMAIN_M, onPick: pick3d });
        view3d.setTheme(isDark());
      });
      view3dLoading.catch(() => { view3dLoading = null; });
    }
    return view3dLoading;
  }
  function updateNorth() {
    const n = $('#hkwNorth');
    n.style.transform = S.view === '3d' && view3d ? `rotate(${view3d.azimuth()}rad)` : '';
  }
  const VIEW_HELP = {
    '2d': 'Plan view, north up. 平面，北向上。',
    '3d': 'Drag to orbit, right-drag or two fingers to pan, scroll or pinch to zoom. 拖曳轉動，雙指移動或縮放。',
  };
  async function setView(v) {
    touch();
    if (v === '3d' && !view3d) {
      $('#hkwViewHelp').textContent = 'Loading 3D… 載入 3D 中…';
      try { await ensure3d(); } catch (e) {
        console.warn('[wind] 3D could not start:', e.message);
        $('#hkwViewHelp').textContent = '3D could not start here (' + e.message + '). The plan view still works. 3D 開唔到。';
        return;
      }
    }
    S.view = v;
    root.querySelectorAll('[data-view]').forEach((b) => b.classList.toggle('on', b.dataset.view === v));
    stage.classList.toggle('is3d', v === '3d');
    $('#hkw3d').hidden = v !== '3d';
    $('#hkwResetView').hidden = v !== '3d';
    $('#hkwViewHelp').textContent = VIEW_HELP[v];
    tip.hidden = true;
    if (v === '3d') {
      view3d.resize(S.size, S.size);
      need3d.user = true;
      if (S.dirIdx !== null) view3d.setWind(S.dirIdx);
      particles.forEach((p) => { p.fresh = true; });
      view3d.showField(S.showField);
      if (solver && S.steps > 0 && S.showField) { paintField(solver.read()); view3d.markField(); }
    } else {
      ctxT.setTransform(1, 0, 0, 1, 0, 0); ctxT.clearRect(0, 0, cvTrail.width, cvTrail.height);
      drawBuildings();
    }
    updateNorth();
  }
  function placeTip(ev, html) {
    if (!html) { tip.hidden = true; return; }
    tip.innerHTML = html; tip.hidden = false;
    const r = stage.getBoundingClientRect();
    let lx = ev.clientX - r.left + 14, ly = ev.clientY - r.top + 14;
    if (lx > S.size - 220) lx = ev.clientX - r.left - 220;
    if (ly > S.size - 60) ly = ev.clientY - r.top - 60;
    tip.style.left = lx + 'px'; tip.style.top = ly + 'px';
  }
  function buildingTip(b) {
    const nm = [b.name, b.name_zh].filter(Boolean).join(' ') || 'Unnamed block 未命名';
    const st = S.removed.has(b.id) ? 'taken away in this sketch' : (blocks(b) ? 'blocks the 2 m slice' : (b.type === 'Open-sided Structure' ? 'open-sided, wind passes under' : 'under 2.5 m, ignored'));
    return `<b>${escapeHTML(nm)}</b><br>${escapeHTML(b.type || 'Block')} · ${b.h} m tall${b.storeys ? `, ${b.storeys} storeys` : ''} · ${st}`;
  }
  function pick3d(hit) {
    touch();
    if (S.mode === 'add') {
      if (hit.kind === 'user') { S.sel = hit.u; }
      else if (hit.kind === 'ground') {
        const h = +$('#hkwNewH').value || 60;
        const u = { id: nextUserId++, cx: Math.round(hit.x), cy: Math.round(hit.y), w: 40, d: 30, rot: 0, h };
        S.user.push(u); S.sel = u; rebuildMask();
      } else { placeTip(hit.ev, 'Click open ground to place a 40 x 30 m block. 撳空地放一座 40 x 30 米嘅樓。'); return; }
    } else if (S.mode === 'remove') {
      if (hit.kind === 'user') { S.user = S.user.filter((v) => v !== hit.u); if (S.sel === hit.u) S.sel = null; }
      else if (hit.kind === 'building') S.removed.add(hit.b.id);
      else if (hit.kind === 'ground') {
        const b = hitBuilding(hit.x, hit.y);       // a taken-away block has no solid: click its outline to bring it back
        if (b && S.removed.has(b.id)) S.removed.delete(b.id);
        else return;
      }
      rebuildMask();
    } else {
      if (hit.kind === 'user') { S.sel = hit.u; placeTip(hit.ev, `<b>Your block 你加嘅樓</b><br>${hit.u.w} m x ${hit.u.d} m, ${hit.u.h} m tall`); }
      else if (hit.kind === 'building') placeTip(hit.ev, buildingTip(hit.b));
      else { S.sel = null; tip.hidden = true; }
    }
    syncSel(); drawBuildings();
  }
  root.addEventListener('click', (ev) => {
    const t = ev.target.closest('[data-view]');
    if (t && root.contains(t)) setView(t.dataset.view);
  });
  $('#hkwResetView').addEventListener('click', () => { if (view3d) view3d.resetView(); });

  function restartFlow() {
    if (!solver) newSolver();
    rebuildMask();
    solver.init(inflow(S.dirIdx));
    S.steps = 0; resetParticles(); drawBuildings(); updateWindText();
    if (view3d) view3d.setWind(S.dirIdx);
    if (!S.running && !S.idlePaused) { /* stay paused */ }
  }
  function restart() { S.field = null; newSolver(); restartFlow(); }
  function setDir(i) { S.dirIdx = i; touch(); drawRose(); restartFlow(); }
  function updateWindText() {
    const [en, zh] = DIRS[S.dirIdx];
    const fr = seasonFreq();
    const pct = fr ? Math.round(fr[S.dirIdx] * 100) : null;
    $('#hkwWindText').innerHTML = `Wind from <b>${en} ${zh}風</b>` +
      (pct !== null ? ` · about ${pct}% of the time in ${SEASONS[S.season]} at 200 m (PlanD simulated data for this 500 m grid cell)` : '');
    const tops = ['winter', 'summer'].map((s) => { const i = topDir(s); return `${SEASONS[s]}: <button type="button" class="hkw-link" data-dir="${i}">${DIRS[i][0]} ${DIRS[i][1]}</button>`; });
    $('#hkwTops').innerHTML = 'Most common here 呢區最常見: ' + tops.join(' · ');
  }

  // wind rose / direction picker (SVG)
  function drawRose() {
    const svg = $('#hkwRose'); if (!svg) return;
    const fr = seasonFreq() || DIRS.map(() => 1 / 16);
    const mx = Math.max(...fr, 0.01);
    const cx = 110, cy = 110, r0 = 20;
    let h = `<circle cx="${cx}" cy="${cy}" r="86" class="ring"/><circle cx="${cx}" cy="${cy}" r="${r0}" class="ring"/>`;
    DIRS.forEach(([en, zh], i) => {
      const a0 = (i * 22.5 - 9) * Math.PI / 180, a1 = (i * 22.5 + 9) * Math.PI / 180;
      const r1 = r0 + 64 * (fr[i] / mx);
      const P = (a, r) => `${(cx + Math.sin(a) * r).toFixed(1)},${(cy - Math.cos(a) * r).toFixed(1)}`;
      const on = i === S.dirIdx ? ' on' : '';
      h += `<path class="petal${on}" d="M${P(a0, r0)} L${P(a0, r1)} A${r1} ${r1} 0 0 1 ${P(a1, r1)} L${P(a1, r0)} A${r0} ${r0} 0 0 0 ${P(a0, r0)}Z"/>`;
      h += `<path class="hit" data-dir="${i}" d="M${P(a0 - 0.02, r0)} L${P(a0 - 0.02, 100)} A100 100 0 0 1 ${P(a1 + 0.02, 100)} L${P(a1 + 0.02, r0)}Z"><title>${en} ${zh}風 ${Math.round(fr[i] * 100)}%</title></path>`;
      if (i % 2 === 0) h += `<text x="${P(i * 22.5 * Math.PI / 180, 100).split(',')[0]}" y="${P(i * 22.5 * Math.PI / 180, 100).split(',')[1]}" class="lbl${on}">${en}</text>`;
    });
    const b = S.dirIdx * 22.5 * Math.PI / 180;   // arrow points downwind
    const tip = [cx - Math.sin(b) * 14, cy + Math.cos(b) * 14], tail = [cx + Math.sin(b) * 14, cy - Math.cos(b) * 14];
    h += `<line x1="${tail[0]}" y1="${tail[1]}" x2="${tip[0]}" y2="${tip[1]}" class="arrow" marker-end="url(#hkwArrowHead)"/>`;
    svg.innerHTML = `<defs><marker id="hkwArrowHead" viewBox="0 0 10 10" refX="6" refY="5" markerWidth="5" markerHeight="5" orient="auto"><path d="M0,0 L10,5 L0,10z" class="arrowhead"/></marker></defs>` + h;
    root.querySelectorAll('[data-season]').forEach((x) => x.classList.toggle('on', x.dataset.season === S.season));
  }

  // ---------------------------------------------------------------- editing
  function toM(ev) {
    const r = cvBld.getBoundingClientRect();
    return [(ev.clientX - r.left) / r.width * DOMAIN_M, (ev.clientY - r.top) / r.height * DOMAIN_M];
  }
  let drag = null, maskTimer = 0;
  function scheduleMask() {
    if (maskTimer) return;
    maskTimer = setTimeout(() => { maskTimer = 0; rebuildMask(); }, 90);
  }
  function syncSel() {
    const box = $('#hkwSel');
    if (!S.sel) { box.hidden = true; return; }
    box.hidden = false;
    $('#hkwSelH').value = S.sel.h; $('#hkwSelW').value = S.sel.w; $('#hkwSelD').value = S.sel.d; $('#hkwSelR').value = S.sel.rot;
    $('#hkwSelRv').textContent = S.sel.rot + '°';
  }
  function setMode(m) {
    S.mode = m;
    root.querySelectorAll('[data-mode]').forEach((b) => b.classList.toggle('on', b.dataset.mode === m));
    cvBld.style.cursor = m === 'add' ? 'crosshair' : (m === 'remove' ? 'not-allowed' : 'default');
    $('#hkwModeHelp').textContent = {
      look: 'Hover a building for its name and height, or a street for how windy it is here.',
      add: 'Drag on an open space to draw a new block. Drag an orange block to move it. 拖曳畫新樓，拖橙色樓可以搬。',
      remove: 'Click an orange block to delete it, or an existing building to take it away (what if it was not there?). Click again to bring it back.',
    }[m];
  }
  cvBld.addEventListener('pointerdown', (ev) => {
    const [x, y] = toM(ev);
    if (S.mode === 'add') {
      const u = hitUser(x, y);
      if (u) { S.sel = u; drag = { kind: 'move', u, ox: x - u.cx, oy: y - u.cy }; }
      else { drag = { kind: 'draw', x0: x, y0: y }; S.draft = { cx: x, cy: y, w: 1, d: 1, rot: 0, h: +$('#hkwNewH').value || 60 }; }
      cvBld.setPointerCapture(ev.pointerId); syncSel(); drawBuildings();
    } else if (S.mode === 'remove') {
      const u = hitUser(x, y);
      if (u) { S.user = S.user.filter((v) => v !== u); if (S.sel === u) S.sel = null; }
      else {
        const b = hitBuilding(x, y);
        if (b) { if (S.removed.has(b.id)) S.removed.delete(b.id); else S.removed.add(b.id); }
      }
      syncSel(); drawBuildings(); rebuildMask();
    } else {
      const u = hitUser(x, y); S.sel = u; syncSel(); drawBuildings();
    }
  });
  cvBld.addEventListener('pointermove', (ev) => {
    const [x, y] = toM(ev);
    if (drag && drag.kind === 'move') {
      drag.u.cx = Math.min(DOMAIN_M, Math.max(0, x - drag.ox)); drag.u.cy = Math.min(DOMAIN_M, Math.max(0, y - drag.oy));
      drawBuildings(); scheduleMask(); return;
    }
    if (drag && drag.kind === 'draw') {
      const d = S.draft; d.cx = (drag.x0 + x) / 2; d.cy = (drag.y0 + y) / 2; d.w = Math.abs(x - drag.x0); d.d = Math.abs(y - drag.y0);
      drawBuildings(); return;
    }
    showTip(ev, x, y);
  });
  const endDrag = () => {
    if (drag && drag.kind === 'draw') {
      const d = S.draft; S.draft = null;
      if (d.w >= 6 && d.d >= 6) {
        const u = { id: nextUserId++, cx: Math.round(d.cx), cy: Math.round(d.cy), w: Math.round(d.w), d: Math.round(d.d), rot: 0, h: d.h };
        S.user.push(u); S.sel = u; rebuildMask();
      }
      syncSel(); drawBuildings();
    } else if (drag && drag.kind === 'move') { rebuildMask(); }
    drag = null;
  };
  cvBld.addEventListener('pointerup', endDrag);
  cvBld.addEventListener('pointercancel', endDrag);
  cvBld.addEventListener('pointerleave', () => { tip.hidden = true; });

  function showTip(ev, x, y) {
    if (S.mode !== 'look' && S.mode !== 'remove') { tip.hidden = true; return; }
    let html = '';
    const u = hitUser(x, y);
    const b = u ? null : hitBuilding(x, y);
    if (u) html = `<b>Your block 你加嘅樓</b><br>${u.w} m x ${u.d} m, ${u.h} m tall`;
    else if (b) {
      const nm = [b.name, b.name_zh].filter(Boolean).join(' ') || 'Unnamed block 未命名';
      const st = S.removed.has(b.id) ? 'taken away in this sketch' : (blocks(b) ? 'blocks the 2 m slice' : (b.type === 'Open-sided Structure' ? 'open-sided, wind passes under' : 'under 2.5 m, ignored'));
      html = `<b>${nm}</b><br>${b.type || 'Block'} · ${b.h} m tall${b.storeys ? `, ${b.storeys} storeys` : ''} · ${st}`;
    } else if (solver && S.steps > 0) {
      const v = sampleU(solver.read(), x, y);
      if (v) {
        const r = Math.hypot(v[0], v[1]) / U0;
        const w = r < 0.35 ? 'Calm 好靜' : r < 0.75 ? 'Calmer than the incoming wind 比來風弱' : r < 1.1 ? 'About the same as the incoming wind 同來風差唔多' : 'Faster than the incoming wind 比來風強';
        html = `<b>${w}</b>${S.steps * U0 / solver.N < 1.5 ? '<br><span class="small">still settling</span>' : ''}`;
      }
    }
    if (!html) { tip.hidden = true; return; }
    tip.innerHTML = html; tip.hidden = false;
    const r = stage.getBoundingClientRect();
    let lx = ev.clientX - r.left + 14, ly = ev.clientY - r.top + 14;
    if (lx > S.size - 220) lx = ev.clientX - r.left - 220;
    if (ly > S.size - 60) ly = ev.clientY - r.top - 60;
    tip.style.left = lx + 'px'; tip.style.top = ly + 'px';
  }

  // ---------------------------------------------------------------- controls
  root.addEventListener('click', (ev) => {
    const t = ev.target.closest('[data-area],[data-dir],[data-season],[data-mode]');
    if (!t || !root.contains(t)) return;
    if (t.dataset.area) selectArea(t.dataset.area);
    else if (t.dataset.dir !== undefined) setDir(+t.dataset.dir);
    else if (t.dataset.season) { S.season = t.dataset.season; drawRose(); updateWindText(); }
    else if (t.dataset.mode) setMode(t.dataset.mode);
  });
  runBtn.addEventListener('click', () => { S.running = !S.running; S.idlePaused = false; touch(); syncRun(); });
  $('#hkwReset').addEventListener('click', () => { touch(); restartFlow(); });
  $('#hkwUndo').addEventListener('click', () => { S.user = []; S.removed = new Set(); S.sel = null; syncSel(); drawBuildings(); rebuildMask(); });
  $('#hkwShowField').addEventListener('change', (e) => { S.showField = e.target.checked; cvField.style.visibility = S.showField ? '' : 'hidden'; if (view3d) view3d.showField(S.showField); });
  $('#hkwShowParticles').addEventListener('change', (e) => { S.showParticles = e.target.checked; });
  $('#hkwQuality').value = S.quality;
  $('#hkwQuality').addEventListener('change', (e) => { S.quality = e.target.value; restart(); });
  const selInput = (id, key, after) => $(id).addEventListener('input', (e) => {
    if (!S.sel) return; const v = +e.target.value; if (!isFinite(v)) return;
    S.sel[key] = v; if (after) after(); drawBuildings(); scheduleMask();
  });
  selInput('#hkwSelH', 'h'); selInput('#hkwSelW', 'w'); selInput('#hkwSelD', 'd');
  selInput('#hkwSelR', 'rot', () => { $('#hkwSelRv').textContent = S.sel.rot + '°'; });
  $('#hkwSelDel').addEventListener('click', () => { S.user = S.user.filter((v) => v !== S.sel); S.sel = null; syncSel(); drawBuildings(); rebuildMask(); });
  window.addEventListener('resize', () => { resize(); resetParticles(); });
  new MutationObserver(() => { drawBuildings(); if (view3d) view3d.setTheme(isDark()); }).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // test hook: add a block programmatically (used by the automated check)
  window.HKW.addBlock = (cx, cy, w, d, h, rot) => {
    const u = { id: nextUserId++, cx, cy, w, d, h: h || 60, rot: rot || 0 };
    S.user.push(u); S.sel = u; syncSel(); drawBuildings(); rebuildMask(); return u.id;
  };
  window.HKW.clearBlocks = () => { S.user = []; S.sel = null; syncSel(); drawBuildings(); rebuildMask(); };
  window.HKW.selectSite = selectSite; window.HKW.selectArea = selectArea; window.HKW.setView = setView;

  // ---------------------------------------------------------------- start
  setMode('look'); syncRun(); resize(); newSolver();
  initMap();
  $('#hkwViewHelp').textContent = VIEW_HELP['2d'];
  const plat = parseFloat(params.get('lat')), plon = parseFloat(params.get('lon'));
  if (isFinite(plat) && isFinite(plon)) { if (params.get('q')) qEl.value = params.get('q'); selectSite(plat, plon, params.get('q')); }
  else selectArea(S.areaId);
  if (params.get('view') === '3d') setView('3d');
  requestAnimationFrame(frame);
})();
