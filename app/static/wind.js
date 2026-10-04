/* HK Wind Sketch v3: a qualitative wind sketch for Hong Kong, in 3D.
 *
 * Main engine: a 3D lattice Boltzmann solver (D3Q19, regularised BGK +
 * Smagorinsky, see solver3d.js) on the graphics card, over a 600 m square of
 * real terrain (LandsD digital terrain model) and LandsD buildings, with a
 * power-law wind profile coming in from the chosen direction. Fallback: the v2
 * 2D slice about 2 m above the ground (D2Q9, WebGL2 or plain JavaScript).
 * Everything runs in the browser: no server, no upload.
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
  // 3D engine: cells across the square, and a cap on the total cell count per quality
  const Q3 = { light: 96, normal: 128, fine: 160 };
  const BUDGET3 = { light: 0.6e6, normal: 1.15e6, fine: 1.75e6 };
  const UREF = 0.1;                // lattice speed of the incoming wind at the reference height
  const ZREF_M = 200;              // reference height = PlanD's 200 m level
  const ALPHA = 0.25;              // power-law profile exponent (assumed, built-up terrain upwind)
  const NU3 = 0.005, CS3 = 0.17;   // lattice viscosity and Smagorinsky constant in 3D
  const READ_MS = 220;             // how often the 3D field is copied back for drawing
  const BATCH_MS = 18;             // GPU time per frame the 3D solver may use (of about 33 ms)
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
  const STOPS = [   // speed relative to the incoming wind at the same height -> colour (blues: pale = calm, deep = fast)
    [0.00, [246, 249, 252]], [0.30, [205, 224, 241]], [0.60, [140, 187, 223]],
    [0.90, [70, 140, 199]], [1.15, [28, 92, 168]], [1.40, [12, 50, 122]], [1.70, [52, 24, 96]],
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

  const isPhone = () => (window.matchMedia && matchMedia('(pointer: coarse)').matches) && Math.min(screen.width, screen.height) < 700;
  const S = {
    data: null, areaId: params.get('lat') ? null : (params.get('area') || 'kwun_tong'), dirIdx: null, season: 'annual',
    view: '2d', loading: false, siteLabel: '',
    quality: params.get('quality') || (isPhone() ? 'light' : 'normal'), running: true, idlePaused: false,
    // engine: '3d' (3D lattice Boltzmann) or '2d' (v2 slice). Phones start on the lighter 2D slice.
    engine: ['2d', 'cpu'].includes(params.get('engine')) ? '2d' : (params.get('engine') === '3d' || !isPhone() ? '3d' : '2d'),
    showField: true, showParticles: true, showLines: true, showSection: false, sectionOff: 0,
    mode: 'look', user: [], removed: new Set(), sel: null, selV: null,
    steps: 0, lastTouch: performance.now(), spf: 20, size: 600, dpr: 1, field: null, fps: 0,
    terrain: null, engineNote: '',
  };
  let solver = null, maskArr = null, particles = [], buildings = [], nextUserId = 1;
  let view3d = null, view3dLoading = null, need3d = { site: true, user: true, terrain: true };
  // 3D engine state
  let s3 = null, s3Key = '', s3Failed = false, g3 = null, field3 = null, field3At = 0, readBusy = false, flowToken = 0;
  let parts3 = [], lines3At = 0, s3Creating = null;
  const slice = { buf: null };
  window.HKW = { state: S, get solver() { return solver; }, get view3d() { return view3d; }, get s3() { return s3; }, get g3() { return g3; }, get field3() { return field3; }, HK80 };

  function forceCPU() { return params.get('engine') === 'cpu'; }
  function engineText() {
    if (solver && solver.kind === '3d') {
      const g = g3;
      return `3D wind on the graphics card (${s3.label}) · ${g.nx} x ${g.ny} x ${g.nz} cells, ${g.dx.toFixed(1)} m each, up to ${Math.round(g.nz * g.dx)} m above the lowest ground` +
        (g.coarsened ? '. Steep or tall site: cells made bigger to fit the air column.' : '') +
        (g.clipped ? '. Some towers or hills here are taller than that: the model cuts them off at its lid.' : '');
    }
    const why = S.engineNote ? ` ${S.engineNote}` : '';
    return (solver.kind === 'gpu'
      ? `2D slice on the graphics card (WebGL2) · ${solver.N} x ${solver.N} cells, ${(DOMAIN_M / solver.N).toFixed(1)} m each`
      : `2D slice on the processor (no WebGL2 float support) · ${solver.N} x ${solver.N} cells, ${(DOMAIN_M / solver.N).toFixed(1)} m each`) + why;
  }
  function syncQualityLabels() {
    const q = $('#hkwQuality'); if (!q) return;
    const lab = S.engine === '3d'
      ? { light: 'Light: 96 x 96 cells across (6.3 m)', normal: 'Normal: 128 x 128 cells across (4.7 m)', fine: 'Fine: 160 x 160 cells across (3.8 m)' }
      : { light: 'Light: 240 x 240 cells (2.5 m)', normal: 'Normal: 300 x 300 cells (2 m)', fine: 'Fine: 400 x 400 cells (1.5 m)' };
    for (const o of q.options) o.textContent = lab[o.value];
    root.querySelectorAll('[data-wengine]').forEach((b) => b.classList.toggle('on', b.dataset.wengine === S.engine));
  }

  // the 2D (v2) solver; the 3D one is made per site in prepare3d()
  function newSolver() {
    if (solver && solver.kind !== '3d') solver.dispose();
    const N = QUALITY[S.quality] || 300;
    solver = forceCPU() ? null : makeGPU(N);
    if (!solver) solver = makeCPU(CPU_N);
    S.spf = solver.kind === 'gpu' ? 20 : 3;
    engineEl.textContent = engineText();
    $('#hkwQuality').disabled = solver.kind === 'cpu';
  }

  // ---------------------------------------------------------------- geometry
  // user buildings are free polygons {id, poly: [[x, y], ...] (local metres, not closed), h}
  const userPoly = (u) => u.poly;
  const blocks = (b) => b.h >= BLOCK_MIN_H && b.type !== 'Open-sided Structure';
  function polyArea(P) { let a = 0; for (let i = 0, j = P.length - 1; i < P.length; j = i++) a += (P[j][0] + P[i][0]) * (P[j][1] - P[i][1]); return a / 2; }
  function polyCentroid(P) {
    let A = 0, cx = 0, cy = 0;
    for (let i = 0, j = P.length - 1; i < P.length; j = i++) { const f = P[j][0] * P[i][1] - P[i][0] * P[j][1]; A += f; cx += (P[j][0] + P[i][0]) * f; cy += (P[j][1] + P[i][1]) * f; }
    if (Math.abs(A) < 1e-9) return P[0].slice();
    return [cx / (3 * A), cy / (3 * A)];
  }
  function segX(a, b, c, d) {
    const o = (p, q, r) => Math.sign((q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0]));
    return o(a, b, c) * o(a, b, d) < 0 && o(c, d, a) * o(c, d, b) < 0;
  }
  function polySelfX(P) {
    const n = P.length;
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) {
      if (j === i + 1 || (i === 0 && j === n - 1)) continue;
      if (segX(P[i], P[(i + 1) % n], P[j], P[(j + 1) % n])) return true;
    }
    return false;
  }
  function pathRing(ctx, ring, k) {
    ctx.moveTo(ring[0][0] * k, ring[0][1] * k);
    for (let i = 1; i < ring.length; i++) ctx.lineTo(ring[i][0] * k, ring[i][1] * k);
    ctx.closePath();
  }
  // terrain: S.terrain = {n, step, z (mPD, rows north to south), min, max, mean, source}
  function groundAt(x, y) {
    const g = S.terrain; if (!g) return 0;
    const fx = Math.min(g.n - 1.0001, Math.max(0, x / g.step)), fy = Math.min(g.n - 1.0001, Math.max(0, y / g.step));
    const i = Math.floor(fx), j = Math.floor(fy), tx = fx - i, ty = fy - j, n = g.n;
    return g.z[j * n + i] * (1 - tx) * (1 - ty) + g.z[j * n + i + 1] * tx * (1 - ty) + g.z[(j + 1) * n + i] * (1 - tx) * ty + g.z[(j + 1) * n + i + 1] * tx * ty;
  }
  // a drawn building stands on the lowest ground under its outline
  function userBase(u) { let m = Infinity; for (const [x, y] of u.poly) m = Math.min(m, groundAt(x, y)); const c = polyCentroid(u.poly); m = Math.min(m, groundAt(c[0], c[1])); return isFinite(m) ? m : 0; }
  function rebuildMask() {
    if (solver && solver.kind === '3d') { rebuildColumns(); need3d.user = true; return; }
    const N = solver.N;
    const oc = rebuildMask.c || (rebuildMask.c = document.createElement('canvas'));
    oc.width = oc.height = N;
    const c = oc.getContext('2d', { willReadFrequently: true });
    c.clearRect(0, 0, N, N); c.fillStyle = '#000';
    const k = N / DOMAIN_M;
    c.beginPath();
    for (const b of buildings) if (blocks(b) && !S.removed.has(b.id)) for (const r of b.rings) pathRing(c, r, k);
    for (const u of S.user) pathRing(c, u.poly, k);
    c.fill('nonzero');
    need3d.user = true;
    const px = c.getImageData(0, 0, N, N).data;
    maskArr = new Uint8Array(N * N);
    for (let i = 0; i < N * N; i++) maskArr[i] = px[i * 4 + 3] >= 128 ? 1 : 0;
    solver.setMask(maskArr);
  }
  // 3D: per ground column, solid cells = terrain + the tallest building standing on it (a 2.5D height field)
  // On steep or very tall sites the air column must be taller; to stay inside the cell budget the
  // cells then get bigger (fewer across) rather than the hill or the towers being cut off.
  function dims3() {
    const t = S.terrain || { min: 0, max: 0 }, budget = BUDGET3[S.quality] || 1.15e6;
    let maxTop = t.max;
    for (const b of buildings) if (blocks(b)) maxTop = Math.max(maxTop, b.top != null ? b.top : t.min + b.h);
    const need = Math.max((maxTop - t.min) + Math.max(60, 0.5 * (maxTop - t.min)), 170);
    const dxMin = DOMAIN_M / (Q3[S.quality] || 128);
    const dxFit = Math.cbrt(DOMAIN_M * DOMAIN_M * need / budget);
    let nx = Math.floor(DOMAIN_M / Math.max(dxMin, dxFit) / 8) * 8;
    nx = Math.max(64, Math.min(Q3[S.quality] || 128, nx));
    const dx = DOMAIN_M / nx, zFloor = t.min - dx;
    const nzMax = Math.max(24, Math.floor(budget / (nx * nx) / 8) * 8);
    const nz = Math.min(nzMax, Math.max(24, Math.ceil((need + dx) / dx / 8) * 8));
    return { nx, ny: nx, nz, dx, zFloor, coarsened: nx < (Q3[S.quality] || 128), clipped: maxTop > zFloor + nz * dx };
  }
  function rebuildColumns() {
    const { nx, ny, nz, dx, zFloor } = g3, NC = nx * ny;
    const gnd = new Float32Array(NC), top = new Float32Array(NC), cover = new Uint8Array(NC);
    for (let j = 0; j < ny; j++) for (let i = 0; i < nx; i++) { const z = groundAt((i + 0.5) * dx, (j + 0.5) * dx); gnd[j * nx + i] = z; top[j * nx + i] = z; }
    const stamp = (rings, bb, zTop) => {
      const i0 = Math.max(0, Math.floor(bb[0] / dx)), i1 = Math.min(nx - 1, Math.ceil(bb[2] / dx));
      const j0 = Math.max(0, Math.floor(bb[1] / dx)), j1 = Math.min(ny - 1, Math.ceil(bb[3] / dx));
      for (let j = j0; j <= j1; j++) for (let i = i0; i <= i1; i++) {
        const x = (i + 0.5) * dx, y = (j + 0.5) * dx;
        if (rings.some((r) => inPoly(x, y, r))) { const k = j * nx + i; cover[k] = 1; if (zTop > top[k]) top[k] = zTop; }
      }
    };
    for (const b of buildings) {
      if (!blocks(b) || S.removed.has(b.id)) continue;
      stamp(b.rings, b.bb, b.top != null ? b.top : groundAt((b.bb[0] + b.bb[2]) / 2, (b.bb[1] + b.bb[3]) / 2) + b.h);
    }
    for (const u of S.user) {
      let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
      for (const [x, y] of u.poly) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      stamp([u.poly], [x0, y0, x1, y1], userBase(u) + u.h);
    }
    const h32 = new Uint32Array(NC), hG = new Uint16Array(NC), hT = new Uint16Array(NC);
    maskArr = new Uint8Array(NC);
    for (let k = 0; k < NC; k++) {
      const g = Math.min(nz - 2, Math.max(1, Math.round((gnd[k] - zFloor) / dx)));
      let t = Math.max(g, Math.round((top[k] - zFloor) / dx));
      if (cover[k] && t <= g) t = g + 1;
      t = Math.min(nz, t);
      hG[k] = g; hT[k] = t; h32[k] = t | (g << 16);
      maskArr[k] = t > g ? 1 : 0;
    }
    g3.hG = hG; g3.hT = hT;
    s3.setHeights(h32);
  }
  function inPoly(x, y, ring) {
    let ins = false;
    for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
      const xi = ring[i][0], yi = ring[i][1], xj = ring[j][0], yj = ring[j][1];
      if ((yi > y) !== (yj > y) && x < (xj - xi) * (y - yi) / (yj - yi) + xi) ins = !ins;
    }
    return ins;
  }
  function hitUser(x, y) { for (let i = S.user.length - 1; i >= 0; i--) if (inPoly(x, y, S.user[i].poly)) return S.user[i]; return null; }
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
  // 5 m terrain contours in plan (10 m drawn heavier), cached per site
  function contourSegs() {
    const g = S.terrain; if (!g) return [];
    if (g._segs) return g._segs;
    const n = g.n, segs = [];
    if (g.max - g.min >= 4) {
      for (let lev = Math.ceil(g.min / 5) * 5; lev <= g.max; lev += 5) {
        const major = lev % 10 === 0, L = [];
        for (let j = 0; j < n - 1; j++) for (let i = 0; i < n - 1; i++) {
          const c = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].map(([a2, b2]) => [a2, b2, g.z[b2 * n + a2]]);
          const pts = [];
          for (let k = 0; k < 4; k++) {
            const p1 = c[k], p2 = c[(k + 1) % 4];
            if ((p1[2] < lev) !== (p2[2] < lev)) { const t = (lev - p1[2]) / (p2[2] - p1[2]); pts.push([(p1[0] + (p2[0] - p1[0]) * t) * g.step, (p1[1] + (p2[1] - p1[1]) * t) * g.step]); }
          }
          if (pts.length >= 2) L.push(pts[0], pts[1]);
          if (pts.length === 4) L.push(pts[2], pts[3]);
        }
        segs.push({ lev, major, L });
      }
    }
    return (g._segs = segs);
  }
  const handleR = () => (window.matchMedia && matchMedia('(pointer: coarse)').matches ? 11 : 6);
  function drawBuildings() {
    need3d.user = true;
    const c = ctxB, k = S.size * S.dpr / DOMAIN_M, dark = isDark();
    c.setTransform(1, 0, 0, 1, 0, 0);
    c.clearRect(0, 0, cvBld.width, cvBld.height);
    c.lineJoin = 'round';
    for (const sg of contourSegs()) {
      c.beginPath();
      for (let i = 0; i < sg.L.length; i += 2) { c.moveTo(sg.L[i][0] * k, sg.L[i][1] * k); c.lineTo(sg.L[i + 1][0] * k, sg.L[i + 1][1] * k); }
      c.strokeStyle = dark ? `rgba(255,255,255,${sg.major ? 0.32 : 0.16})` : `rgba(40,38,34,${sg.major ? 0.34 : 0.16})`;
      c.lineWidth = (sg.major ? 0.9 : 0.6) * S.dpr; c.stroke();
    }
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
      c.beginPath(); pathRing(c, u.poly, k);
      c.fillStyle = 'rgba(217,119,87,.92)'; c.fill();
      c.strokeStyle = u === S.sel ? '#111' : '#7a3a22'; c.lineWidth = (u === S.sel ? 2.2 : 1.2) * S.dpr; c.stroke();
      const ce = polyCentroid(u.poly);
      c.fillStyle = '#fff'; c.font = `600 ${11 * S.dpr}px ui-sans-serif,system-ui,sans-serif`; c.textAlign = 'center'; c.textBaseline = 'middle';
      c.fillText(`${u.h} m`, ce[0] * k, ce[1] * k);
    }
    // handles on the selected drawn building: squares on corners, small circles at edge midpoints (drag one to add a corner)
    if (S.sel && S.mode === 'look') {
      const P = S.sel.poly, r = handleR() * S.dpr;
      for (let i = 0; i < P.length; i++) {
        const a = P[i], b = P[(i + 1) % P.length], mx = (a[0] + b[0]) / 2 * k, my = (a[1] + b[1]) / 2 * k;
        if (Math.hypot((b[0] - a[0]) * k, (b[1] - a[1]) * k) >= 26 * S.dpr) {
          c.beginPath(); c.arc(mx, my, r * 0.55, 0, Math.PI * 2); c.fillStyle = 'rgba(255,255,255,.9)'; c.fill(); c.strokeStyle = '#7a3a22'; c.lineWidth = 1 * S.dpr; c.stroke();
        }
      }
      P.forEach(([x, y], i) => {
        c.fillStyle = i === S.selV ? '#111' : '#fff'; c.strokeStyle = '#111'; c.lineWidth = 1.2 * S.dpr;
        c.fillRect(x * k - r / 2, y * k - r / 2, r, r); c.strokeRect(x * k - r / 2, y * k - r / 2, r, r);
      });
    }
    // the outline being drawn: placed points, rubber band to the cursor, close hint on the first point
    const D = S.draft;
    if (D) {
      const pts = D.rect ? rectPts(D.pts[0], D.cur) : D.pts.concat([D.cur]);
      c.beginPath(); c.moveTo(pts[0][0] * k, pts[0][1] * k);
      for (let i = 1; i < pts.length; i++) c.lineTo(pts[i][0] * k, pts[i][1] * k);
      if (D.rect) c.closePath();
      c.fillStyle = 'rgba(217,119,87,.18)'; if (pts.length >= 3) c.fill();
      c.setLineDash([5 * S.dpr, 4 * S.dpr]); c.strokeStyle = '#c4552f'; c.lineWidth = 2 * S.dpr; c.stroke(); c.setLineDash([]);
      if (!D.rect && D.pts.length >= 2) {   // faint closing edge
        c.beginPath(); c.moveTo(D.cur[0] * k, D.cur[1] * k); c.lineTo(D.pts[0][0] * k, D.pts[0][1] * k);
        c.setLineDash([2 * S.dpr, 4 * S.dpr]); c.strokeStyle = 'rgba(196,85,47,.55)'; c.lineWidth = 1.2 * S.dpr; c.stroke(); c.setLineDash([]);
      }
      const r = handleR() * S.dpr;
      D.pts.forEach(([x, y], i) => {
        c.fillStyle = i === 0 && D.closeHint ? '#c4552f' : '#fff'; c.strokeStyle = '#c4552f'; c.lineWidth = 1.4 * S.dpr;
        c.fillRect(x * k - r / 2, y * k - r / 2, r, r); c.strokeRect(x * k - r / 2, y * k - r / 2, r, r);
      });
      if (D.snapped) { c.beginPath(); c.arc(D.cur[0] * k, D.cur[1] * k, r, 0, Math.PI * 2); c.strokeStyle = '#1d6fb8'; c.lineWidth = 1.4 * S.dpr; c.stroke(); }
      if (D.guide) {   // alignment guide to the first point
        c.beginPath(); c.moveTo(D.guide[0][0] * k, D.guide[0][1] * k); c.lineTo(D.guide[1][0] * k, D.guide[1][1] * k);
        c.setLineDash([3 * S.dpr, 3 * S.dpr]); c.strokeStyle = '#1d6fb8'; c.lineWidth = 1 * S.dpr; c.stroke(); c.setLineDash([]);
      }
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
    // 3D cells are coarse: give solid cells their fluid neighbours' colour so no pale halo shows round the footprints
    if (solver.kind === '3d') {
      for (let pass = 0; pass < 2; pass++) {
        const src = new Uint8ClampedArray(d);
        for (let y = 0; y < N; y++) for (let x = 0; x < N; x++) {
          const o = (y * N + x) * 4; if (src[o + 3]) continue;
          let r = 0, g = 0, b = 0, n = 0;
          for (let dy = -1; dy <= 1; dy++) for (let dx = -1; dx <= 1; dx++) {
            const xx = x + dx, yy = y + dy; if (xx < 0 || yy < 0 || xx >= N || yy >= N) continue;
            const q = (yy * N + xx) * 4; if (!src[q + 3]) continue;
            r += src[q]; g += src[q + 1]; b += src[q + 2]; n++;
          }
          if (n) { d[o] = r / n; d[o + 1] = g / n; d[o + 2] = b / n; d[o + 3] = 255; }
        }
      }
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
    const n = solver.kind === 'gpu' ? 2600 : solver.kind === '3d' ? 2000 : 1400;
    particles = Array.from({ length: n }, () => seedParticle({}));
    parts3 = [];
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
    if (view3d && S.view === '3d') {
      if (solver.kind === '3d') move3d(steps);
      else view3d.pushParticles(particles, S.showParticles);
    }
  }

  // ---------------------------------------------------------------- 3D field helpers
  const prof3 = (zaglCells) => (zaglCells >= ZREF_M / g3.dx ? 1 : Math.pow(Math.max(zaglCells, 0.05) / (ZREF_M / g3.dx), ALPHA));
  // pedestrian slice: the first fluid cell above the terrain (centre about dx/2 = 2-3 m up), as a v2-style 2D buffer,
  // scaled so that the undisturbed incoming wind at that height reads 1 (same meaning as the v2 colours)
  function fillSlice(f) {
    const { nx, ny } = g3, NC = nx * ny;
    if (!slice.buf || slice.buf.length !== NC * 4) slice.buf = new Float32Array(NC * 4);
    const b = slice.buf, u = inflow(S.dirIdx), sc = U0 / (UREF * prof3(0.5));
    for (let k = 0; k < NC; k++) {
      const o = k * 4;
      if (maskArr && maskArr[k]) { b[o + 1] = 0; b[o + 2] = 0; b[o + 3] = -1; continue; }
      if (!f) { b[o + 1] = u[0]; b[o + 2] = u[1]; b[o + 3] = 1; continue; }
      const q = (g3.hG[k] * ny * nx + k) * 4;          // ((z*ny + y)*nx + x) with k = y*nx + x
      b[o + 1] = f[q] * sc; b[o + 2] = f[q + 1] * sc; b[o + 3] = f[q + 3] < 0 ? -1 : 1;
    }
  }
  function solid3(x, y, zl) {     // zl = level in mPD
    const { nx, dx, zFloor } = g3;
    const i = Math.floor(x / dx), j = Math.floor(y / dx);
    if (i < 0 || j < 0 || i >= nx || j >= g3.ny) return true;
    return Math.floor((zl - zFloor) / dx) < g3.hT[j * nx + i];
  }
  // trilinear velocity (lattice units) at local metres x, y and level zl (mPD); null inside solids or outside
  function f3at(x, y, zl) {
    const f = field3; if (!f) return null;
    const { nx, ny, nz, dx, zFloor } = g3;
    if (x < 0 || y < 0 || x >= DOMAIN_M || y >= DOMAIN_M || solid3(x, y, zl)) return null;
    const gx = x / dx - 0.5, gy = y / dx - 0.5, gz = Math.min(nz - 1.001, (zl - zFloor) / dx - 0.5);
    if (gz < 0) return null;
    const i0 = Math.max(0, Math.min(nx - 2, Math.floor(gx))), j0 = Math.max(0, Math.min(ny - 2, Math.floor(gy))), k0 = Math.max(0, Math.min(nz - 2, Math.floor(gz)));
    const tx = Math.min(1, Math.max(0, gx - i0)), ty = Math.min(1, Math.max(0, gy - j0)), tz = Math.min(1, Math.max(0, gz - k0));
    let ux = 0, uy = 0, uz = 0;
    for (let c = 0; c < 8; c++) {
      const di = c & 1, dj = (c >> 1) & 1, dk = c >> 2;
      const w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty) * (dk ? tz : 1 - tz);
      const o = (((k0 + dk) * ny + j0 + dj) * nx + i0 + di) * 4;
      if (f[o + 3] < 0) continue;                       // wall: no-slip, counts as zero velocity
      ux += f[o] * w; uy += f[o + 1] * w; uz += f[o + 2] * w;
    }
    return [ux, uy, uz];
  }
  const LINE_SLOW = [150, 196, 232], LINE_FAST = [8, 44, 118];
  function lineColour(sp) {
    const t = Math.min(1, Math.max(0, sp / (UREF * 1.15)));
    return [0, 1, 2].map((i) => (LINE_SLOW[i] + (LINE_FAST[i] - LINE_SLOW[i]) * t) / 255);
  }
  function seed3(p) {
    for (let tries = 0; tries < 30; tries++) {
      const x = Math.random() * DOMAIN_M, y = Math.random() * DOMAIN_M;
      const h = Math.random() < 0.6 ? 2 + Math.random() * 30 : 30 + Math.random() * 120;
      const zl = groundAt(x, y) + h;
      if (!solid3(x, y, zl)) { p.x = x; p.y = y; p.zl = zl; p.age = Math.floor(Math.random() * 90); p.fresh = true; return p; }
    }
    p.x = Math.random() * DOMAIN_M; p.y = Math.random() * DOMAIN_M; p.zl = (S.terrain ? S.terrain.max : 0) + 60; p.age = 0; p.fresh = true; return p;
  }
  function move3d(steps) {
    const n = S.quality === 'light' || isPhone() ? 900 : 1600;
    if (parts3.length !== n) parts3 = Array.from({ length: n }, () => seed3({}));
    if (!field3) return;
    const gain = steps * g3.dx * 0.6, zRef = S.terrain ? S.terrain.min : 0;
    for (const p of parts3) {
      const u = f3at(p.x, p.y, p.zl);
      p.age++;
      if (!u || p.age > 160) { seed3(p); p.age = 0; }
      else {
        p.x += u[0] * gain; p.y += u[1] * gain; p.zl += u[2] * gain;
        if (p.x < 0 || p.y < 0 || p.x >= DOMAIN_M || p.y >= DOMAIN_M) { seed3(p); p.age = 0; }
        else p.c = lineColour(Math.hypot(u[0], u[1], u[2]));
      }
      p.z = p.zl - zRef;
    }
    view3d.pushParticles(parts3, S.showParticles);
  }
  // streamlines from a rake along the upwind edge, several heights above the ground
  function traceLines() {
    if (!field3 || !view3d) return;
    const b = S.dirIdx * 22.5 * Math.PI / 180, dw = [-Math.sin(b), Math.cos(b)], pp = [-dw[1], dw[0]], H = DOMAIN_M / 2;
    const inside = (x, y) => x >= 0.5 && y >= 0.5 && x <= DOMAIN_M - 0.5 && y <= DOMAIN_M - 0.5;
    const pos = [], col = [], zRef = S.terrain ? S.terrain.min : 0, ds = g3.dx * 0.6;
    const lateral = S.quality === 'light' ? 18 : 26, heights = [3, 9, 18, 32, 55, 100];
    for (let a = 0; a < lateral; a++) {
      const sOff = -H * 0.96 + (a + 0.5) * (1.92 * H / lateral);
      let x = H + pp[0] * sOff, y = H + pp[1] * sOff;
      if (!inside(x, y)) continue;
      while (inside(x - dw[0] * 4, y - dw[1] * 4)) { x -= dw[0] * 4; y -= dw[1] * 4; }
      x += dw[0] * g3.dx * 8; y += dw[1] * g3.dx * 8;           // start just inside the sponge layer
      for (const h of heights) {
        let px = x, py = y, pz = groundAt(x, y) + h;
        for (let k = 0; k < 520; k++) {
          const u1 = f3at(px, py, pz); if (!u1) break;
          const s1 = Math.hypot(u1[0], u1[1], u1[2]); if (s1 < 0.003) break;
          const mx = px + u1[0] / s1 * ds * 0.5, my = py + u1[1] / s1 * ds * 0.5, mz = pz + u1[2] / s1 * ds * 0.5;
          const u2 = f3at(mx, my, mz) || u1; const s2 = Math.hypot(u2[0], u2[1], u2[2]) || s1;
          const nx2 = px + u2[0] / s2 * ds, ny2 = py + u2[1] / s2 * ds, nz2 = pz + u2[2] / s2 * ds;
          if (!inside(nx2, ny2)) break;
          const c = lineColour(s2);
          pos.push(px - H, pz - zRef, py - H, nx2 - H, nz2 - zRef, ny2 - H);
          col.push(c[0], c[1], c[2], 0.9, c[0], c[1], c[2], 0.9);
          px = nx2; py = ny2; pz = nz2;
        }
      }
    }
    S.lines3 = new Float32Array(pos);
    view3d.setStreamlines(S.lines3, new Float32Array(col), S.showLines);
    S.lineSegs = pos.length / 6;
  }
  // for the automated check: streamline segments that pass over a footprint (above its roof) or beside it (below the roof)
  window.HKW.overAround = (poly, roof, near) => {
    const L = S.lines3; if (!L) return null;
    const H = DOMAIN_M / 2, zRef = S.terrain ? S.terrain.min : 0, nr = near || 12;
    const distOut = (x, y) => { let d = Infinity; for (let i = 0; i < poly.length; i++) { const a = poly[i], b = poly[(i + 1) % poly.length], vx = b[0] - a[0], vy = b[1] - a[1], t = Math.max(0, Math.min(1, ((x - a[0]) * vx + (y - a[1]) * vy) / (vx * vx + vy * vy))); d = Math.min(d, Math.hypot(a[0] + vx * t - x, a[1] + vy * t - y)); } return d; };
    let over = 0, around = 0, inside = 0;
    for (let i = 0; i < L.length; i += 6) {
      const x = (L[i] + L[i + 3]) / 2 + H, z = (L[i + 1] + L[i + 4]) / 2 + zRef, y = (L[i + 2] + L[i + 5]) / 2 + H;
      const ins = inPoly(x, y, poly);
      if (ins && z > roof) over++;
      else if (ins) inside++;
      else if (z < roof - 3 && distOut(x, y) < nr) around++;
    }
    return { over, around, insideBelowRoof: inside, segments: L.length / 6 };
  };
  // vertical section along the wind through the middle (moved sideways by the slider)
  const secCanvas = document.createElement('canvas');
  function paintSection() {
    if (!view3d) return;
    if (!S.showSection || !field3 || !solver || solver.kind !== '3d') { view3d.setSection(null); return; }
    const b = S.dirIdx * 22.5 * Math.PI / 180, dw = [-Math.sin(b), Math.cos(b)], pp = [-dw[1], dw[0]], H = DOMAIN_M / 2;
    const cx = H + pp[0] * S.sectionOff, cy = H + pp[1] * S.sectionOff;
    let t0 = -1e9, t1 = 1e9;                     // where the line through (cx, cy) along the wind leaves the square
    for (const [p0, d] of [[cx, dw[0]], [cy, dw[1]]]) {
      if (Math.abs(d) < 1e-9) { if (p0 < 0 || p0 > DOMAIN_M) { view3d.setSection(null); return; } continue; }
      const a1 = (0 - p0) / d, a2 = (DOMAIN_M - p0) / d; t0 = Math.max(t0, Math.min(a1, a2)); t1 = Math.min(t1, Math.max(a1, a2));
    }
    if (t1 - t0 < 50) { view3d.setSection(null); return; }
    const W = 200, Hh = g3.nz * 2, zRef = S.terrain ? S.terrain.min : 0;
    secCanvas.width = W; secCanvas.height = Hh;
    const ctx = secCanvas.getContext('2d'), img = ctx.createImageData(W, Hh), d = img.data;
    const top = g3.zFloor + g3.nz * g3.dx;
    for (let i = 0; i < W; i++) {
      const t = t0 + (t1 - t0) * (i + 0.5) / W, x = cx + dw[0] * t, y = cy + dw[1] * t, gz = groundAt(x, y);
      for (let r = 0; r < Hh; r++) {
        const zl = top - (r + 0.5) / Hh * (top - g3.zFloor), o = (r * W + i) * 4;
        const u = f3at(x, y, zl);
        if (!u) { if (zl < top - g3.dx) { d[o] = 58; d[o + 1] = 56; d[o + 2] = 52; d[o + 3] = 235; } else d[o + 3] = 0; continue; }
        const ratio = Math.hypot(u[0], u[1], u[2]) / (UREF * prof3(Math.max(0.5, (zl - gz) / g3.dx)));
        const li = Math.min(255, (ratio * 255 / LUT_MAX) | 0) * 3;
        d[o] = LUT[li]; d[o + 1] = LUT[li + 1]; d[o + 2] = LUT[li + 2]; d[o + 3] = 225;
      }
    }
    ctx.putImageData(img, 0, 0);
    view3d.setSection({ canvas: secCanvas, a: [cx + dw[0] * t0, cy + dw[1] * t0], b: [cx + dw[0] * t1, cy + dw[1] * t1], y0: g3.zFloor - zRef, y1: top - zRef });
  }
  // numbers for the automated check
  window.HKW.stats3 = () => {
    const f = field3; if (!f || !g3 || !s3) return null;
    const { nx, ny, nz } = g3, NC = nx * ny * nz;
    let nan = 0, solidVel = 0, nSolid = 0, mass = 0, nf = 0, maxU = 0, wrongSolid = 0;
    for (let c = 0; c < NC; c++) {
      const o = c * 4, r = f[o + 3];
      const x = c % nx, y = Math.floor(c / nx) % ny, z = Math.floor(c / (nx * ny));
      const shouldSolid = z < g3.hT[y * nx + x];
      if (r < 0) { nSolid++; if (!shouldSolid) wrongSolid++; if (f[o] || f[o + 1] || f[o + 2]) solidVel++; continue; }
      if (shouldSolid) wrongSolid++;
      if (!isFinite(r) || !isFinite(f[o]) || !isFinite(f[o + 1]) || !isFinite(f[o + 2])) { nan++; continue; }
      mass += r; nf++; const sp = Math.hypot(f[o], f[o + 1], f[o + 2]); if (sp > maxU) maxU = sp;
    }
    return { kind: s3.kind, nx, ny, nz, dx: +g3.dx.toFixed(2), steps: S.steps, cells: NC, solidCells: nSolid, wrongSolid, solidWithVelocity: solidVel, nan,
      meanRho: +(mass / Math.max(1, nf)).toFixed(5), maxSpeedRatio: +(maxU / UREF).toFixed(2), lineSegs: S.lineSegs || 0, clipped: !!g3.clipped };
  };
  // velocity at a point for tests: [ux, uy, uz] relative to the reference wind
  window.HKW.vel3 = (x, y, hAboveGround) => { const u = f3at(x, y, groundAt(x, y) + hAboveGround); return u ? u.map((v) => +(v / UREF).toFixed(3)) : null; };
  window.HKW.groundAt = (x, y) => groundAt(x, y);
  window.HKW.debug3 = () => ({ readBusy, readAgeMs: Math.round(performance.now() - field3At), field: !!field3, busy: s3 ? s3.busy : null, steps: S.steps, spf: S.spf, preparing: !!S.preparing, token: flowToken });

  // ---------------------------------------------------------------- loop
  let lastFrame = 0, fpsAcc = 0, fpsN = 0, fpsT = performance.now(), stepsT = 0, stepEMA = 8;
  // How long one solver step takes on this device: time a short and a long batch twice, so the fixed
  // delay before the browser reports "done" cancels out. Run once per new solver.
  async function calibrate3d(s) {
    const sleep = (ms) => new Promise((r) => setTimeout(r, ms));
    const timed = async (n) => {
      while (s.busy) await sleep(2);
      const t0 = performance.now(); s.step(n);
      while (s.busy) await sleep(1);
      return performance.now() - t0;
    };
    await timed(2);
    const a1 = await timed(4), b1 = await timed(20), a2 = await timed(4), b2 = await timed(20);
    const per = Math.max(0.2, ((b1 + b2) - (a1 + a2)) / 32);
    s.msPerStep = per;
    return { per, steps: 50 };
  }
  let busySkips = 0, busyFrames = 0;
  function readBack3d(now) {
    if (readBusy || now - field3At <= READ_MS || S.steps <= 0) return;
    readBusy = true; const tok = flowToken, t0 = now;
    readT0 = t0;
    s3.readField().then((f) => {
      readBusy = false;
      if (!f || tok !== flowToken) return;
      field3 = f; field3At = performance.now(); onField3();
    }).catch((e) => { readBusy = false; console.warn('[wind] reading the 3D field failed:', e.message); });
  }
  let readT0 = 0;
  function step3d(now) {
    readBack3d(now);                    // copy the field back for drawing, whether or not the solver is busy
    if (s3.busy) {                      // the GPU has not finished the last batch: skip, and ease off if it keeps happening
      busySkips++; busyFrames++;
      if (busyFrames > 15 && busySkips / busyFrames > 0.3) { S.spf = Math.max(2, Math.floor(S.spf * 0.85)); busySkips = busyFrames = 0; }
      moveParticles(slice.buf, Math.max(1, stepEMA)); stepEMA *= 0.9;
      return;
    }
    busyFrames++;
    const n = s3.step(S.spf);
    if (n) S.steps += n;
    stepEMA = stepEMA * 0.85 + n * 0.15;
    moveParticles(slice.buf, Math.max(1, stepEMA));
  }
  function onField3() {
    fillSlice(field3);
    if (S.showField) { paintField(slice.buf); if (view3d) view3d.markField(); }
    if (view3d && S.view === '3d') {
      const now = performance.now();
      if (now - lines3At > 450) { traceLines(); lines3At = now; }
      if (S.showSection) paintSection();
    }
  }
  function frame(now) {
    requestAnimationFrame(frame);
    if (document.hidden) return;                       // pause when the tab is hidden
    if (now - lastFrame < 1000 / FPS_CAP - 2) return;  // frame cap
    lastFrame = now;
    if (S.running && now - S.lastTouch > IDLE_PAUSE_MS) { S.running = false; S.idlePaused = true; syncRun(); }
    if (solver && S.data && !S.loading && !S.preparing && S.running) {
      if (solver.lost) {
        if (solver.kind === '3d') {
          console.warn('[wind] the 3D solver stopped (graphics device lost); showing the 2D slice');
          s3Failed = true; s3 = null; S.engine = '2d'; S.engineNote = '(The graphics card stopped the 3D solver, so this shows the 2D slice.)';
          syncQualityLabels(); restartFlow(); return;
        }
        console.warn('[wind] GPU context lost, switching to CPU'); params.set('engine', 'cpu'); restart(); return;
      }
      if (solver.kind === '3d') { step3d(now); fpsN++; }
      else {
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
    }
    if (view3d && S.view === '3d') { sync3d(); view3d.frame(); updateNorth(); }
    if (now - fpsT > 1000) {
      S.fps = fpsN * 1000 / (now - fpsT); S.sps = (S.steps - stepsT) * 1000 / (now - fpsT);
      if (S.sps < 0) S.sps = 0;
      stepsT = S.steps; fpsN = 0; fpsT = now; if (solver && !S.preparing) updateStatus();
    }
  }
  function updateStatus() {
    // how many times the wind has crossed the square (3D: the slower air near the ground takes longer)
    const flow = solver.kind === '3d' ? S.steps * UREF * 0.7 / solver.N : S.steps * U0 / solver.N;
    const pct = Math.min(100, Math.round(flow / 1.5 * 100));
    settleEl.textContent = pct < 100 ? `Settling 形成中 ${pct}%` : 'Settled 已穩定 (eddies keep moving, as real wind does)';
    settleEl.classList.toggle('ok', pct >= 100);
    statusEl.textContent = solver.kind === '3d'
      ? `${S.fps.toFixed(0)} fps · ${Math.round(S.sps || 0)} solver steps per second`
      : `${S.fps.toFixed(0)} fps · ${S.spf} steps per frame`;
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
    prepBuildings(d);
    await withTerrain(d);
    return (cache[id] = d);
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
  // terrain: LandsD 5 m digital terrain model, cut into 1 km tiles kept in this repository (prep_terrain.py)
  let dtmIdx = null;
  const dtmTiles = new Map();
  function dtmTile(name) {
    if (dtmTiles.has(name)) return dtmTiles.get(name);
    const p = (async () => {
      if (!dtmIdx.set.has(name)) return null;                   // no tile = sea, 0 mPD
      const r = await fetch(DATA_URL + 'dtm/' + name + '.bin');
      if (!r.ok) throw new Error(`terrain tile ${name}: ${r.status}`);
      const buf = await new Response(r.body.pipeThrough(new DecompressionStream('deflate'))).arrayBuffer();
      return new Int16Array(buf);
    })();
    dtmTiles.set(name, p); p.catch(() => dtmTiles.delete(name));
    return p;
  }
  async function fetchTerrain(lat, lon) {
    if (typeof DecompressionStream === 'undefined') throw new Error('this browser cannot unpack the terrain tiles');
    if (!dtmIdx) { const d = await getJSON(DATA_URL + 'dtm/index.json'); d.set = new Set(d.tiles); dtmIdx = d; }
    const I = dtmIdx, T = I.tile, n = 121, step = DOMAIN_M / (n - 1), kx = mPerDegLon(lat);
    const cxy = new Float64Array(n * n * 2), need = new Set();
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) {
      const la = lat - (j * step - DOMAIN_M / 2) / M_PER_DEG_LAT, lo = lon + (i * step - DOMAIN_M / 2) / kx;
      const [E, N] = HK80.toGrid(la, lo);
      const ci = (E - I.x0) / I.cell - 0.5, cj = (I.north - N) / I.cell - 0.5;     // cell-centre coordinates, rows from the north
      cxy[(j * n + i) * 2] = ci; cxy[(j * n + i) * 2 + 1] = cj;
      for (const a of [Math.floor(ci), Math.floor(ci) + 1]) for (const b of [Math.floor(cj), Math.floor(cj) + 1]) need.add(`${Math.floor(a / T)}_${Math.floor(b / T)}`);
    }
    const tiles = {};
    await Promise.all([...need].map(async (nm) => { tiles[nm] = await dtmTile(nm); }));
    const cell = (ci, cj) => {
      const ti = Math.floor(ci / T), tj = Math.floor(cj / T), t = tiles[`${ti}_${tj}`];
      return t ? t[(cj - tj * T) * T + (ci - ti * T)] : 0;
    };
    const z = new Float32Array(n * n);
    let mn = Infinity, mx = -Infinity, sum = 0;
    for (let k = 0; k < n * n; k++) {
      const ci = cxy[k * 2], cj = cxy[k * 2 + 1], i0 = Math.floor(ci), j0 = Math.floor(cj), tx = ci - i0, ty = cj - j0;
      const v = cell(i0, j0) * (1 - tx) * (1 - ty) + cell(i0 + 1, j0) * tx * (1 - ty) + cell(i0, j0 + 1) * (1 - tx) * ty + cell(i0 + 1, j0 + 1) * tx * ty;
      z[k] = v; mn = Math.min(mn, v); mx = Math.max(mx, v); sum += v;
    }
    return { n, step, z, min: mn, max: mx, mean: sum / (n * n), source: 'dtm' };
  }
  // if the tiles cannot load, a flat ground at the lowest LandsD building base level (said so on the page)
  function flatTerrain(blds) {
    let m = Infinity; for (const b of blds) if (b.base != null) m = Math.min(m, b.base);
    if (!isFinite(m)) m = 0;
    return { n: 2, step: DOMAIN_M, z: new Float32Array([m, m, m, m]), min: m, max: m, mean: m, source: 'flat' };
  }
  async function withTerrain(d) {
    if (d.terrain) return d;
    try { d.terrain = await fetchTerrain(d.centre_wgs84[0], d.centre_wgs84[1]); }
    catch (e) { console.info('[wind] terrain not loaded:', e.message); d.terrain = flatTerrain(d.buildings); d.terrainErr = e.message; }
    return d;
  }
  async function loadSite(lat, lon, label, onCount) {
    const key = lat.toFixed(5) + ',' + lon.toFixed(5);
    if (cache[key]) return cache[key];
    const [feats, wind] = await Promise.all([fetchBuildings(lat, lon, onCount), fetchSiteWind(lat, lon)]);
    const d = prepBuildings({ id: key, name_en: label || 'Picked on the map', name_zh: label ? '' : '地圖揀點', live: true,
      centre_wgs84: [lat, lon], size_m: DOMAIN_M, buildings: convertFeatures(feats, lat, lon), site_wind: wind,
      fetched: new Date().toISOString().slice(0, 10) });
    await withTerrain(d);
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
    S.data = d; buildings = d.buildings; S.user = []; S.removed = new Set(); S.sel = null; S.selV = null; S.draft = null; syncSel();
    S.terrain = d.terrain || flatTerrain(buildings);
    need3d.site = true; need3d.user = true; need3d.terrain = true;
    const nb = buildings.filter(blocks).length;
    const sw = d.site_wind || {};
    const tr = S.terrain;
    const terrNote = tr.source === 'dtm'
      ? ` Ground ${tr.min.toFixed(0)} to ${tr.max.toFixed(0)} mPD (LandsD 5 m terrain model).`
      : ` Terrain could not load (${escapeHTML(d.terrainErr || 'unknown')}): ground drawn flat at ${tr.min.toFixed(1)} mPD. 地形載入唔到，當平地計。`;
    const windNote = (sw.grid_id ? `Wind table: PlanD grid cell ${escapeHTML(sw.grid_id)}.` : 'No PlanD wind table for this point (sea or outside the grid): all 16 directions are shown as equally likely. 呢點冇 PlanD 風資料。') + terrNote;
    if (d.live) {
      const [la, lo] = d.centre_wgs84;
      $('#hkwAreaNote').innerHTML = `<b>${escapeHTML(d.name_en)}${d.name_zh ? ' ' + escapeHTML(d.name_zh) : ''}</b> · ${la.toFixed(5)} N, ${lo.toFixed(5)} E · ${d.size_m} m x ${d.size_m} m · ` +
        `${nb} building blocks stand in the 2 m slice (of ${buildings.length} in LandsD data, fetched live ${d.fetched}). ${windNote}` +
        (buildings.length === 0 ? '<br>No buildings here: open ground or sea. 呢度冇建築物。' : '');
    } else {
      $('#hkwAreaNote').innerHTML = `<b>${d.name_en} ${d.name_zh}</b> · ${d.size_m} m x ${d.size_m} m · ${nb} building blocks stand in the 2 m slice (of ${buildings.length} in LandsD data).${terrNote}<br>${d.note_en}<br>${d.note_zh}`;
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
    if (need3d.terrain && S.terrain) {
      const t = S.terrain;
      view3d.setTerrain({ n: t.n, step: t.step, z: t.z, zRef: t.min, zMax: t.max, zMean: t.mean });
      need3d.terrain = false; need3d.site = true; need3d.user = true;
      if (S.dirIdx !== null) view3d.setWind(S.dirIdx);
    }
    const key = [...S.removed].join(',');
    if (need3d.site || key !== lastRemovedKey) { view3d.rebuild(stateFor3d(), need3d.site); need3d.site = false; lastRemovedKey = key; }
    if (need3d.user) { for (const u of S.user) u.base = userBase(u); view3d.rebuildUser(stateFor3d()); need3d.user = false; }
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
      sync3d();
      if (S.dirIdx !== null) view3d.setWind(S.dirIdx);
      particles.forEach((p) => { p.fresh = true; });
      parts3.forEach((p) => { p.fresh = true; });
      view3d.showField(S.showField);
      if (solver && solver.kind !== '3d' && S.steps > 0 && S.showField) { paintField(solver.read()); view3d.markField(); }
      if (solver && solver.kind === '3d') { view3d.setStreamlines(null); if (field3) { onField3(); traceLines(); lines3At = performance.now(); } }
      syncEngineUI();
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
    if (S.mode === 'draw') {
      if (hit.kind === 'user') { S.sel = hit.u; }
      else if (hit.kind === 'ground') {
        const h = +$('#hkwNewH').value || 60, x = Math.round(hit.x), y = Math.round(hit.y);
        const u = { id: nextUserId++, poly: [[x - 20, y - 15], [x + 20, y - 15], [x + 20, y + 15], [x - 20, y + 15]], h };
        S.user.push(u); S.sel = u; rebuildMask();
        placeTip(hit.ev, 'Placed a 40 x 30 m block. Draw any outline in the Plan view. 放咗一座 40 x 30 米嘅樓；喺平面可以畫任何形狀。');
      } else { placeTip(hit.ev, 'Click open ground to place a 40 x 30 m block, or draw any outline in the Plan view. 撳空地放樓，或者喺平面畫。'); return; }
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
      if (hit.kind === 'user') { S.sel = hit.u; placeTip(hit.ev, `<b>Your building 你加嘅樓</b><br>${Math.round(Math.abs(polyArea(hit.u.poly)))} m², ${hit.u.h} m tall`); }
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

  // 3D solver for the current site and quality (re-made only when the grid size changes)
  async function prepare3d() {
    if (!window.HKWSolver3D) return false;
    if (s3Creating) { try { await s3Creating; } catch (e) { /* handled by its caller */ } }
    const d = dims3(), key = `${d.nx}x${d.nz}`;
    if (!s3 || s3.lost || s3Key !== key) {
      if (s3) s3.dispose();
      s3 = null; s3Key = ''; g3 = null;
      const pref = ['webgpu', 'webgl2'].includes(params.get('engine')) ? params.get('engine') : 'auto';
      s3Creating = HKWSolver3D.create({ nx: d.nx, ny: d.ny, nz: d.nz, prefer: pref, zref: ZREF_M / d.dx, alpha: ALPHA,
        sponge: Math.max(4, Math.round(d.nx * 0.06)), tau0: 0.5 + 3 * NU3, cs2: CS3 * CS3, umax: 0.25 });
      try { s3 = await s3Creating; } finally { s3Creating = null; }
      if (!s3) { s3Failed = true; S.engineNote = '(This device has no WebGPU, and no WebGL2 with five float render targets, so it shows the 2D slice only.)'; return false; }
      s3Key = key;
    }
    g3 = d;
    return true;
  }
  function adapter3() {
    return { kind: '3d', N: g3.nx, get lost() { return !s3 || s3.lost; }, init() {}, setMask() {}, step(n) { return s3.step(n); }, read() { return slice.buf; }, dispose() {} };
  }
  function restartFlow() {
    if (window.HKW_TRACE) console.log('[trace] restartFlow engine=' + S.engine + ' failed=' + s3Failed + ' ' + String(new Error().stack).split(String.fromCharCode(10)).slice(2, 6).map((l) => l.trim()).join(' | '));
    const tok = ++flowToken;
    field3 = null; parts3 = []; S.lineSegs = 0;
    if (view3d) { view3d.setStreamlines(null); view3d.setSection(null); }
    if (S.engine === '3d' && !s3Failed && window.HKWSolver3D && S.data && S.dirIdx !== null) {
      S.preparing = true;
      settleEl.textContent = 'Preparing the 3D wind 準備緊三維風場…'; settleEl.classList.remove('ok');
      prepare3d().then((ok) => {
        if (tok !== flowToken) return;
        S.preparing = false;
        if (!ok) { S.engine = '2d'; syncQualityLabels(); restartFlow(); return; }
        solver = adapter3();
        rebuildMask();
        const u = inflow(S.dirIdx);
        s3.setInflow(u[0] * UREF / U0, u[1] * UREF / U0); s3.init();
        fillSlice(null);
        const go = () => {
          if (tok !== flowToken) return;
          S.spf = Math.max(2, Math.min(60, Math.round(BATCH_MS / s3.msPerStep)));
          busySkips = busyFrames = 0;
          afterRestart();
        };
        if (s3.msPerStep) { go(); return; }
        S.preparing = true;
        const sCal = s3;
        calibrate3d(sCal).then(() => { if (tok !== flowToken) return; S.preparing = false; sCal.init(); go(); })
          .catch((e) => { if (tok !== flowToken) return; S.preparing = false; sCal.msPerStep = 10; console.info('[wind] calibration skipped:', e.message); go(); });
      }).catch((e) => {
        if (tok !== flowToken) return;
        console.warn('[wind] 3D could not start:', e.message);
        S.preparing = false; s3Failed = true; s3 = null; S.engine = '2d';
        S.engineNote = '(3D could not start here, so this shows the 2D slice.)';
        syncQualityLabels(); restartFlow();
      });
      return;
    }
    if (!solver || solver.kind === '3d') newSolver();
    rebuildMask();
    solver.init(inflow(S.dirIdx));
    afterRestart();
  }
  function afterRestart() {
    S.steps = 0; stepsT = 0; S.field = null; resetParticles(); drawBuildings(); updateWindText();
    if (S.showField && solver.kind === '3d') paintField(slice.buf);
    engineEl.textContent = engineText();
    $('#hkwQuality').disabled = solver.kind === 'cpu';
    if (view3d) { view3d.setWind(S.dirIdx); if (S.showField) view3d.markField(); }
    syncEngineUI();
    updateStatus();
  }
  function syncEngineUI() {
    const is3 = !!solver && solver.kind === '3d';
    root.querySelectorAll('.hkw-only3d').forEach((el) => { el.hidden = !is3; });
    syncQualityLabels();
    stage.classList.toggle('engine3d', is3);
  }
  function restart() {
    S.field = null;
    if (S.engine === '3d' && !s3Failed) { restartFlow(); return; }
    newSolver(); restartFlow();
  }
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
  // Drawing works like the HK landscape tool: click corner by corner; double-click, click the first
  // corner, press Enter or right-click to finish; press and drag for a quick rectangle; Esc cancels.
  // Snapping: corners of drawn buildings, right angles to the previous edge (Shift flips it), 1 m grid.
  const SNAP = { grid: true, square: true };
  const coarse = () => !!(window.matchMedia && matchMedia('(pointer: coarse)').matches);
  function toM(ev) {
    const r = cvBld.getBoundingClientRect();
    return [Math.min(DOMAIN_M, Math.max(0, (ev.clientX - r.left) / r.width * DOMAIN_M)), Math.min(DOMAIN_M, Math.max(0, (ev.clientY - r.top) / r.height * DOMAIN_M))];
  }
  function toPx(ev) { const r = cvBld.getBoundingClientRect(); return [ev.clientX - r.left, ev.clientY - r.top]; }
  const kPx = () => S.size / DOMAIN_M;                      // CSS pixels per metre
  const rectPts = (a, b) => [[a[0], a[1]], [b[0], a[1]], [b[0], b[1]], [a[0], b[1]]];
  const r1 = (v) => Math.round(v * 100) / 100;
  let drag = null, maskTimer = 0, lastDown = null;
  function scheduleMask() {
    if (maskTimer) return;
    maskTimer = setTimeout(() => { maskTimer = 0; rebuildMask(); }, 90);
  }
  function drawMsg(t) { const el = $('#hkwDrawMsg'); if (el) el.textContent = t || ''; }
  function syncSel() {
    const box = $('#hkwSel');
    if (!S.sel) { box.hidden = true; return; }
    box.hidden = false;
    $('#hkwSelH').value = S.sel.h;
    const a = Math.abs(polyArea(S.sel.poly)), base = userBase(S.sel);
    $('#hkwSelInfo').textContent = `${S.sel.poly.length} corners · ${Math.round(a).toLocaleString()} m² · stands on the ground at ${base.toFixed(1)} mPD, roof at ${(base + S.sel.h).toFixed(1)} mPD. 企喺地面 ${base.toFixed(1)} mPD。`;
  }
  const MODE_HELP = {
    look: 'Hover a building for its name and height. Click one of your orange buildings to select it: drag it to move, drag a white corner to reshape, drag a small circle on an edge to add a corner, right-click a corner to delete it. Delete removes the selected corner or building. 揀橙色樓：拖嚟搬，拉白角改形，拉邊上圓點加角。',
    draw: 'Click corner by corner. Double-click, click the first corner or press Enter to finish; Esc cancels, Backspace takes back the last corner. Press and drag for a quick rectangle. Shift flips right-angle snapping. 逐點撳畫外形，雙擊或撳返第一點完成；按住拖可畫長方形。',
    remove: 'Click an orange building to delete it, or an existing building to take it away (what if it was not there?). Click again to bring it back.',
  };
  function setMode(m) {
    if (m !== 'draw') S.draft = null;
    S.mode = m;
    if (m !== 'look') S.selV = null;
    root.querySelectorAll('[data-mode]').forEach((b) => b.classList.toggle('on', b.dataset.mode === m));
    cvBld.style.cursor = m === 'draw' ? 'crosshair' : (m === 'remove' ? 'not-allowed' : 'default');
    $('#hkwModeHelp').textContent = MODE_HELP[m];
    drawMsg('');
    drawBuildings();
  }
  const prevDirOf = (P) => (P.length >= 2 ? Math.atan2(P[P.length - 1][1] - P[P.length - 2][1], P[P.length - 1][0] - P[P.length - 2][0]) : null);
  function snapPt(x, y, from, prevDir, ev, skip) {
    const r = (handleR() + 3) / kPx();
    for (const u of S.user) for (let i = 0; i < u.poly.length; i++) {
      if (skip && skip.u === u && skip.i === i) continue;
      const q = u.poly[i]; if (Math.hypot(q[0] - x, q[1] - y) < r) return { p: q.slice(), snapped: true };
    }
    const square = SNAP.square !== !!(ev && ev.shiftKey);
    if (square && from) {
      let a = Math.atan2(y - from[1], x - from[0]), L = Math.hypot(x - from[0], y - from[1]);
      a = prevDir == null ? Math.round(a / (Math.PI / 12)) * (Math.PI / 12) : prevDir + Math.round((a - prevDir) / (Math.PI / 2)) * (Math.PI / 2);
      if (SNAP.grid) L = Math.max(1, Math.round(L));
      const d = [Math.cos(a), Math.sin(a)];
      let p = [from[0] + d[0] * L, from[1] + d[1] * L], guide = null, snapped = false;
      const D = S.draft;
      if (D && D.pts.length >= 2) {   // stop where the edge lines up with the first corner, so the outline can close square
        const f = D.pts[0], t = (f[0] - from[0]) * d[0] + (f[1] - from[1]) * d[1];
        if (t > 1 && Math.abs(t - L) < r * 1.6) { p = [from[0] + d[0] * t, from[1] + d[1] * t]; guide = [p, f]; snapped = true; }
      }
      return { p: p.map(r1), snapped, guide };
    }
    if (SNAP.grid) return { p: [Math.round(x), Math.round(y)], snapped: false };
    return { p: [r1(x), r1(y)], snapped: false };
  }
  function draftTip(ev) {
    const D = S.draft; if (!D) { tip.hidden = true; return; }
    let html;
    if (D.rect) {
      const a = D.pts[0], b = D.cur, w = Math.abs(b[0] - a[0]), d = Math.abs(b[1] - a[1]);
      html = `${w.toFixed(0)} m x ${d.toFixed(0)} m · ${Math.round(w * d)} m²`;
    } else {
      const last = D.pts[D.pts.length - 1], L = Math.hypot(D.cur[0] - last[0], D.cur[1] - last[1]);
      const pd = prevDirOf(D.pts);
      let turn = '';
      if (pd != null && L > 0.01) { let t = (Math.atan2(D.cur[1] - last[1], D.cur[0] - last[0]) - pd) * 180 / Math.PI; t = ((t + 540) % 360) - 180; turn = ` · turn ${Math.abs(t).toFixed(0)}°`; }
      html = `Edge ${L.toFixed(1)} m${turn}` + (D.closeHint ? '<br><b>Click to finish 撳一下完成</b>' : '');
    }
    placeTip(ev, html);
  }
  function draftDown(ev, x, y, px, py) {
    const now = performance.now(), D = S.draft;
    if (!D) {
      const p = SNAP.grid ? [Math.round(x), Math.round(y)] : [r1(x), r1(y)];
      S.draft = { pts: [p], cur: p, down: [px, py], pressed: true, rect: false };
      lastDown = { t: now, px, py }; drawMsg(''); return;
    }
    D.pressed = true;
    const k = kPx(), f = D.pts[0];
    const near1 = Math.hypot(f[0] * k - px, f[1] * k - py) < (coarse() ? 16 : 10);
    const dbl = lastDown && now - lastDown.t < 380 && Math.hypot(lastDown.px - px, lastDown.py - py) < 6;
    lastDown = { t: now, px, py };
    if ((near1 || dbl) && D.pts.length >= 3) { finishDraft(); return; }
    if (dbl) return;
    const last = D.pts[D.pts.length - 1], s = snapPt(x, y, last, prevDirOf(D.pts), ev);
    if (Math.hypot(s.p[0] - last[0], s.p[1] - last[1]) < 0.5) return;
    D.pts.push(s.p); D.cur = s.p;
  }
  function draftMove(ev, x, y, px, py) {
    const D = S.draft; if (!D) return;
    if (D.pressed && D.pts.length === 1 && Math.hypot(px - D.down[0], py - D.down[1]) > 6) D.rect = true;
    const last = D.pts[D.pts.length - 1];
    const s = D.rect ? { p: SNAP.grid ? [Math.round(x), Math.round(y)] : [r1(x), r1(y)] } : snapPt(x, y, last, prevDirOf(D.pts), ev);
    D.cur = s.p; D.snapped = !!s.snapped; D.guide = s.guide || null;
    const f = D.pts[0], k = kPx();
    D.closeHint = !D.rect && D.pts.length >= 3 && Math.hypot(f[0] * k - px, f[1] * k - py) < (coarse() ? 16 : 10);
  }
  function finishDraft() {
    const D = S.draft; if (!D) return;
    S.draft = null; tip.hidden = true;
    let P = D.rect ? rectPts(D.pts[0], D.cur) : D.pts.map((p) => p.slice());
    if (P.length < 3 || Math.abs(polyArea(P)) < 20) { drawMsg('Too small (under 20 m²), not added. 太細，冇加到。'); drawBuildings(); return; }
    if (polySelfX(P)) { drawMsg('The edges cross. Draw it again. 邊交叉咗，請重畫。'); drawBuildings(); return; }
    if (polyArea(P) < 0) P = P.reverse();
    const u = { id: nextUserId++, poly: P, h: +$('#hkwNewH').value || 60 };
    S.user.push(u); S.sel = u; S.selV = null;
    setMode('look'); syncSel(); drawBuildings(); rebuildMask();
    drawMsg(`Added a ${P.length}-corner building, ${Math.round(Math.abs(polyArea(P)))} m². Set its height below. 加咗一座 ${P.length} 角嘅樓，喺下面改高度。`);
  }
  function deleteVertex(u, i) {
    if (!u || u.poly.length <= 3) { drawMsg('A building needs at least 3 corners. 最少要 3 個角。'); return; }
    u.poly.splice(i, 1); S.selV = null; syncSel(); drawBuildings(); rebuildMask();
  }
  function deleteUser(u) { S.user = S.user.filter((v) => v !== u); if (S.sel === u) { S.sel = null; S.selV = null; } syncSel(); drawBuildings(); rebuildMask(); }
  cvBld.addEventListener('pointerdown', (ev) => {
    if (ev.button === 2) return;
    const [x, y] = toM(ev), [px, py] = toPx(ev);
    const capture = () => { try { cvBld.setPointerCapture(ev.pointerId); } catch (e) { /* synthetic events */ } };
    if (S.mode === 'draw') { draftDown(ev, x, y, px, py); capture(); ev.preventDefault(); drawBuildings(); draftTip(ev); return; }
    if (S.mode === 'remove') {
      const u = hitUser(x, y);
      if (u) deleteUser(u);
      else {
        const b = hitBuilding(x, y);
        if (b) { if (S.removed.has(b.id)) S.removed.delete(b.id); else S.removed.add(b.id); }
        syncSel(); drawBuildings(); rebuildMask();
      }
      return;
    }
    const k = kPx(), R = ev.pointerType === 'touch' || coarse() ? 16 : 9;
    if (S.sel) {
      const P = S.sel.poly;
      const vi = P.findIndex((q) => Math.hypot(q[0] * k - px, q[1] * k - py) < R);
      if (vi >= 0) { S.selV = vi; drag = { kind: 'vertex', u: S.sel, i: vi, orig: P.map((p) => p.slice()) }; capture(); drawBuildings(); return; }
      for (let i = 0; i < P.length; i++) {
        const a = P[i], b = P[(i + 1) % P.length], mx = (a[0] + b[0]) / 2, my = (a[1] + b[1]) / 2;
        if (Math.hypot((b[0] - a[0]) * k, (b[1] - a[1]) * k) >= 26 && Math.hypot(mx * k - px, my * k - py) < R) {
          const orig = P.map((p) => p.slice());
          P.splice(i + 1, 0, [r1(mx), r1(my)]); S.selV = i + 1;
          drag = { kind: 'vertex', u: S.sel, i: i + 1, orig, moved: true }; capture(); drawBuildings(); return;
        }
      }
    }
    const u = hitUser(x, y);
    if (u) { S.sel = u; S.selV = null; drag = { kind: 'move', u, x0: x, y0: y, orig: u.poly.map((p) => p.slice()) }; capture(); syncSel(); drawBuildings(); return; }
    if (S.sel) { S.sel = null; S.selV = null; syncSel(); drawBuildings(); }
  });
  cvBld.addEventListener('pointermove', (ev) => {
    const [x, y] = toM(ev), [px, py] = toPx(ev);
    if (S.mode === 'draw') { if (S.draft) { draftMove(ev, x, y, px, py); drawBuildings(); draftTip(ev); } else tip.hidden = true; return; }
    if (drag && drag.kind === 'vertex') {
      const P = drag.u.poly, n = P.length, s = snapPt(x, y, null, null, ev, { u: drag.u, i: drag.i });
      P[drag.i] = s.p; drag.moved = true;
      const prev = P[(drag.i + n - 1) % n];
      placeTip(ev, `Corner ${drag.i + 1} · edge ${Math.hypot(s.p[0] - prev[0], s.p[1] - prev[1]).toFixed(1)} m`);
      drawBuildings(); scheduleMask(); return;
    }
    if (drag && drag.kind === 'move') {
      let dx = x - drag.x0, dy = y - drag.y0; if (SNAP.grid) { dx = Math.round(dx); dy = Math.round(dy); }
      if (Math.abs(dx) + Math.abs(dy) > 0) drag.moved = true;
      drag.u.poly = drag.orig.map(([a, b]) => [r1(a + dx), r1(b + dy)]);
      drawBuildings(); scheduleMask(); return;
    }
    showTip(ev, x, y);
  });
  const endDrag = () => {
    if (S.draft) { S.draft.pressed = false; if (S.draft.rect) finishDraft(); drawBuildings(); return; }
    if (!drag) return;
    if (drag.kind === 'vertex' && drag.moved && polySelfX(drag.u.poly)) { drag.u.poly = drag.orig; drawMsg('That would make the edges cross, so the corner went back. 邊會交叉，角已還原。'); }
    if (drag.moved) { rebuildMask(); syncSel(); }
    tip.hidden = true; drag = null; drawBuildings();
  };
  cvBld.addEventListener('pointerup', endDrag);
  cvBld.addEventListener('pointercancel', endDrag);
  cvBld.addEventListener('pointerleave', () => { if (!drag) tip.hidden = true; });
  cvBld.addEventListener('dblclick', (ev) => ev.preventDefault());
  cvBld.addEventListener('contextmenu', (ev) => {     // right-click: finish the outline, or delete the corner under the cursor
    ev.preventDefault();
    if (S.draft) { if (S.draft.pts.length >= 3) finishDraft(); drawBuildings(); return; }
    if (S.sel) {
      const [px, py] = toPx(ev), k = kPx();
      const i = S.sel.poly.findIndex((q) => Math.hypot(q[0] * k - px, q[1] * k - py) < 10);
      if (i >= 0) deleteVertex(S.sel, i);
    }
  });
  document.addEventListener('keydown', (ev) => {
    const tg = ev.target && ev.target.tagName;
    if (tg === 'INPUT' || tg === 'TEXTAREA' || tg === 'SELECT') return;
    if (S.draft) {
      if (ev.key === 'Enter') { ev.preventDefault(); if (S.draft.pts.length >= 3) finishDraft(); else drawMsg('Place at least 3 corners. 最少 3 個角。'); }
      else if (ev.key === 'Escape') { S.draft = null; tip.hidden = true; drawMsg('Cancelled. 取消咗。'); }
      else if (ev.key === 'Backspace') { ev.preventDefault(); S.draft.pts.pop(); if (!S.draft.pts.length) S.draft = null; }
      else return;
      drawBuildings(); return;
    }
    if (S.sel && (ev.key === 'Delete' || ev.key === 'Backspace')) {
      ev.preventDefault();
      if (S.selV != null && S.sel.poly[S.selV]) deleteVertex(S.sel, S.selV); else deleteUser(S.sel);
    } else if (S.sel && ev.key === 'Escape') { S.sel = null; S.selV = null; syncSel(); drawBuildings(); }
  });

  function showTip(ev, x, y) {
    if (S.mode !== 'look' && S.mode !== 'remove') { tip.hidden = true; return; }
    let html = '';
    const u = hitUser(x, y);
    const b = u ? null : hitBuilding(x, y);
    if (u) html = `<b>Your building 你加嘅樓</b><br>${Math.round(Math.abs(polyArea(u.poly)))} m², ${u.h} m tall`;
    else if (b) html = buildingTip(b);
    else if (solver && S.steps > 0 && solver.read()) {
      const v = sampleU(solver.read(), x, y);
      if (v) {
        const r = Math.hypot(v[0], v[1]) / U0;
        const w = r < 0.35 ? 'Calm 好靜' : r < 0.75 ? 'Calmer than the incoming wind 比來風弱' : r < 1.1 ? 'About the same as the incoming wind 同來風差唔多' : 'Faster than the incoming wind 比來風強';
        const settling = solver.kind === '3d' ? S.steps * UREF * 0.7 / solver.N < 1.5 : S.steps * U0 / solver.N < 1.5;
        html = `<b>${w}</b>${settling ? '<br><span class="small">still settling</span>' : ''}<br><span class="small">ground ${groundAt(x, y).toFixed(0)} mPD</span>`;
      }
    }
    if (!html) { tip.hidden = true; return; }
    placeTip(ev, html);
  }

  // ---------------------------------------------------------------- controls
  root.addEventListener('click', (ev) => {
    const t = ev.target.closest('[data-area],[data-dir],[data-season],[data-mode],[data-wengine]');
    if (!t || !root.contains(t)) return;
    if (t.dataset.area) selectArea(t.dataset.area);
    else if (t.dataset.dir !== undefined) setDir(+t.dataset.dir);
    else if (t.dataset.season) { S.season = t.dataset.season; drawRose(); updateWindText(); }
    else if (t.dataset.mode) setMode(t.dataset.mode);
    else if (t.dataset.wengine === '3d' || t.dataset.wengine === '2d') { touch(); S.engine = t.dataset.wengine; if (S.engine === '3d') { s3Failed = false; S.engineNote = ''; } syncQualityLabels(); restartFlow(); }
  });
  runBtn.addEventListener('click', () => { S.running = !S.running; S.idlePaused = false; touch(); syncRun(); });
  $('#hkwReset').addEventListener('click', () => { touch(); restartFlow(); });
  $('#hkwUndo').addEventListener('click', () => { S.user = []; S.removed = new Set(); S.sel = null; S.selV = null; S.draft = null; syncSel(); drawBuildings(); rebuildMask(); });
  $('#hkwShowField').addEventListener('change', (e) => { S.showField = e.target.checked; cvField.style.visibility = S.showField ? '' : 'hidden'; if (view3d) view3d.showField(S.showField); if (S.showField && solver && solver.kind === '3d' && slice.buf) paintField(slice.buf); });
  $('#hkwShowParticles').addEventListener('change', (e) => { S.showParticles = e.target.checked; });
  $('#hkwShowLines').addEventListener('change', (e) => { S.showLines = e.target.checked; if (view3d) { if (S.showLines) traceLines(); else view3d.setStreamlines(null); } });
  $('#hkwShowSection').addEventListener('change', (e) => { S.showSection = e.target.checked; $('#hkwSecOff').disabled = !S.showSection; paintSection(); });
  $('#hkwSecOff').addEventListener('input', (e) => { S.sectionOff = +e.target.value; $('#hkwSecOffV').textContent = `${S.sectionOff > 0 ? '+' : ''}${S.sectionOff} m`; paintSection(); });
  $('#hkwSnapGrid').addEventListener('change', (e) => { SNAP.grid = e.target.checked; });
  $('#hkwSnapSquare').addEventListener('change', (e) => { SNAP.square = e.target.checked; });
  $('#hkwQuality').value = S.quality;
  $('#hkwQuality').addEventListener('change', (e) => { S.quality = e.target.value; restart(); });
  $('#hkwSelH').addEventListener('input', (e) => {
    if (!S.sel) return; const v = +e.target.value; if (!isFinite(v) || v < 1) return;
    S.sel.h = Math.min(500, v); syncSel(); drawBuildings(); scheduleMask();
  });
  $('#hkwSelDel').addEventListener('click', () => { if (S.sel) deleteUser(S.sel); });
  window.addEventListener('resize', () => { resize(); resetParticles(); });
  new MutationObserver(() => { drawBuildings(); if (view3d) view3d.setTheme(isDark()); }).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // test hooks (used by the automated check): add a rectangle (optionally rotated) or any outline
  window.HKW.addBlock = (cx, cy, w, d, h, rot) => {
    const a = (rot || 0) * Math.PI / 180, c = Math.cos(a), s = Math.sin(a);
    const poly = [[-1, -1], [1, -1], [1, 1], [-1, 1]].map(([px, py]) => { const x = px * w / 2, y = py * d / 2; return [r1(cx + x * c - y * s), r1(cy + x * s + y * c)]; });
    const u = { id: nextUserId++, poly, h: h || 60 };
    S.user.push(u); S.sel = u; syncSel(); drawBuildings(); rebuildMask(); return u.id;
  };
  window.HKW.addPoly = (poly, h) => { const u = { id: nextUserId++, poly: poly.map((p) => p.slice()), h: h || 60 }; S.user.push(u); S.sel = u; syncSel(); drawBuildings(); rebuildMask(); return u.id; };
  window.HKW.clearBlocks = () => { S.user = []; S.sel = null; S.selV = null; syncSel(); drawBuildings(); rebuildMask(); };
  window.HKW.selectSite = selectSite; window.HKW.selectArea = selectArea; window.HKW.setView = setView; window.HKW.setDir = setDir;
  window.HKW.setMode = setMode; window.HKW.paintSection = paintSection; window.HKW.traceLines = traceLines;

  // ---------------------------------------------------------------- start
  setMode('look'); syncRun(); resize(); newSolver(); syncQualityLabels();
  initMap();
  $('#hkwViewHelp').textContent = VIEW_HELP['2d'];
  const plat = parseFloat(params.get('lat')), plon = parseFloat(params.get('lon'));
  if (isFinite(plat) && isFinite(plon)) { if (params.get('q')) qEl.value = params.get('q'); selectSite(plat, plon, params.get('q')); }
  else selectArea(S.areaId);
  if (params.get('view') === '3d') setView('3d');
  requestAnimationFrame(frame);
})();
