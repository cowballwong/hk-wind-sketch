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
    data: null, areaId: params.get('area') || 'kwun_tong', dirIdx: null, season: 'annual',
    quality: params.get('quality') || 'normal', running: true, idlePaused: false,
    showField: true, showParticles: true, mode: 'look', user: [], removed: new Set(), sel: null,
    steps: 0, lastTouch: performance.now(), spf: 20, size: 600, dpr: 1, field: null, fps: 0,
  };
  let solver = null, maskArr = null, particles = [], buildings = [], nextUserId = 1;
  window.HKW = { state: S, get solver() { return solver; } };   // for tests and the curious

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
      if (!maskArr || !maskArr[cy * N + cx]) { p.x = x; p.y = y; p.age = Math.floor(Math.random() * 80); return p; }
    }
    p.x = Math.random() * DOMAIN_M; p.y = Math.random() * DOMAIN_M; p.age = 0; return p;
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
    if (!S.showParticles) return;
    c.strokeStyle = S.showField ? 'rgba(255,255,255,0.78)' : (isDark() ? 'rgba(235,232,225,.8)' : 'rgba(30,30,30,.7)');
    c.lineWidth = 1.1 * S.dpr; c.beginPath();
    const gain = steps * dx * 0.55;
    for (const p of particles) {
      const u = sampleU(buf, p.x, p.y);
      p.age++;
      if (!u || p.age > 140) { seedParticle(p); p.age = 0; continue; }
      const nx = p.x + u[0] * gain, ny = p.y + u[1] * gain;
      if (nx < 0 || ny < 0 || nx >= DOMAIN_M || ny >= DOMAIN_M) { seedParticle(p); p.age = 0; continue; }
      c.moveTo(p.x * k, p.y * k); c.lineTo(nx * k, ny * k);
      p.x = nx; p.y = ny;
    }
    c.stroke();
  }

  // ---------------------------------------------------------------- loop
  let lastFrame = 0, fpsAcc = 0, fpsN = 0, fpsT = performance.now();
  function frame(now) {
    requestAnimationFrame(frame);
    if (!solver || !S.data) return;
    if (document.hidden) return;
    if (now - lastFrame < 1000 / FPS_CAP - 2) return;
    lastFrame = now;
    if (S.running && now - S.lastTouch > IDLE_PAUSE_MS) { S.running = false; S.idlePaused = true; syncRun(); }
    if (!S.running) return;
    if (solver.lost) { console.warn('[wind] GPU context lost, switching to CPU'); params.set('engine', 'cpu'); restart(); return; }
    const t0 = performance.now();
    solver.step(S.spf);
    S.steps += S.spf;
    const buf = solver.read();
    if (S.showField) paintField(buf);
    moveParticles(buf, S.spf);
    const work = performance.now() - t0;
    // keep the page light: aim for well under one frame of work
    if (work > 24 && S.spf > 2) S.spf = Math.max(2, Math.floor(S.spf * 0.8));
    else if (work < 10 && S.spf < (solver.kind === 'gpu' ? 36 : 8)) S.spf += 1;
    fpsN++;
    if (now - fpsT > 1000) { S.fps = fpsN * 1000 / (now - fpsT); fpsN = 0; fpsT = now; updateStatus(); }
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
  async function loadArea(id) {
    if (cache[id]) return cache[id];
    let d;
    if (window.HKW_DATA && window.HKW_DATA[id]) d = window.HKW_DATA[id];
    else {
      const base = window.HKW_DATA_URL || 'data/';
      const r = await fetch(base + id + '.json');
      if (!r.ok) throw new Error('could not load ' + id);
      d = await r.json();
    }
    for (const b of d.buildings) {
      let x0 = 1e9, y0 = 1e9, x1 = -1e9, y1 = -1e9;
      for (const r of b.rings) for (const [x, y] of r) { x0 = Math.min(x0, x); y0 = Math.min(y0, y); x1 = Math.max(x1, x); y1 = Math.max(y1, y); }
      b.bb = [x0, y0, x1, y1];
    }
    cache[id] = d; return d;
  }
  function seasonFreq() {
    const s = S.data && S.data.site_wind && S.data.site_wind.seasons[S.season];
    if (!s) return null;
    const m = {}; s.dirs.forEach((d, i) => { m[d] = s.freq[i]; });
    return DIRS.map(([d]) => m[d] || 0);
  }
  function topDir(season) {
    const s = S.data.site_wind.seasons[season]; if (!s) return 4;
    let bi = 0; s.freq.forEach((f, i) => { if (f > s.freq[bi]) bi = i; });
    return Math.max(0, DIRS.findIndex(([d]) => d === s.dirs[bi]));
  }
  async function selectArea(id) {
    S.areaId = id;
    root.querySelectorAll('[data-area]').forEach((b) => b.classList.toggle('on', b.dataset.area === id));
    $('#hkwAreaNote').textContent = 'Loading 載入中…';
    let d;
    try { d = await loadArea(id); } catch (e) { $('#hkwAreaNote').textContent = 'Could not load this area: ' + e.message; return; }
    S.data = d; buildings = d.buildings; S.user = []; S.removed = new Set(); S.sel = null; syncSel();
    const nb = buildings.filter(blocks).length;
    $('#hkwAreaNote').innerHTML = `<b>${d.name_en} ${d.name_zh}</b> · ${d.size_m} m x ${d.size_m} m · ${nb} building blocks stand in the 2 m slice (of ${buildings.length} in LandsD data).<br>${d.note_en}<br>${d.note_zh}`;
    if (S.dirIdx === null) {
      const want = params.get('dir');
      const wi = want ? DIRS.findIndex(([e]) => e === want.toUpperCase()) : -1;
      S.dirIdx = wi >= 0 ? wi : topDir('annual');
    }
    drawRose(); restartFlow();
  }
  function restartFlow() {
    if (!solver) newSolver();
    rebuildMask();
    solver.init(inflow(S.dirIdx));
    S.steps = 0; resetParticles(); drawBuildings(); updateWindText();
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
  $('#hkwShowField').addEventListener('change', (e) => { S.showField = e.target.checked; cvField.style.visibility = S.showField ? 'visible' : 'hidden'; });
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
  new MutationObserver(() => drawBuildings()).observe(document.documentElement, { attributes: true, attributeFilter: ['data-theme'] });

  // test hook: add a block programmatically (used by the automated check)
  window.HKW.addBlock = (cx, cy, w, d, h, rot) => {
    const u = { id: nextUserId++, cx, cy, w, d, h: h || 60, rot: rot || 0 };
    S.user.push(u); S.sel = u; syncSel(); drawBuildings(); rebuildMask(); return u.id;
  };
  window.HKW.clearBlocks = () => { S.user = []; S.sel = null; syncSel(); drawBuildings(); rebuildMask(); };

  // ---------------------------------------------------------------- start
  setMode('look'); syncRun(); resize(); newSolver();
  selectArea(S.areaId);
  requestAnimationFrame(frame);
})();
