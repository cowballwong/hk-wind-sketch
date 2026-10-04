/* HK Wind Sketch v3: 3D lattice Boltzmann wind solver (D3Q19) for the browser.
 *
 * Scheme: D3Q19, regularised BGK collision with a Smagorinsky eddy viscosity,
 * pull streaming, half-way bounce-back on solids, free-slip top, and an inflow
 * equilibrium with a power-law profile on the four side boundaries plus a thin
 * sponge layer (so any of the 16 wind directions works on a fixed grid).
 *
 * The solid geometry is a 2.5D height field: per ground column, the number of
 * solid cells from the bottom (terrain plus any building standing there) and the
 * number that are terrain only (used for "height above ground" in the inflow
 * profile). Overhangs and open ground floors are not represented.
 *
 * Two back ends with one interface:
 *   WebGPU  compute shaders, populations in storage buffers (structure of arrays)
 *   WebGL2  fragment shaders, z-slices tiled into a 2D atlas, 5 float targets
 * Shader sources are generated with the 19 directions unrolled.
 *
 * Own code, MIT licence. No code from FluidX3D or any other LBM package.
 */
(function () {
  'use strict';
  const C = [[0, 0, 0], [1, 0, 0], [-1, 0, 0], [0, 1, 0], [0, -1, 0], [0, 0, 1], [0, 0, -1],
    [1, 1, 0], [-1, -1, 0], [1, -1, 0], [-1, 1, 0], [1, 0, 1], [-1, 0, -1], [1, 0, -1], [-1, 0, 1],
    [0, 1, 1], [0, -1, -1], [0, 1, -1], [0, -1, 1]];
  const Q = 19;
  const W = C.map((c, i) => (i === 0 ? 1 / 3 : i < 7 ? 1 / 18 : 1 / 36));
  const find = (v) => C.findIndex((d) => d[0] === v[0] && d[1] === v[1] && d[2] === v[2]);
  const OPP = C.map((c) => find([-c[0], -c[1], -c[2]]));
  const MZ = C.map((c) => find([c[0], c[1], -c[2]]));
  const COLS = [];                                   // distinct (cx, cy) neighbour columns
  C.forEach((c) => { if (!COLS.some((k) => k[0] === c[0] && k[1] === c[1])) COLS.push([c[0], c[1]]); });
  const colOf = (i) => COLS.findIndex((k) => k[0] === C[i][0] && k[1] === C[i][1]);
  const fl = (v) => { const s = String(+v.toPrecision(9)); return /[.e]/.test(s) ? s : s + '.0'; };

  // ---- shared shader text (both languages accept these expressions) --------
  // feq_i(rho, u) with u = (ux, uy, uz); returns an expression string
  const feqExpr = (i, rho, ux, uy, uz, usq) => {
    const [cx, cy, cz] = C[i];
    const terms = [];
    if (cx) terms.push((cx > 0 ? '' : '-') + ux);
    if (cy) terms.push((cy > 0 ? '' : '-') + uy);
    if (cz) terms.push((cz > 0 ? '' : '-') + uz);
    const eu = terms.length ? '(' + terms.join('+').replace(/\+-/g, '-') + ')' : '0.0';
    return `${fl(W[i])}*${rho}*(1.0+3.0*${eu}+4.5*${eu}*${eu}-1.5*${usq})`;
  };
  // collision body: given f0..f18 (vars), uin (vec3 inflow here), sponge s; defines g0..g18
  function collideBody(lang) {
    const v3 = lang === 'wgsl' ? 'vec3<f32>' : 'vec3';
    const decl = lang === 'wgsl' ? 'var' : 'float';
    const declc = lang === 'wgsl' ? 'let' : 'float';
    let s = '';
    s += `  ${decl} rho = ${C.map((c, i) => 'f' + i).join('+')};\n`;
    const mom = (a) => C.map((c, i) => (c[a] === 0 ? '' : (c[a] > 0 ? '+' : '-') + 'f' + i)).join('').replace(/^\+/, '');
    s += `  ${decl} mx = ${mom(0)};\n  ${decl} my = ${mom(1)};\n  ${decl} mz = ${mom(2)};\n`;
    s += `  ${decl} bad = 0.0;\n`;
    s += `  if (!(rho > 0.3 && rho < 3.0)) { rho = 1.0; mx = uin.x; my = uin.y; mz = 0.0; bad = 1.0; }\n`;
    s += `  ${lang === 'wgsl' ? 'var' : 'vec3'} u = ${v3}(mx, my, mz) / rho;\n`;
    s += `  ${declc} sp = length(u);\n  if (sp > UMAX) { u = u * (UMAX / sp); }\n`;
    s += `  ${declc} usq = dot(u, u);\n`;
    // non-equilibrium stress
    const comps = [['xx', 0, 0], ['yy', 1, 1], ['zz', 2, 2], ['xy', 0, 1], ['xz', 0, 2], ['yz', 1, 2]];
    for (let i = 0; i < Q; i++) s += `  ${declc} n${i} = f${i} - ${feqExpr(i, 'rho', 'u.x', 'u.y', 'u.z', 'usq')};\n`;
    for (const [nm, a, b] of comps) {
      const t = C.map((c, i) => { const p = c[a] * c[b]; return p === 0 ? '' : (p > 0 ? '+' : '-') + 'n' + i; }).join('').replace(/^\+/, '');
      s += `  ${declc} P${nm} = ${t || '0.0'};\n`;
    }
    s += `  ${declc} Qn = sqrt(Pxx*Pxx + Pyy*Pyy + Pzz*Pzz + 2.0*(Pxy*Pxy + Pxz*Pxz + Pyz*Pyz));\n`;
    s += `  ${declc} tau = 0.5*(TAU0 + sqrt(TAU0*TAU0 + 18.0*1.41421356*CS2*Qn/rho));\n`;
    s += `  ${declc} keep = (1.0 - 1.0/tau) * (1.0 - bad);\n`;
    for (let i = 0; i < Q; i++) {
      const [cx, cy, cz] = C[i];
      const q = (c1, c2) => c1 * c2;
      const t = [];
      const add = (coef, name) => { if (coef) t.push(`${fl(coef)}*P${name}`); };
      add(q(cx, cx) - 1 / 3, 'xx'); add(q(cy, cy) - 1 / 3, 'yy'); add(q(cz, cz) - 1 / 3, 'zz');
      add(2 * cx * cy, 'xy'); add(2 * cx * cz, 'xz'); add(2 * cy * cz, 'yz');
      const f1 = `${fl(W[i] * 4.5)}*(${t.join('+').replace(/\+-/g, '-')})`;
      s += `  ${declc} g${i} = mix(${feqExpr(i, 'rho', 'u.x', 'u.y', 'u.z', 'usq')} + keep*${f1}, ${feqExpr(i, '1.0', 'uin.x', 'uin.y', 'uin.z', 'uisq')}, spg);\n`;
    }
    return s;
  }

  // ======================================================================= WebGPU
  function wgslStep() {
    let s = `
struct Prm { nx: i32, ny: i32, nz: i32, init: i32, ux: f32, uy: f32, zref: f32, alpha: f32, sponge: f32, p0: f32, p1: f32, p2: f32 };
@group(0) @binding(0) var<uniform> P: Prm;
@group(0) @binding(1) var<storage, read> A: array<f32>;
@group(0) @binding(2) var<storage, read_write> B: array<f32>;
@group(0) @binding(3) var<storage, read> H: array<u32>;
override TAU0: f32 = 0.515;
override CS2: f32 = 0.0289;
override UMAX: f32 = 0.25;
@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let x = i32(gid.x); let y = i32(gid.y); let z = i32(gid.z);
  if (x >= P.nx || y >= P.ny || z >= P.nz) { return; }
  let NC = u32(P.nx * P.ny * P.nz);
  let ci = (z * P.ny + y) * P.nx + x;
  let c = u32(ci);
  let h = H[u32(y * P.nx + x)];
  let top = i32(h & 0xffffu); let gnd = i32(h >> 16u);
  if (z < top) {
${C.map((cc, i) => `    B[${i}u * NC + c] = ${fl(W[i])};`).join('\n')}
    return;
  }
  let zagl = f32(z - gnd) + 0.5;
  let prof = select(pow(max(zagl, 0.05) / P.zref, P.alpha), 1.0, zagl >= P.zref);
  let uin = vec3<f32>(P.ux * prof, P.uy * prof, 0.0);
  let uisq = dot(uin, uin);
  if (P.init == 1) {
${C.map((cc, i) => `    B[${i}u * NC + c] = ${feqExpr(i, '1.0', 'uin.x', 'uin.y', 'uin.z', 'uisq')};`).join('\n')}
    return;
  }
`;
    COLS.forEach(([cx, cy], k) => {
      s += `  let qx${k} = x - (${cx}); let qy${k} = y - (${cy});\n`;
      s += `  let in${k} = qx${k} >= 0 && qx${k} < P.nx && qy${k} >= 0 && qy${k} < P.ny;\n`;
      s += `  var t${k} = 0; if (in${k}) { t${k} = i32(H[u32(qy${k} * P.nx + qx${k})] & 0xffffu); }\n`;
    });
    for (let i = 0; i < Q; i++) {
      const [cx, cy, cz] = C[i], k = colOf(i);
      const off = `(${-cz} * P.ny * P.nx + ${-cy} * P.nx + ${-cx})`;
      s += `  var f${i}: f32;\n`;
      s += `  if (!in${k}) { f${i} = ${feqExpr(i, '1.0', 'uin.x', 'uin.y', 'uin.z', 'uisq')}; }\n`;
      if (cz < 0) {
        // source above the lid: free slip, take the mirrored population from the same layer
        s += `  else if (z + 1 >= P.nz) { if (z < t${k}) { f${i} = A[${OPP[i]}u * NC + c]; } else { f${i} = A[${MZ[i]}u * NC + u32(ci + (${-cy} * P.nx + ${-cx}))]; } }\n`;
      }
      if (cz > 0) s += `  else if (z - 1 < 0 || z - 1 < t${k}) { f${i} = A[${OPP[i]}u * NC + c]; }\n`;
      else if (cz < 0) s += `  else if (z + 1 < t${k}) { f${i} = A[${OPP[i]}u * NC + c]; }\n`;
      else if (i !== 0) s += `  else if (z < t${k}) { f${i} = A[${OPP[i]}u * NC + c]; }\n`;
      s += `  else { f${i} = A[${i}u * NC + u32(ci + ${off})]; }\n`;
    }
    s += `  let dd = f32(min(min(x, y), min(P.nx - 1 - x, P.ny - 1 - y)));
  var spg = 0.0; if (dd < P.sponge) { let tt = 1.0 - dd / P.sponge; spg = 0.25 * tt * tt; }
`;
    s += collideBody('wgsl');
    for (let i = 0; i < Q; i++) s += `  B[${i}u * NC + c] = g${i};\n`;
    s += '}\n';
    return s;
  }
  function wgslMoments() {
    return `
struct Prm { nx: i32, ny: i32, nz: i32, init: i32, ux: f32, uy: f32, zref: f32, alpha: f32, sponge: f32, p0: f32, p1: f32, p2: f32 };
@group(0) @binding(0) var<uniform> P: Prm;
@group(0) @binding(1) var<storage, read> A: array<f32>;
@group(0) @binding(2) var<storage, read_write> V: array<vec4<f32>>;
@group(0) @binding(3) var<storage, read> H: array<u32>;
@compute @workgroup_size(64, 1, 1)
fn main(@builtin(global_invocation_id) gid: vec3<u32>) {
  let x = i32(gid.x); let y = i32(gid.y); let z = i32(gid.z);
  if (x >= P.nx || y >= P.ny || z >= P.nz) { return; }
  let NC = u32(P.nx * P.ny * P.nz);
  let c = u32((z * P.ny + y) * P.nx + x);
  let top = i32(H[u32(y * P.nx + x)] & 0xffffu);
  if (z < top) { V[c] = vec4<f32>(0.0, 0.0, 0.0, -1.0); return; }
${C.map((cc, i) => `  let f${i} = A[${i}u * NC + c];`).join('\n')}
  let rho = ${C.map((c, i) => 'f' + i).join('+')};
  let mx = ${C.map((c, i) => (c[0] === 0 ? '' : (c[0] > 0 ? '+' : '-') + 'f' + i)).join('').replace(/^\+/, '')};
  let my = ${C.map((c, i) => (c[1] === 0 ? '' : (c[1] > 0 ? '+' : '-') + 'f' + i)).join('').replace(/^\+/, '')};
  let mz = ${C.map((c, i) => (c[2] === 0 ? '' : (c[2] > 0 ? '+' : '-') + 'f' + i)).join('').replace(/^\+/, '')};
  V[c] = vec4<f32>(mx / rho, my / rho, mz / rho, rho);
}
`;
  }

  async function makeWebGPU(o) {
    if (!navigator.gpu) return null;
    const adapter = await navigator.gpu.requestAdapter(gpuOpts());
    if (!adapter) return null;
    const { nx, ny, nz } = o, NC = nx * ny * nz;
    const fBytes = NC * Q * 4, vBytes = NC * 16;
    if (fBytes > adapter.limits.maxStorageBufferBindingSize || fBytes > adapter.limits.maxBufferSize) return null;
    const device = await adapter.requestDevice({ requiredLimits: {
      maxStorageBufferBindingSize: Math.max(fBytes, vBytes), maxBufferSize: Math.max(fBytes, vBytes) } });
    let lost = false;
    let disposed = false;
    device.lost.then((info) => { lost = true; if (!disposed) console.warn('[wind3d] WebGPU device lost:', info.message); });
    const consts = { TAU0: o.tau0, CS2: o.cs2, UMAX: o.umax };
    const mkPipe = async (code, constants) => {
      const module = device.createShaderModule({ code });
      const info = await module.getCompilationInfo();
      const errs = info.messages.filter((m) => m.type === 'error');
      if (errs.length) throw new Error('WGSL: ' + errs[0].message + ' (line ' + errs[0].lineNum + ')');
      return device.createComputePipelineAsync({ layout: 'auto', compute: { module, entryPoint: 'main', constants } });
    };
    let stepPipe, momPipe;
    device.pushErrorScope('validation');
    try {
      [stepPipe, momPipe] = await Promise.all([mkPipe(wgslStep(), consts), mkPipe(wgslMoments(), {})]);
    } catch (e) { device.destroy(); throw e; }
    const U = GPUBufferUsage;
    const prm = device.createBuffer({ size: 48, usage: U.UNIFORM | U.COPY_DST });
    const fA = device.createBuffer({ size: fBytes, usage: U.STORAGE });
    const fB = device.createBuffer({ size: fBytes, usage: U.STORAGE });
    const hBuf = device.createBuffer({ size: Math.max(16, nx * ny * 4), usage: U.STORAGE | U.COPY_DST });
    const vBuf = device.createBuffer({ size: vBytes, usage: U.STORAGE | U.COPY_SRC });
    const stage = device.createBuffer({ size: vBytes, usage: U.MAP_READ | U.COPY_DST });
    const err = await device.popErrorScope();
    if (err) { device.destroy(); throw new Error(err.message); }
    const bg = (pipe, a, b) => device.createBindGroup({ layout: pipe.getBindGroupLayout(0), entries: [
      { binding: 0, resource: { buffer: prm } }, { binding: 1, resource: { buffer: a } },
      { binding: 2, resource: { buffer: b } }, { binding: 3, resource: { buffer: hBuf } }] });
    const stepBG = [bg(stepPipe, fA, fB), bg(stepPipe, fB, fA)];
    const momBG = [bg(momPipe, fA, vBuf), bg(momPipe, fB, vBuf)];
    let cur = 0;                         // index of the buffer holding the latest populations (0 = fA)
    const P = new ArrayBuffer(48), Pi = new Int32Array(P), Pf = new Float32Array(P);
    const wg = [Math.ceil(nx / 64), ny, nz];
    let inflight = false, mapping = false, lastBatch = { n: 0, ms: 0 };
    function writeP(init) {
      Pi[0] = nx; Pi[1] = ny; Pi[2] = nz; Pi[3] = init ? 1 : 0;
      Pf[4] = o.ux; Pf[5] = o.uy; Pf[6] = o.zref; Pf[7] = o.alpha; Pf[8] = o.sponge;
      device.queue.writeBuffer(prm, 0, P);
    }
    function dispatch(n) {
      const enc = device.createCommandEncoder();
      const pass = enc.beginComputePass();
      pass.setPipeline(stepPipe);
      for (let k = 0; k < n; k++) { pass.setBindGroup(0, stepBG[cur]); pass.dispatchWorkgroups(...wg); cur = 1 - cur; }
      pass.end();
      device.queue.submit([enc.finish()]);
    }
    return {
      kind: 'webgpu', nx, ny, nz, label: 'WebGPU',
      get lost() { return lost; },
      get busy() { return inflight; },
      lastBatch() { return lastBatch; },
      setHeights(h32) { device.queue.writeBuffer(hBuf, 0, h32); },
      setInflow(ux, uy) { o.ux = ux; o.uy = uy; writeP(false); },
      init() { writeP(true); dispatch(1); writeP(false); },
      // queue n steps unless the previous batch is still running (keeps the GPU from backing up)
      step(n) {
        if (inflight || lost) return 0;
        inflight = true; const t0 = performance.now();
        dispatch(n);
        device.queue.onSubmittedWorkDone().then(() => { inflight = false; lastBatch = { n, ms: performance.now() - t0 }; });
        return n;
      },
      // velocity + density per cell, linear order ((z*ny + y)*nx + x)*4: ux, uy, uz, rho (rho = -1 in solids)
      async readField() {
        if (mapping || lost) return null;
        mapping = true;
        try {
          const enc = device.createCommandEncoder();
          const pass = enc.beginComputePass();
          pass.setPipeline(momPipe); pass.setBindGroup(0, momBG[cur]); pass.dispatchWorkgroups(...wg); pass.end();
          enc.copyBufferToBuffer(vBuf, 0, stage, 0, vBytes);
          device.queue.submit([enc.finish()]);
          await stage.mapAsync(GPUMapMode.READ);
          const out = new Float32Array(stage.getMappedRange().slice(0));
          stage.unmap();
          return out;
        } finally { mapping = false; }
      },
      dispose() { disposed = true; lost = true; try { device.destroy(); } catch (e) { /* already gone */ } },
    };
  }

  // ======================================================================= WebGL2
  const NT = 5;                                    // 19 populations in 5 RGBA32F targets
  const tex = (i) => `T${Math.floor(i / 4)}`, comp = (i) => 'xyzw'[i % 4];
  function glslHead() {
    return `#version 300 es
precision highp float; precision highp int; precision highp sampler2D;
uniform sampler2D T0; uniform sampler2D T1; uniform sampler2D T2; uniform sampler2D T3; uniform sampler2D T4;
uniform sampler2D HM;
uniform int nx; uniform int ny; uniform int nz; uniform int tilesX; uniform int uInit;
uniform vec2 uIn; uniform float zref; uniform float alpha; uniform float sponge;
const float TAU0 = TAU0_V; const float CS2 = CS2_V; const float UMAX = UMAX_V;
ivec2 at(int x, int y, int z) { return ivec2((z % tilesX) * nx + x, (z / tilesX) * ny + y); }
`;
  }
  function glslStep() {
    let s = glslHead() + `
layout(location=0) out vec4 o0; layout(location=1) out vec4 o1; layout(location=2) out vec4 o2; layout(location=3) out vec4 o3; layout(location=4) out vec4 o4;
int topAt(int x, int y) { return int(texelFetch(HM, ivec2(x, y), 0).x + 0.5); }
void main() {
  ivec2 pp = ivec2(gl_FragCoord.xy);
  int tx = pp.x / nx, ty = pp.y / ny;
  int x = pp.x - tx * nx, y = pp.y - ty * ny, z = ty * tilesX + tx;
  if (z >= nz) { o0 = o1 = o2 = o3 = o4 = vec4(0.0); return; }
  vec2 hh = texelFetch(HM, ivec2(x, y), 0).xy;
  int top = int(hh.x + 0.5), gnd = int(hh.y + 0.5);
  float zagl = float(z - gnd) + 0.5;
  float prof = zagl >= zref ? 1.0 : pow(max(zagl, 0.05) / zref, alpha);
  vec3 uin = vec3(uIn * prof, 0.0);
  float uisq = dot(uin, uin);
  if (z < top || uInit == 1) {
    vec3 ue = z < top ? vec3(0.0) : uin; float usq0 = dot(ue, ue);
${C.map((cc, i) => `    float e${i} = ${feqExpr(i, '1.0', 'ue.x', 'ue.y', 'ue.z', 'usq0')};`).join('\n')}
    o0 = vec4(e0, e1, e2, e3); o1 = vec4(e4, e5, e6, e7); o2 = vec4(e8, e9, e10, e11); o3 = vec4(e12, e13, e14, e15); o4 = vec4(e16, e17, e18, 0.0);
    return;
  }
  ivec2 me = pp;
`;
    COLS.forEach(([cx, cy], k) => {
      s += `  int qx${k} = x - (${cx}), qy${k} = y - (${cy});\n`;
      s += `  bool in${k} = qx${k} >= 0 && qx${k} < nx && qy${k} >= 0 && qy${k} < ny;\n`;
      s += `  int t${k} = in${k} ? topAt(qx${k}, qy${k}) : 0;\n`;
    });
    for (let i = 0; i < Q; i++) {
      const [cx, cy, cz] = C[i], k = colOf(i);
      const own = (j) => `texelFetch(${tex(j)}, me, 0).${comp(j)}`;
      s += `  float f${i};\n`;
      s += `  if (!in${k}) { f${i} = ${feqExpr(i, '1.0', 'uin.x', 'uin.y', 'uin.z', 'uisq')}; }\n`;
      if (cz < 0) s += `  else if (z + 1 >= nz) { f${i} = z < t${k} ? ${own(OPP[i])} : texelFetch(${tex(MZ[i])}, at(qx${k}, qy${k}, z), 0).${comp(MZ[i])}; }\n`;
      if (cz > 0) s += `  else if (z - 1 < 0 || z - 1 < t${k}) { f${i} = ${own(OPP[i])}; }\n`;
      else if (cz < 0) s += `  else if (z + 1 < t${k}) { f${i} = ${own(OPP[i])}; }\n`;
      else if (i !== 0) s += `  else if (z < t${k}) { f${i} = ${own(OPP[i])}; }\n`;
      s += `  else { f${i} = texelFetch(${tex(i)}, at(qx${k}, qy${k}, z - (${cz})), 0).${comp(i)}; }\n`;
    }
    s += `  float dd = float(min(min(x, y), min(nx - 1 - x, ny - 1 - y)));
  float spg = 0.0; if (dd < sponge) { float tt = 1.0 - dd / sponge; spg = 0.25 * tt * tt; }
`;
    s += collideBody('glsl');
    s += '  o0 = vec4(g0, g1, g2, g3); o1 = vec4(g4, g5, g6, g7); o2 = vec4(g8, g9, g10, g11); o3 = vec4(g12, g13, g14, g15); o4 = vec4(g16, g17, g18, 0.0);\n}\n';
    return s;
  }
  function glslMoments() {
    return glslHead() + `
out vec4 o;
void main() {
  ivec2 pp = ivec2(gl_FragCoord.xy);
  int tx = pp.x / nx, ty = pp.y / ny;
  int x = pp.x - tx * nx, y = pp.y - ty * ny, z = ty * tilesX + tx;
  int top = int(texelFetch(HM, ivec2(x, y), 0).x + 0.5);
  if (z >= nz || z < top) { o = vec4(0.0, 0.0, 0.0, -1.0); return; }
  vec4 a = texelFetch(T0, pp, 0), b = texelFetch(T1, pp, 0), c = texelFetch(T2, pp, 0), d = texelFetch(T3, pp, 0), e = texelFetch(T4, pp, 0);
  float f[19] = float[19](a.x, a.y, a.z, a.w, b.x, b.y, b.z, b.w, c.x, c.y, c.z, c.w, d.x, d.y, d.z, d.w, e.x, e.y, e.z);
  float rho = 0.0; vec3 m = vec3(0.0);
${C.map((cc, i) => `  rho += f[${i}]; m += vec3(${cc.map((v) => fl(v)).join(', ')}) * f[${i}];`).join('\n')}
  o = vec4(m / rho, rho);
}
`;
  }

  function makeWebGL2(o) {
    const cv = document.createElement('canvas'); cv.width = cv.height = 1;
    let gl;
    try { gl = cv.getContext('webgl2', { antialias: false, depth: false, stencil: false, powerPreference: 'low-power' }); } catch (e) { gl = null; }
    if (!gl || !gl.getExtension('EXT_color_buffer_float')) return null;
    if (gl.getParameter(gl.MAX_DRAW_BUFFERS) < NT || gl.getParameter(gl.MAX_COLOR_ATTACHMENTS) < NT) return null;
    const { nx, ny, nz } = o;
    const maxT = gl.getParameter(gl.MAX_TEXTURE_SIZE);
    let tilesX = Math.ceil(Math.sqrt(nz * ny / nx));
    while (tilesX * nx > maxT && tilesX > 1) tilesX--;
    const tilesY = Math.ceil(nz / tilesX);
    const AW = tilesX * nx, AH = tilesY * ny;
    if (AW > maxT || AH > maxT) return null;
    const sub = (src) => src.replace('TAU0_V', fl(o.tau0)).replace('CS2_V', fl(o.cs2)).replace('UMAX_V', fl(o.umax));
    const sh = (type, src) => {
      const s = gl.createShader(type); gl.shaderSource(s, src); gl.compileShader(s);
      if (!gl.getShaderParameter(s, gl.COMPILE_STATUS)) throw new Error(gl.getShaderInfoLog(s));
      return s;
    };
    const VS = '#version 300 es\nin vec2 p; void main(){ gl_Position = vec4(p, 0.0, 1.0); }';
    const link = (fs) => {
      const pr = gl.createProgram();
      gl.attachShader(pr, sh(gl.VERTEX_SHADER, VS)); gl.attachShader(pr, sh(gl.FRAGMENT_SHADER, sub(fs)));
      gl.bindAttribLocation(pr, 0, 'p'); gl.linkProgram(pr);
      if (!gl.getProgramParameter(pr, gl.LINK_STATUS)) throw new Error(gl.getProgramInfoLog(pr));
      return pr;
    };
    const stepProg = link(glslStep()), momProg = link(glslMoments());
    const vb = gl.createBuffer(); gl.bindBuffer(gl.ARRAY_BUFFER, vb);
    gl.bufferData(gl.ARRAY_BUFFER, new Float32Array([-1, -1, 3, -1, -1, 3]), gl.STATIC_DRAW);
    gl.enableVertexAttribArray(0); gl.vertexAttribPointer(0, 2, gl.FLOAT, false, 0, 0);
    const mk = (w, h, internal, fmt, type) => {
      const t = gl.createTexture(); gl.bindTexture(gl.TEXTURE_2D, t);
      for (const [k, v] of [[gl.TEXTURE_MIN_FILTER, gl.NEAREST], [gl.TEXTURE_MAG_FILTER, gl.NEAREST], [gl.TEXTURE_WRAP_S, gl.CLAMP_TO_EDGE], [gl.TEXTURE_WRAP_T, gl.CLAMP_TO_EDGE]]) gl.texParameteri(gl.TEXTURE_2D, k, v);
      gl.texImage2D(gl.TEXTURE_2D, 0, internal, w, h, 0, fmt, type, null);
      return t;
    };
    const sets = [0, 1].map(() => {
      const t = Array.from({ length: NT }, () => mk(AW, AH, gl.RGBA32F, gl.RGBA, gl.FLOAT));
      const fb = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, fb);
      t.forEach((tx, k) => gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0 + k, gl.TEXTURE_2D, tx, 0));
      gl.drawBuffers(t.map((_, k) => gl.COLOR_ATTACHMENT0 + k));
      return { t, fb, ok: gl.checkFramebufferStatus(gl.FRAMEBUFFER) === gl.FRAMEBUFFER_COMPLETE };
    });
    if (!sets[0].ok || !sets[1].ok) return null;
    const velT = mk(AW, AH, gl.RGBA32F, gl.RGBA, gl.FLOAT);
    const velFB = gl.createFramebuffer(); gl.bindFramebuffer(gl.FRAMEBUFFER, velFB);
    gl.framebufferTexture2D(gl.FRAMEBUFFER, gl.COLOR_ATTACHMENT0, gl.TEXTURE_2D, velT, 0);
    const hm = mk(nx, ny, gl.RG32F, gl.RG, gl.FLOAT);
    const uni = (pr) => { const u = (n) => gl.getUniformLocation(pr, n); return { T: [0, 1, 2, 3, 4].map((k) => u('T' + k)), HM: u('HM'), nx: u('nx'), ny: u('ny'), nz: u('nz'), tilesX: u('tilesX'), uInit: u('uInit'), uIn: u('uIn'), zref: u('zref'), alpha: u('alpha'), sponge: u('sponge') }; };
    const US = uni(stepProg), UM = uni(momProg);
    const setCommon = (pr, u) => {
      gl.useProgram(pr);
      u.T.forEach((l, k) => gl.uniform1i(l, k)); gl.uniform1i(u.HM, 5);
      gl.uniform1i(u.nx, nx); gl.uniform1i(u.ny, ny); gl.uniform1i(u.nz, nz); gl.uniform1i(u.tilesX, tilesX);
      if (u.uIn) gl.uniform2f(u.uIn, o.ux, o.uy);
      if (u.zref) gl.uniform1f(u.zref, o.zref);
      if (u.alpha) gl.uniform1f(u.alpha, o.alpha);
      if (u.sponge) gl.uniform1f(u.sponge, o.sponge);
    };
    setCommon(momProg, UM); setCommon(stepProg, US);
    let cur = 0, lost = false, fence = null, readJob = null;
    cv.addEventListener('webglcontextlost', (e) => { e.preventDefault(); lost = true; });
    gl.viewport(0, 0, AW, AH);
    gl.pixelStorei(gl.UNPACK_ALIGNMENT, 1);
    function bindSrc(set) {
      for (let k = 0; k < NT; k++) { gl.activeTexture(gl.TEXTURE0 + k); gl.bindTexture(gl.TEXTURE_2D, set.t[k]); }
      gl.activeTexture(gl.TEXTURE5); gl.bindTexture(gl.TEXTURE_2D, hm);
    }
    function pass(init) {
      gl.useProgram(stepProg);
      gl.uniform1i(US.uInit, init ? 1 : 0);
      const src = sets[cur], dst = sets[1 - cur];
      gl.bindFramebuffer(gl.FRAMEBUFFER, dst.fb);
      bindSrc(src);
      gl.drawArrays(gl.TRIANGLES, 0, 3);
      cur = 1 - cur;
    }
    const pbo = gl.createBuffer();
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo);
    gl.bufferData(gl.PIXEL_PACK_BUFFER, AW * AH * 16, gl.STREAM_READ);
    gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
    const atlas = new Float32Array(AW * AH * 4);
    const fenceDone = (f) => { const r = gl.clientWaitSync(f, 0, 0); return r === gl.ALREADY_SIGNALED || r === gl.CONDITION_SATISFIED; };
    let lastBatch = { n: 0, ms: 0 }, batchT0 = 0, batchN = 0;
    return {
      kind: 'webgl2', nx, ny, nz, label: 'WebGL2', atlas: [AW, AH, tilesX],
      get lost() { return lost || gl.isContextLost(); },
      get busy() { if (fence && fenceDone(fence)) { gl.deleteSync(fence); fence = null; lastBatch = { n: batchN, ms: performance.now() - batchT0 }; } return !!fence; },
      lastBatch() { return lastBatch; },
      setHeights(h32) {
        const a = new Float32Array(nx * ny * 2);
        for (let i = 0; i < nx * ny; i++) { a[i * 2] = h32[i] & 0xffff; a[i * 2 + 1] = h32[i] >>> 16; }
        gl.activeTexture(gl.TEXTURE5); gl.bindTexture(gl.TEXTURE_2D, hm);
        gl.texSubImage2D(gl.TEXTURE_2D, 0, 0, 0, nx, ny, gl.RG, gl.FLOAT, a);
      },
      setInflow(ux, uy) { o.ux = ux; o.uy = uy; gl.useProgram(stepProg); gl.uniform2f(US.uIn, ux, uy); },
      init() { pass(true); },
      step(n) {
        if (this.busy || this.lost) return 0;
        batchT0 = performance.now(); batchN = n;
        for (let k = 0; k < n; k++) pass(false);
        fence = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0);
        gl.flush();
        return n;
      },
      readField() {
        if (this.lost) return Promise.resolve(null);
        if (readJob) return Promise.resolve(null);
        gl.useProgram(momProg);
        gl.bindFramebuffer(gl.FRAMEBUFFER, velFB);
        bindSrc(sets[cur]);
        gl.drawArrays(gl.TRIANGLES, 0, 3);
        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo);
        gl.readBuffer(gl.COLOR_ATTACHMENT0);
        gl.readPixels(0, 0, AW, AH, gl.RGBA, gl.FLOAT, 0);
        gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
        const f = gl.fenceSync(gl.SYNC_GPU_COMMANDS_COMPLETE, 0); gl.flush();
        readJob = new Promise((resolve) => {
          const poll = () => {
            if (lost) { gl.deleteSync(f); readJob = null; resolve(null); return; }
            if (!fenceDone(f)) { setTimeout(poll, 8); return; }
            gl.deleteSync(f);
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, pbo);
            gl.getBufferSubData(gl.PIXEL_PACK_BUFFER, 0, atlas);
            gl.bindBuffer(gl.PIXEL_PACK_BUFFER, null);
            const out = new Float32Array(nx * ny * nz * 4);       // atlas -> linear ((z*ny + y)*nx + x)*4
            for (let z = 0; z < nz; z++) {
              const ox = (z % tilesX) * nx, oy = Math.floor(z / tilesX) * ny;
              for (let y = 0; y < ny; y++) {
                const s0 = ((oy + y) * AW + ox) * 4;
                out.set(atlas.subarray(s0, s0 + nx * 4), ((z * ny + y) * nx) * 4);
              }
            }
            readJob = null; resolve(out);
          };
          setTimeout(poll, 4);
        });
        return readJob;
      },
      dispose() { const x = gl.getExtension('WEBGL_lose_context'); if (x) x.loseContext(); lost = true; },
    };
  }

  // Chrome on Windows ignores powerPreference for WebGPU and logs a warning when it is passed
  const isWin = () => /Windows/.test(navigator.userAgent);
  const gpuOpts = () => (isWin() ? {} : { powerPreference: 'low-power' });
  const VENDORS = ['intel', 'nvidia', 'amd', 'apple', 'qualcomm', 'arm', 'microsoft'];
  const vendorOf = (t) => { t = (t || '').toLowerCase(); if (t.includes('radeon') || t.includes('ati ')) return 'amd'; if (t.includes('geforce') || t.includes('quadro') || t.includes('rtx')) return 'nvidia'; if (t.includes('adreno')) return 'qualcomm'; if (t.includes('mali')) return 'arm'; return VENDORS.find((v) => t.includes(v)) || ''; };
  async function gpuVendors() {
    const out = { gpu: '', gl: '' };
    if (navigator.gpu) {
      const a = await navigator.gpu.requestAdapter(gpuOpts());
      if (a) { const i = a.info || (a.requestAdapterInfo ? await a.requestAdapterInfo() : {}); out.gpu = vendorOf(i.vendor) || vendorOf(i.description); }
    }
    const cv = document.createElement('canvas'); cv.width = cv.height = 1;
    const gl = cv.getContext('webgl2', { powerPreference: 'low-power' });
    if (gl) {
      const d = gl.getExtension('WEBGL_debug_renderer_info');
      out.gl = vendorOf(d ? gl.getParameter(d.UNMASKED_RENDERER_WEBGL) : gl.getParameter(gl.RENDERER));
      const x = gl.getExtension('WEBGL_lose_context'); if (x) x.loseContext();
    }
    return out;
  }

  /* create({nx, ny, nz, ux, uy, zref, alpha, sponge, tau0, cs2, umax, prefer}) -> solver or null.
     prefer: 'auto' (WebGPU, then WebGL2), 'webgpu', 'webgl2'. */
  async function create(o) {
    const opt = Object.assign({ tau0: 0.5 + 3 * 0.005, cs2: 0.17 * 0.17, umax: 0.25, ux: 0, uy: 0, zref: 40, alpha: 0.25, sponge: 6 }, o);
    let tries = opt.prefer === 'webgl2' ? ['webgl2'] : opt.prefer === 'webgpu' ? ['webgpu'] : ['webgpu', 'webgl2'];
    // Chrome on Windows ignores WebGPU's powerPreference, so on a laptop with two GPUs WebGPU can land on the
    // big one. WebGL2 does honour 'low-power'. If the two disagree about the vendor, stay on WebGL2.
    if (opt.prefer === 'auto' || !opt.prefer) {
      try {
        const v = await gpuVendors();
        opt.vendors = v;
        if (v.gpu && v.gl && v.gpu !== v.gl) tries = ['webgl2', 'webgpu'];
      } catch (e) { /* keep the default order */ }
    }
    let lastErr = null;
    for (const k of tries) {
      try {
        const s = k === 'webgpu' ? await makeWebGPU(Object.assign({}, opt)) : makeWebGL2(Object.assign({}, opt));
        if (s) { s.vendors = opt.vendors || null; return s; }
      } catch (e) { lastErr = e; console.warn('[wind3d] ' + k + ' solver failed:', e.message); }
    }
    if (lastErr) console.info('[wind3d] no 3D solver on this device');
    return null;
  }

  window.HKWSolver3D = { create, C, W, OPP, MZ, _src: { wgslStep, wgslMoments, glslStep, glslMoments } };
})();
