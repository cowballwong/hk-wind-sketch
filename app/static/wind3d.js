/* HK Wind Sketch v3: 3D view (ES module, loaded only when someone opens 3D).
 *
 * three.js r160 from a pinned CDN (see the import map in the page). The ground is
 * a terrain mesh from the LandsD digital terrain model, built like a site model
 * (a plinth with cut sides). Buildings are extruded from LandsD footprints and sit
 * on their own base level; the bottom is taken down to the lowest ground under the
 * footprint so nothing floats on a slope. The coloured skin on the ground is the
 * pedestrian-level wind field (about 2 m above the terrain).
 *
 * When the 3D solver runs, the page also passes 3D streamlines, moving particles
 * and an optional vertical section; with the 2D fallback only the draped field and
 * flat streaks are shown.
 *
 * World axes: X = east, Y = up, Z = south (north is -Z). Page local metres
 * (x east, y south, origin top-left of the square) map to X = x - D/2, Z = y - D/2,
 * and Y = level (mPD) - zRef, where zRef is the lowest ground in the square.
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';

const TAIL = 14;              // points per particle trail
const PLINTH = 14;            // depth of the site-model base below the lowest ground (m)

export function create(opts) {
  const { container, fieldCanvas, domain: D, onPick } = opts;
  const H = D / 2;
  const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'low-power', alpha: false });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.5));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.shadowMap.autoUpdate = false;
  container.appendChild(renderer.domElement);
  renderer.domElement.className = 'hkw-3dcanvas';

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(32, 1, 5, 9000);
  const HOME = [520, 600, 800];
  camera.position.set(...HOME);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 0, 30);
  controls.enableDamping = true; controls.dampingFactor = 0.09;
  controls.maxPolarAngle = 1.5; controls.minDistance = 60; controls.maxDistance = 2800;
  controls.screenSpacePanning = true;
  let dirty = true;
  controls.addEventListener('change', () => { dirty = true; });

  const hemi = new THREE.HemisphereLight(0xf3f6f9, 0x8d8b84, 1.6);
  scene.add(hemi);
  const sun = new THREE.DirectionalLight(0xfff6ea, 2.2);
  sun.position.set(-330, 560, 420);
  sun.castShadow = true;
  sun.shadow.mapSize.set(2048, 2048);
  Object.assign(sun.shadow.camera, { left: -H - 90, right: H + 90, top: H + 90, bottom: -H - 90, near: 50, far: 2600 });
  sun.shadow.bias = -0.0004; sun.shadow.normalBias = 0.6;
  scene.add(sun); scene.add(sun.target);

  // ---------------------------------------------------------------- terrain
  const T = { n: 2, step: D, z: new Float32Array(4), zRef: 0, zMax: 0 };
  function groundY(x, y) {           // local metres -> Y (level above zRef), bilinear on the terrain grid
    const g = T, fx = Math.min(g.n - 1.0001, Math.max(0, x / g.step)), fy = Math.min(g.n - 1.0001, Math.max(0, y / g.step));
    const i = Math.floor(fx), j = Math.floor(fy), tx = fx - i, ty = fy - j, n = g.n;
    const z = g.z[j * n + i] * (1 - tx) * (1 - ty) + g.z[j * n + i + 1] * tx * (1 - ty) + g.z[(j + 1) * n + i] * (1 - tx) * ty + g.z[(j + 1) * n + i + 1] * tx * ty;
    return z - g.zRef;
  }
  const fieldTex = new THREE.CanvasTexture(fieldCanvas);
  fieldTex.colorSpace = THREE.SRGBColorSpace;
  fieldTex.minFilter = THREE.LinearFilter; fieldTex.magFilter = THREE.LinearFilter; fieldTex.generateMipmaps = false;
  // lit for the relief, plus some of the field as glow so the colours stay true to the legend
  const groundMat = new THREE.MeshStandardMaterial({ color: 0xffffff, map: fieldTex, emissive: 0xffffff, emissiveMap: fieldTex, emissiveIntensity: 0.42, roughness: 1, metalness: 0 });
  const plainGround = new THREE.Color(0xd9d6cf);
  const sideMat = new THREE.MeshStandardMaterial({ color: 0xc9c4b9, roughness: 1 });
  const contourMat = new THREE.LineBasicMaterial({ color: 0x2c2b28, transparent: true, opacity: 0.22 });
  let ground = null, sides = null, rim = null, contours = null;
  function setTerrain(t) {
    Object.assign(T, t);
    for (const o of [ground, sides, rim, contours]) if (o) { scene.remove(o); o.geometry.dispose(); }
    const n = T.n, segs = n - 1;
    const g = new THREE.PlaneGeometry(D, D, segs, segs); g.rotateX(-Math.PI / 2);
    const p = g.attributes.position;
    for (let k = 0; k < p.count; k++) {
      const i = k % n, j = Math.floor(k / n);      // PlaneGeometry rows run from -Z (north) to +Z after rotation
      p.setY(k, T.z[j * n + i] - T.zRef);
    }
    g.computeVertexNormals();
    ground = new THREE.Mesh(g, groundMat); ground.receiveShadow = true; scene.add(ground);
    // cut sides of the site model
    const sv = [], base = -PLINTH;
    const edge = (pts) => { for (let k = 0; k < pts.length - 1; k++) { const [a, b] = [pts[k], pts[k + 1]]; sv.push(a[0], a[1], a[2], b[0], b[1], b[2], b[0], base, b[2], a[0], a[1], a[2], b[0], base, b[2], a[0], base, a[2]); } };
    const P = (i, j) => [i * T.step - H, T.z[j * n + i] - T.zRef, j * T.step - H];
    const N = [], S2 = [], Wv = [], E = [];
    for (let i = 0; i < n; i++) { N.push(P(n - 1 - i, 0)); S2.push(P(i, n - 1)); }
    for (let j = 0; j < n; j++) { Wv.push(P(0, j)); E.push(P(n - 1, n - 1 - j)); }
    edge(N); edge(S2); edge(Wv); edge(E);
    const sg = new THREE.BufferGeometry(); sg.setAttribute('position', new THREE.Float32BufferAttribute(sv, 3)); sg.computeVertexNormals();
    sideMat.side = THREE.DoubleSide;
    sides = new THREE.Mesh(sg, sideMat); scene.add(sides);
    const rv = [];
    for (const L of [N, S2, Wv, E]) for (let k = 0; k < L.length - 1; k++) rv.push(...L[k], ...L[k + 1]);
    for (const [i, j] of [[0, 0], [n - 1, 0], [0, n - 1], [n - 1, n - 1]]) { const q = P(i, j); rv.push(q[0], q[1], q[2], q[0], base, q[2]); }
    const rg = new THREE.BufferGeometry(); rg.setAttribute('position', new THREE.Float32BufferAttribute(rv, 3));
    rim = new THREE.LineSegments(rg, new THREE.LineBasicMaterial({ color: 0x3b3a36, transparent: true, opacity: 0.6 })); scene.add(rim);
    contours = makeContours(); if (contours) scene.add(contours);
    // default camera: on the downhill side of a sloping site (looking up the slope), else from the south-east
    let gx = 0, gz = 0;
    for (let j = 0; j < n; j++) for (let i = 0; i < n; i++) { const w = T.z[j * n + i] - (T.zMean || T.zRef); gx += w * (i * T.step - H); gz += w * (j * T.step - H); }
    const relief = T.zMax - T.zRef;
    home.az = relief > 30 && Math.hypot(gx, gz) > 0 ? Math.atan2(-gx, -gz) : 0.576;
    home.polar = relief > 30 ? 0.9 : 0.97;
    home.dist = 1150 + Math.min(500, relief * 0.8);
    home.ty = (T.zMean || T.zRef) - T.zRef;
    goHome();
    renderer.shadowMap.needsUpdate = true; dirty = true;
  }
  const home = { az: 0.576, polar: 0.97, dist: 1150, ty: 0 };
  function goHome() {
    controls.target.set(0, home.ty, 0);
    camera.position.set(home.dist * Math.sin(home.polar) * Math.sin(home.az), home.ty + home.dist * Math.cos(home.polar), home.dist * Math.sin(home.polar) * Math.cos(home.az));
    controls.update(); dirty = true;
  }
  // 10 m contours by marching squares (drawn faintly, like a site model's terraces)
  function makeContours() {
    const n = T.n, v = [];
    const lo = Math.ceil(T.zRef / 10) * 10, hi = T.zMax;
    if (hi - T.zRef < 8) return null;
    for (let lev = lo; lev <= hi; lev += 10) {
      for (let j = 0; j < n - 1; j++) for (let i = 0; i < n - 1; i++) {
        const c = [[i, j], [i + 1, j], [i + 1, j + 1], [i, j + 1]].map(([a, b]) => [a, b, T.z[b * n + a]]);
        const pts = [];
        for (let k = 0; k < 4; k++) {
          const a = c[k], b = c[(k + 1) % 4];
          if ((a[2] < lev) !== (b[2] < lev)) { const t = (lev - a[2]) / (b[2] - a[2]); pts.push([(a[0] + (b[0] - a[0]) * t) * T.step - H, lev - T.zRef + 0.25, (a[1] + (b[1] - a[1]) * t) * T.step - H]); }
        }
        if (pts.length >= 2) v.push(...pts[0], ...pts[1]);
        if (pts.length === 4) v.push(...pts[2], ...pts[3]);
      }
    }
    if (!v.length) return null;
    const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(v, 3));
    return new THREE.LineSegments(g, contourMat);
  }

  // ---------------------------------------------------------------- buildings
  const bldMat = new THREE.MeshStandardMaterial({ color: 0xf6f5f1, roughness: 0.8, metalness: 0 });
  const ghostMat = new THREE.MeshStandardMaterial({ color: 0xf4f2ed, roughness: 0.9, transparent: true, opacity: 0.35, depthWrite: false });
  const userMat = new THREE.MeshStandardMaterial({ color: 0xd97757, roughness: 0.7 });
  const selMat = new THREE.MeshStandardMaterial({ color: 0xe8956f, roughness: 0.7, emissive: 0x3a1408 });
  const edgeMat = new THREE.LineBasicMaterial({ color: 0x2c2b28, transparent: true, opacity: 0.5 });
  const removedMat = new THREE.LineDashedMaterial({ color: 0xc0392b, dashSize: 4, gapSize: 3 });
  const geoCache = new Map();
  let solid = null, ghosts = null, edges = null, removedLines = null;
  const userGroup = new THREE.Group();
  let ranges = [];
  scene.add(userGroup);

  // bottom: the lower of the LandsD base level and the lowest ground under the footprint
  function seat(rings, base) {
    let gmin = Infinity;
    for (const r of rings) for (const [x, y] of r) gmin = Math.min(gmin, groundY(x, y));
    const b = base == null ? gmin : Math.min(base - T.zRef, gmin);
    return isFinite(b) ? b - 0.3 : 0;
  }
  function ringShape(ring) {
    const s = new THREE.Shape();
    ring.forEach(([x, y], i) => (i ? s.lineTo(x - H, -(y - H)) : s.moveTo(x - H, -(y - H))));
    return s;
  }
  function extrudeRings(rings, y0, y1) {
    const shapes = rings.filter((r) => r.length >= 3).map(ringShape);
    if (!shapes.length) return null;
    const g = new THREE.ExtrudeGeometry(shapes, { depth: Math.max(0.6, y1 - y0), bevelEnabled: false, curveSegments: 1 });
    g.rotateX(-Math.PI / 2); g.translate(0, y0, 0); g.deleteAttribute('uv');
    return g;
  }
  function extrude(b) {
    if (geoCache.has(b.id)) return geoCache.get(b.id);
    let g = null, y0 = 0, y1 = 0;
    try {
      y0 = seat(b.rings, b.base);
      y1 = b.top != null ? b.top - T.zRef : y0 + b.h;
      g = extrudeRings(b.rings, y0, y1);
    } catch (e) { g = null; }
    const r = g ? { g, y0, y1 } : null;
    geoCache.set(b.id, r);
    return r;
  }
  function pushEdges(arr, ring, y0, y1) {
    const closed = ring[0][0] === ring[ring.length - 1][0] && ring[0][1] === ring[ring.length - 1][1];
    const n = closed ? ring.length - 1 : ring.length;
    for (let i = 0; i < n; i++) {
      const [x0, z0] = ring[i], [x1, z1] = ring[(i + 1) % n];
      arr.push(x0 - H, y1, z0 - H, x1 - H, y1, z1 - H, x0 - H, y0 + 0.3, z0 - H, x1 - H, y0 + 0.3, z1 - H);
      const [xp, zp] = ring[(i - 1 + n) % n];
      const a1 = Math.atan2(z0 - zp, x0 - xp), a2 = Math.atan2(z1 - z0, x1 - x0);
      let d = Math.abs(a2 - a1); if (d > Math.PI) d = 2 * Math.PI - d;
      if (d > 0.35) arr.push(x0 - H, y0 + 0.3, z0 - H, x0 - H, y1, z0 - H);
    }
  }
  function disposeObj(o) { if (o) { scene.remove(o); o.geometry.dispose(); } }

  function rebuild(st, newSite) {
    if (newSite) { for (const r of geoCache.values()) if (r) r.g.dispose(); geoCache.clear(); }
    disposeObj(solid); disposeObj(ghosts); disposeObj(edges); disposeObj(removedLines);
    solid = ghosts = edges = removedLines = null;
    const sg = [], gg = [], ev = [], rv = [];
    ranges = []; let v = 0;
    for (const b of st.buildings) {
      const r = extrude(b); if (!r) continue;
      if (st.removed.has(b.id)) { for (const ring of b.rings) pushEdges(rv, ring, r.y0, r.y1); continue; }
      if (st.blocks(b)) { sg.push(r.g); ranges.push([v, b]); v += r.g.attributes.position.count; for (const ring of b.rings) pushEdges(ev, ring, r.y0, r.y1); }
      else gg.push(r.g);
    }
    if (sg.length) { solid = new THREE.Mesh(mergeGeometries(sg, false), bldMat); solid.castShadow = solid.receiveShadow = true; scene.add(solid); }
    if (gg.length) { ghosts = new THREE.Mesh(mergeGeometries(gg, false), ghostMat); scene.add(ghosts); }
    const lg = (arr, mat) => { const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3)); const l = new THREE.LineSegments(g, mat); scene.add(l); return l; };
    if (ev.length) edges = lg(ev, edgeMat);
    if (rv.length) { removedLines = lg(rv, removedMat); removedLines.computeLineDistances(); }
    renderer.shadowMap.needsUpdate = true; dirty = true;
  }

  // user buildings: free polygons {poly, h, base}
  function rebuildUser(st) {
    for (const c of userGroup.children.slice()) { userGroup.remove(c); c.geometry.dispose(); c.children.forEach((e) => e.geometry.dispose()); }
    for (const u of st.user) {
      const y0 = (u.base != null ? u.base - T.zRef : seat([u.poly], null)) - 0.3;
      const y1 = y0 + 0.3 + u.h;
      const g = extrudeRings([u.poly], y0, y1); if (!g) continue;
      const m = new THREE.Mesh(g, u === st.sel ? selMat : userMat);
      m.castShadow = m.receiveShadow = true; m.userData.u = u;
      const ev = []; pushEdges(ev, u.poly, y0, y1);
      const eg = new THREE.BufferGeometry(); eg.setAttribute('position', new THREE.Float32BufferAttribute(ev, 3));
      m.add(new THREE.LineSegments(eg, edgeMat));
      userGroup.add(m);
    }
    renderer.shadowMap.needsUpdate = true; dirty = true;
  }

  // ---------------------------------------------------------------- wind arrow
  const arrow = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 110, 0x1d3f6e, 34, 20);
  scene.add(arrow);
  const label = makeLabel('Wind 風');
  scene.add(label);
  let windIdx = 0;
  function setWind(dirIdx) {
    windIdx = dirIdx;
    const b = dirIdx * 22.5 * Math.PI / 180;
    const down = new THREE.Vector3(-Math.sin(b), 0, Math.cos(b));
    const start = down.clone().multiplyScalar(-(H + 150));
    start.y = Math.max(40, (T.zMax - T.zRef) * 0.5 + 40);
    arrow.position.copy(start); arrow.setDirection(down);
    label.position.copy(start).add(new THREE.Vector3(0, 30, 0));
    dirty = true;
  }
  function makeLabel(text) {
    const c = document.createElement('canvas'); c.width = 256; c.height = 64;
    const x = c.getContext('2d'); x.font = '600 34px ui-sans-serif,system-ui,"Microsoft JhengHei",sans-serif';
    x.fillStyle = '#1d3f6e'; x.textAlign = 'center'; x.textBaseline = 'middle'; x.fillText(text, 128, 34);
    const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace;
    const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: t, depthTest: false, transparent: true }));
    s.scale.set(80, 20, 1); s.renderOrder = 5; return s;
  }

  // ---------------------------------------------------------------- streamlines (3D solver)
  const lineMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.9, depthWrite: false, toneMapped: false });
  let lines = null;
  function setStreamlines(pos, col, show) {
    if (lines) { scene.remove(lines); lines.geometry.dispose(); lines = null; }
    if (pos && pos.length && show) {
      const g = new THREE.BufferGeometry();
      g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
      g.setAttribute('color', new THREE.BufferAttribute(col, 4));
      lines = new THREE.LineSegments(g, lineMat); lines.frustumCulled = false; lines.renderOrder = 3;
      scene.add(lines);
    }
    dirty = true;
  }

  // ---------------------------------------------------------------- particles (3D) or flat streaks (2D fallback)
  let hist = null, nPart = 0, streaks = null;
  const streakMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, depthWrite: false, toneMapped: false });
  function ensureStreaks(n) {
    if (n === nPart && streaks) return;
    if (streaks) { scene.remove(streaks); streaks.geometry.dispose(); }
    nPart = n; hist = new Float32Array(n * TAIL * 3);
    const segs = n * (TAIL - 1);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(segs * 2 * 3), 3).setUsage(THREE.DynamicDrawUsage));
    g.setAttribute('color', new THREE.BufferAttribute(new Float32Array(segs * 2 * 4), 4).setUsage(THREE.DynamicDrawUsage));
    streaks = new THREE.LineSegments(g, streakMat); streaks.frustumCulled = false; streaks.renderOrder = 4;
    scene.add(streaks);
  }
  // parts: {x, y, z?, fresh, c:[r,g,b]} ; z = height above zRef (3D) or undefined (flat 2D streaks just above the ground)
  function pushParticles(parts, show) {
    ensureStreaks(parts.length);
    streaks.visible = show;
    if (!show) return;
    for (let i = 0; i < parts.length; i++) {
      const p = parts[i], o = i * TAIL * 3;
      const zz = p.z != null ? p.z : groundY(p.x, p.y) + 3;
      if (p.fresh) { for (let k = 0; k < TAIL; k++) { hist[o + k * 3] = p.x; hist[o + k * 3 + 1] = p.y; hist[o + k * 3 + 2] = zz; } p.fresh = false; continue; }
      hist.copyWithin(o + 3, o, o + (TAIL - 1) * 3);
      hist[o] = p.x; hist[o + 1] = p.y; hist[o + 2] = zz;
    }
    const pos = streaks.geometry.attributes.position.array, col = streaks.geometry.attributes.color.array;
    let w = 0, cw = 0;
    for (let i = 0; i < parts.length; i++) {
      const o = i * TAIL * 3, c = parts[i].c || [1, 1, 1];
      for (let k = 0; k < TAIL - 1; k++) {
        pos[w++] = hist[o + k * 3] - H; pos[w++] = hist[o + k * 3 + 2]; pos[w++] = hist[o + k * 3 + 1] - H;
        pos[w++] = hist[o + k * 3 + 3] - H; pos[w++] = hist[o + k * 3 + 5]; pos[w++] = hist[o + k * 3 + 4] - H;
        const a0 = 0.95 * (1 - k / (TAIL - 1)), a1 = 0.95 * (1 - (k + 1) / (TAIL - 1));
        col[cw++] = c[0]; col[cw++] = c[1]; col[cw++] = c[2]; col[cw++] = a0;
        col[cw++] = c[0]; col[cw++] = c[1]; col[cw++] = c[2]; col[cw++] = a1;
      }
    }
    streaks.geometry.attributes.position.needsUpdate = true;
    streaks.geometry.attributes.color.needsUpdate = true;
    dirty = true;
  }

  // ---------------------------------------------------------------- vertical section
  let section = null, sectionTex = null, sectionEdge = null;
  function setSection(s) {
    if (section) { scene.remove(section); section.geometry.dispose(); section = null; }
    if (sectionEdge) { scene.remove(sectionEdge); sectionEdge.geometry.dispose(); sectionEdge = null; }
    if (!s) { dirty = true; return; }
    if (!sectionTex || sectionTex.image !== s.canvas || sectionTex.userData.w !== s.canvas.width || sectionTex.userData.h !== s.canvas.height) {
      if (sectionTex) sectionTex.dispose();
      sectionTex = new THREE.CanvasTexture(s.canvas); sectionTex.colorSpace = THREE.SRGBColorSpace;
      sectionTex.minFilter = THREE.LinearFilter; sectionTex.magFilter = THREE.LinearFilter; sectionTex.generateMipmaps = false;
      sectionTex.userData.w = s.canvas.width; sectionTex.userData.h = s.canvas.height;
    }
    sectionTex.needsUpdate = true;
    const [ax, ay] = s.a, [bx, by] = s.b, y0 = s.y0, y1 = s.y1;
    const g = new THREE.BufferGeometry();
    const P = [ax - H, y0, ay - H, bx - H, y0, by - H, bx - H, y1, by - H, ax - H, y1, ay - H];
    g.setAttribute('position', new THREE.Float32BufferAttribute(P, 3));
    g.setAttribute('uv', new THREE.Float32BufferAttribute([0, 0, 1, 0, 1, 1, 0, 1], 2));
    g.setIndex([0, 1, 2, 0, 2, 3]);
    section = new THREE.Mesh(g, new THREE.MeshBasicMaterial({ map: sectionTex, transparent: true, side: THREE.DoubleSide, depthWrite: false, toneMapped: false, opacity: 0.92 }));
    section.renderOrder = 2; scene.add(section);
    const eg = new THREE.BufferGeometry(); eg.setAttribute('position', new THREE.Float32BufferAttribute([...P.slice(0, 3), ...P.slice(3, 6), ...P.slice(3, 6), ...P.slice(6, 9), ...P.slice(6, 9), ...P.slice(9, 12), ...P.slice(9, 12), ...P.slice(0, 3)], 3));
    sectionEdge = new THREE.LineSegments(eg, new THREE.LineBasicMaterial({ color: 0x1d3f6e, transparent: true, opacity: 0.7 }));
    scene.add(sectionEdge);
    dirty = true;
  }

  function setTheme(dark) {
    const bg = dark ? 0x22221f : 0xeef1f4;
    scene.background = new THREE.Color(bg);
    scene.fog = new THREE.Fog(bg, 1600, 5000);
    bldMat.color.set(dark ? 0xe2ded6 : 0xf6f5f1);
    sideMat.color.set(dark ? 0x4a4843 : 0xc9c4b9);
    plainGround.set(dark ? 0x3d3b37 : 0xd9d6cf);
    if (!groundMat.map) groundMat.color.copy(plainGround);
    dirty = true;
  }

  // ---------------------------------------------------------------- picking
  const ray = new THREE.Raycaster(), ndc = new THREE.Vector2();
  let down = null;
  renderer.domElement.addEventListener('pointerdown', (e) => { down = { x: e.clientX, y: e.clientY, t: performance.now() }; });
  renderer.domElement.addEventListener('pointerup', (e) => {
    if (!down) return;
    const moved = Math.hypot(e.clientX - down.x, e.clientY - down.y), dt = performance.now() - down.t;
    down = null;
    if (moved > 6 || dt > 600 || !onPick) return;
    const r = renderer.domElement.getBoundingClientRect();
    ndc.set(((e.clientX - r.left) / r.width) * 2 - 1, -((e.clientY - r.top) / r.height) * 2 + 1);
    ray.setFromCamera(ndc, camera);
    const ui = ray.intersectObjects(userGroup.children, false)[0];
    const hit = solid ? ray.intersectObject(solid, false)[0] : null;
    const gh = ground ? ray.intersectObject(ground, false)[0] : null;
    const best = [ui && { k: 'user', d: ui.distance, h: ui }, hit && { k: 'building', d: hit.distance, h: hit }, gh && { k: 'ground', d: gh.distance, h: gh }]
      .filter(Boolean).sort((a, b) => a.d - b.d)[0];
    if (!best) return;
    if (best.k === 'user') return onPick({ kind: 'user', u: best.h.object.userData.u, ev: e });
    if (best.k === 'building') {
      const vi = best.h.face.a; let lo = 0, hi = ranges.length - 1;
      while (lo < hi) { const m = (lo + hi + 1) >> 1; if (ranges[m][0] <= vi) lo = m; else hi = m - 1; }
      return onPick({ kind: 'building', b: ranges[lo][1], ev: e });
    }
    const gp = best.h.point;
    if (Math.abs(gp.x) <= H && Math.abs(gp.z) <= H) return onPick({ kind: 'ground', x: gp.x + H, y: gp.z + H, ev: e });
  });

  function resize(w, h) {
    renderer.setSize(w, h, false);
    renderer.domElement.style.width = w + 'px'; renderer.domElement.style.height = h + 'px';
    camera.aspect = w / h; camera.updateProjectionMatrix(); dirty = true;
  }
  let fieldDirty = false;
  return {
    THREE,
    setTerrain, groundY, rebuild, rebuildUser, setWind, pushParticles, setStreamlines, setSection, setTheme, resize,
    markField() { fieldDirty = true; dirty = true; },
    showField(on) {
      groundMat.map = on ? fieldTex : null; groundMat.emissiveMap = on ? fieldTex : null; groundMat.emissiveIntensity = on ? 0.42 : 0;
      groundMat.color.copy(on ? new THREE.Color(0xffffff) : plainGround);
      groundMat.needsUpdate = true; dirty = true;
    },
    frame() {
      controls.update();
      if (fieldDirty) {
        // a new site can change the grid size: let three.js allocate the texture again instead of copying into the old size
        if (fieldTex.userData.w !== fieldCanvas.width || fieldTex.userData.h !== fieldCanvas.height) { fieldTex.dispose(); fieldTex.userData.w = fieldCanvas.width; fieldTex.userData.h = fieldCanvas.height; }
        fieldTex.needsUpdate = true; fieldDirty = false;
      }
      if (!dirty) return false;
      renderer.render(scene, camera); dirty = false; return true;
    },
    azimuth() { return controls.getAzimuthalAngle(); },
    camera, controls, renderer,
    stats() {
      const i = renderer.info;
      return { calls: i.render.calls, triangles: i.render.triangles, solid: !!solid, user: userGroup.children.length,
        lines: lines ? lines.geometry.attributes.position.count / 2 : 0, section: !!section, terrainN: T.n,
        cam: camera.position.toArray().map((v) => +v.toFixed(1)) };
    },
    // where a building's extrusion starts and ends (for the seating check)
    seatOf(b) { const r = extrude(b); return r ? { y0: r.y0 + T.zRef, y1: r.y1 + T.zRef } : null; },
    resetView() { goHome(); },
    homeAzimuth() { return home.az; },
    orbitTo(az, polar, dist) {
      const t = controls.target, r = dist || camera.position.distanceTo(t);
      camera.position.set(t.x + r * Math.sin(polar) * Math.sin(az), t.y + r * Math.cos(polar), t.z + r * Math.sin(polar) * Math.cos(az));
      controls.update(); dirty = true;
    },
    dispose() { renderer.dispose(); renderer.forceContextLoss(); renderer.domElement.remove(); },
  };
}
