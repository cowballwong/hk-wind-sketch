/* HK Wind Sketch: 3D view (ES module, loaded only when someone opens 3D).
 *
 * three.js r160 from a pinned CDN (see the import map in the page). Buildings are
 * extruded from LandsD footprints to their real height (TopHeight - BaseHeight),
 * standing on a flat site model: the sketch has no terrain. The coloured ground
 * layer is the SAME 2D slice the plan view shows, draped at about 2 m; the 3D view
 * does not run a 3D simulation and changes nothing in the numbers.
 *
 * World axes: X = east, Y = up, Z = south (so north is -Z). The page's local metres
 * (x east, y south, origin top-left of the square) map to X = x - D/2, Z = y - D/2.
 */
import * as THREE from 'three';
import { OrbitControls } from 'three/addons/controls/OrbitControls.js';
import { mergeGeometries } from 'three/addons/utils/BufferGeometryUtils.js';

const TAIL = 11;            // points per streak
const STREAK_Y = 3.0;      // streaks float just above the 2 m slice
const FIELD_Y = 2.0;

export function create(opts) {
  const { container, fieldCanvas, domain: D, onPick } = opts;
  const H = D / 2;
  const renderer = new THREE.WebGLRenderer({ antialias: true, powerPreference: 'low-power', alpha: false });
  renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 1.5));
  renderer.outputColorSpace = THREE.SRGBColorSpace;
  renderer.shadowMap.enabled = true;
  renderer.shadowMap.type = THREE.PCFSoftShadowMap;
  renderer.shadowMap.autoUpdate = false;            // the sun does not move: redraw shadows only on edits
  container.appendChild(renderer.domElement);
  renderer.domElement.className = 'hkw-3dcanvas';

  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(32, 1, 5, 8000);
  const HOME = [500, 640, 830];
  camera.position.set(...HOME);
  const controls = new OrbitControls(camera, renderer.domElement);
  controls.target.set(0, 0, 30);
  controls.enableDamping = true; controls.dampingFactor = 0.09;
  controls.maxPolarAngle = 1.45; controls.minDistance = 90; controls.maxDistance = 2600;
  controls.screenSpacePanning = true;
  let dirty = true;
  controls.addEventListener('change', () => { dirty = true; });

  // light: soft sky plus a low-ish afternoon sun from the south-west
  const hemi = new THREE.HemisphereLight(0xf3f5f7, 0x8d877c, 1.55);
  scene.add(hemi);
  const sun = new THREE.DirectionalLight(0xfff3e2, 2.3);
  sun.position.set(-330, 520, 420);
  sun.castShadow = true;
  sun.shadow.mapSize.set(2048, 2048);
  Object.assign(sun.shadow.camera, { left: -H - 80, right: H + 80, top: H + 80, bottom: -H - 80, near: 50, far: 2000 });
  sun.shadow.bias = -0.0004; sun.shadow.normalBias = 0.6;
  scene.add(sun);

  // site model: plinth, its edge, and a shadow-catcher over the field
  const plinthMat = new THREE.MeshStandardMaterial({ color: 0xdedad2, roughness: 0.95 });
  const plinth = new THREE.Mesh(new THREE.BoxGeometry(D + 8, 8, D + 8), plinthMat);
  plinth.position.y = -4; plinth.receiveShadow = true; scene.add(plinth);
  const frame = new THREE.LineSegments(new THREE.EdgesGeometry(new THREE.BoxGeometry(D + 8, 8, D + 8)),
    new THREE.LineBasicMaterial({ color: 0x3b3a36, transparent: true, opacity: 0.6 }));
  frame.position.y = -4; scene.add(frame);

  const fieldTex = new THREE.CanvasTexture(fieldCanvas);
  fieldTex.colorSpace = THREE.SRGBColorSpace;
  fieldTex.minFilter = THREE.LinearFilter; fieldTex.magFilter = THREE.LinearFilter; fieldTex.generateMipmaps = false;
  const fieldMat = new THREE.MeshBasicMaterial({ map: fieldTex, transparent: true, toneMapped: false, depthWrite: false });
  const fieldGeo = new THREE.PlaneGeometry(D, D); fieldGeo.rotateX(-Math.PI / 2);
  const field = new THREE.Mesh(fieldGeo, fieldMat); field.position.y = FIELD_Y; field.renderOrder = 1; scene.add(field);
  const shadowCatch = new THREE.Mesh(fieldGeo, new THREE.ShadowMaterial({ opacity: 0.2, depthWrite: false }));
  shadowCatch.position.y = FIELD_Y + 0.05; shadowCatch.receiveShadow = true; shadowCatch.renderOrder = 2; scene.add(shadowCatch);
  const groundGeo = new THREE.PlaneGeometry(D, D); groundGeo.rotateX(-Math.PI / 2);
  const groundMat = new THREE.MeshStandardMaterial({ color: 0xcfcac0, roughness: 1 });
  const ground = new THREE.Mesh(groundGeo, groundMat); ground.position.y = 0.02; ground.receiveShadow = true; scene.add(ground);

  // buildings
  const bldMat = new THREE.MeshStandardMaterial({ color: 0xf4f2ed, roughness: 0.82, metalness: 0 });
  const ghostMat = new THREE.MeshStandardMaterial({ color: 0xf4f2ed, roughness: 0.9, transparent: true, opacity: 0.35, depthWrite: false });
  const userMat = new THREE.MeshStandardMaterial({ color: 0xd97757, roughness: 0.7 });
  const selMat = new THREE.MeshStandardMaterial({ color: 0xe8956f, roughness: 0.7, emissive: 0x3a1408 });
  const edgeMat = new THREE.LineBasicMaterial({ color: 0x2c2b28, transparent: true, opacity: 0.55 });
  const removedMat = new THREE.LineDashedMaterial({ color: 0xc0392b, dashSize: 4, gapSize: 3 });
  const geoCache = new Map();      // building id -> extruded geometry (world coords)
  let solid = null, ghosts = null, edges = null, removedLines = null, userGroup = new THREE.Group();
  let ranges = [];                 // [startVertex, building] for picking in the merged mesh
  scene.add(userGroup);

  function ringShape(ring) {
    const s = new THREE.Shape();
    ring.forEach(([x, y], i) => (i ? s.lineTo(x - H, -(y - H)) : s.moveTo(x - H, -(y - H))));
    return s;
  }
  function extrude(b) {
    if (geoCache.has(b.id)) return geoCache.get(b.id);
    let g = null;
    try {
      const shapes = b.rings.filter((r) => r.length >= 4).map(ringShape);
      if (shapes.length) {
        g = new THREE.ExtrudeGeometry(shapes, { depth: Math.max(0.6, b.h), bevelEnabled: false, curveSegments: 1 });
        g.rotateX(-Math.PI / 2);
        g.deleteAttribute('uv');
      }
    } catch (e) { g = null; }
    geoCache.set(b.id, g);
    return g;
  }
  // outline edges: top and bottom loops plus verticals where the footprint turns
  function pushEdges(arr, ring, h) {
    const n = ring.length - 1;   // rings are closed (last == first)
    for (let i = 0; i < n; i++) {
      const [x0, y0] = ring[i], [x1, y1] = ring[i + 1];
      arr.push(x0 - H, h, y0 - H, x1 - H, h, y1 - H, x0 - H, 0.05, y0 - H, x1 - H, 0.05, y1 - H);
      const [xp, yp] = ring[(i - 1 + n) % n];
      const a1 = Math.atan2(y0 - yp, x0 - xp), a2 = Math.atan2(y1 - y0, x1 - x0);
      let d = Math.abs(a2 - a1); if (d > Math.PI) d = 2 * Math.PI - d;
      if (d > 0.35) arr.push(x0 - H, 0.05, y0 - H, x0 - H, h, y0 - H);
    }
  }
  function disposeObj(o) { if (o) { scene.remove(o); o.geometry.dispose(); } }

  function rebuild(st, newSite) {
    if (newSite) { for (const g of geoCache.values()) if (g) g.dispose(); geoCache.clear(); }
    disposeObj(solid); disposeObj(ghosts); disposeObj(edges); disposeObj(removedLines);
    solid = ghosts = edges = removedLines = null;
    const sg = [], gg = [], ev = [], rv = [];
    ranges = []; let v = 0;
    for (const b of st.buildings) {
      if (st.removed.has(b.id)) { for (const r of b.rings) pushEdges(rv, r, Math.max(0.6, b.h)); continue; }
      const g = extrude(b); if (!g) continue;
      if (st.blocks(b)) { sg.push(g); ranges.push([v, b]); v += g.attributes.position.count; for (const r of b.rings) pushEdges(ev, r, b.h); }
      else gg.push(g);
    }
    if (sg.length) { solid = new THREE.Mesh(mergeGeometries(sg, false), bldMat); solid.castShadow = solid.receiveShadow = true; scene.add(solid); }
    if (gg.length) { ghosts = new THREE.Mesh(mergeGeometries(gg, false), ghostMat); scene.add(ghosts); }
    const lg = (arr, mat) => { const g = new THREE.BufferGeometry(); g.setAttribute('position', new THREE.Float32BufferAttribute(arr, 3)); const l = new THREE.LineSegments(g, mat); scene.add(l); return l; };
    if (ev.length) edges = lg(ev, edgeMat);
    if (rv.length) { removedLines = lg(rv, removedMat); removedLines.computeLineDistances(); }
    renderer.shadowMap.needsUpdate = true; dirty = true;
  }

  function rebuildUser(st) {
    for (const c of userGroup.children.slice()) { userGroup.remove(c); c.geometry.dispose(); }
    for (const u of st.user) {
      const g = new THREE.BoxGeometry(u.w, u.h, u.d); g.translate(0, u.h / 2, 0);
      const m = new THREE.Mesh(g, u === st.sel ? selMat : userMat);
      m.position.set(u.cx - H, 0, u.cy - H); m.rotation.y = -u.rot * Math.PI / 180;
      m.castShadow = m.receiveShadow = true; m.userData.u = u;
      const e = new THREE.LineSegments(new THREE.EdgesGeometry(g), edgeMat); m.add(e);
      userGroup.add(m);
    }
    renderer.shadowMap.needsUpdate = true; dirty = true;
  }

  // wind arrow at the upwind edge
  const arrow = new THREE.ArrowHelper(new THREE.Vector3(1, 0, 0), new THREE.Vector3(), 110, 0x2a2926, 34, 20);
  arrow.line.material.linewidth = 2; scene.add(arrow);
  const label = makeLabel('Wind 風');
  scene.add(label);
  function setWind(dirIdx) {
    const b = dirIdx * 22.5 * Math.PI / 180;
    const down = new THREE.Vector3(-Math.sin(b), 0, Math.cos(b));      // towards where the wind goes (x east, z south)
    const start = down.clone().multiplyScalar(-(H + 150)); start.y = 40;
    arrow.position.copy(start); arrow.setDirection(down);
    label.position.copy(start).add(new THREE.Vector3(0, 30, 0));
    dirty = true;
  }
  function makeLabel(text) {
    const c = document.createElement('canvas'); c.width = 256; c.height = 64;
    const x = c.getContext('2d'); x.font = '600 34px ui-sans-serif,system-ui,"Microsoft JhengHei",sans-serif';
    x.fillStyle = '#2a2926'; x.textAlign = 'center'; x.textBaseline = 'middle'; x.fillText(text, 128, 34);
    const t = new THREE.CanvasTexture(c); t.colorSpace = THREE.SRGBColorSpace;
    const s = new THREE.Sprite(new THREE.SpriteMaterial({ map: t, depthTest: false, transparent: true }));
    s.scale.set(80, 20, 1); s.renderOrder = 5; return s;
  }

  // streaks: each particle keeps its last TAIL positions
  let hist = null, nPart = 0, streaks = null;
  const streakMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, depthWrite: false, toneMapped: false });
  function ensureStreaks(n) {
    if (n === nPart && streaks) return;
    if (streaks) { scene.remove(streaks); streaks.geometry.dispose(); }
    nPart = n; hist = new Float32Array(n * TAIL * 2);
    const segs = n * (TAIL - 1);
    const g = new THREE.BufferGeometry();
    g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(segs * 2 * 3), 3).setUsage(THREE.DynamicDrawUsage));
    const col = new Float32Array(segs * 2 * 4);
    for (let i = 0; i < n; i++) for (let k = 0; k < TAIL - 1; k++) {
      const o = ((i * (TAIL - 1) + k) * 2) * 4;
      const a0 = 0.95 * (1 - k / (TAIL - 1)), a1 = 0.95 * (1 - (k + 1) / (TAIL - 1));
      col.set([1, 1, 1, a0, 1, 1, 1, a1], o);
    }
    g.setAttribute('color', new THREE.BufferAttribute(col, 4));
    streaks = new THREE.LineSegments(g, streakMat); streaks.frustumCulled = false; streaks.renderOrder = 3;
    scene.add(streaks);
  }
  function pushParticles(parts, show) {
    ensureStreaks(parts.length);
    streaks.visible = show;
    if (!show) return;
    for (let i = 0; i < parts.length; i++) {
      const p = parts[i], o = i * TAIL * 2;
      if (p.fresh) { for (let k = 0; k < TAIL; k++) { hist[o + k * 2] = p.x; hist[o + k * 2 + 1] = p.y; } p.fresh = false; continue; }
      hist.copyWithin(o + 2, o, o + (TAIL - 1) * 2);
      hist[o] = p.x; hist[o + 1] = p.y;
    }
    const pos = streaks.geometry.attributes.position.array;
    let w = 0;
    for (let i = 0; i < parts.length; i++) {
      const o = i * TAIL * 2;
      for (let k = 0; k < TAIL - 1; k++) {
        pos[w++] = hist[o + k * 2] - H; pos[w++] = STREAK_Y; pos[w++] = hist[o + k * 2 + 1] - H;
        pos[w++] = hist[o + k * 2 + 2] - H; pos[w++] = STREAK_Y; pos[w++] = hist[o + k * 2 + 3] - H;
      }
    }
    streaks.geometry.attributes.position.needsUpdate = true;
    dirty = true;
  }

  function setTheme(dark) {
    const bg = dark ? 0x22221f : 0xe9ecee;
    scene.background = new THREE.Color(bg);
    scene.fog = new THREE.Fog(bg, 1400, 4200);
    bldMat.color.set(dark ? 0xd9d4ca : 0xf4f2ed);
    plinthMat.color.set(dark ? 0x4a4843 : 0xdedad2);
    groundMat.color.set(dark ? 0x3d3b37 : 0xcfcac0);
    dirty = true;
  }

  // clicks (not drags) pick the ground or a block
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
    if (ui) return onPick({ kind: 'user', u: ui.object.userData.u, ev: e });
    const hit = solid ? ray.intersectObject(solid, false)[0] : null;
    const gp = new THREE.Vector3();
    const onGround = ray.ray.intersectPlane(new THREE.Plane(new THREE.Vector3(0, 1, 0), 0), gp);
    if (hit && (!onGround || hit.distance <= ray.ray.origin.distanceTo(gp) + 0.01)) {
      const vi = hit.face.a; let lo = 0, hi = ranges.length - 1;
      while (lo < hi) { const m = (lo + hi + 1) >> 1; if (ranges[m][0] <= vi) lo = m; else hi = m - 1; }
      return onPick({ kind: 'building', b: ranges[lo][1], ev: e });
    }
    // removed buildings have no solid: let a ground click inside one bring it back
    if (onGround && Math.abs(gp.x) <= H && Math.abs(gp.z) <= H) return onPick({ kind: 'ground', x: gp.x + H, y: gp.z + H, ev: e });
  });

  function resize(w, h) {
    renderer.setSize(w, h, false);
    renderer.domElement.style.width = w + 'px'; renderer.domElement.style.height = h + 'px';
    camera.aspect = w / h; camera.updateProjectionMatrix(); dirty = true;
  }
  let fieldDirty = false;
  return {
    THREE,
    rebuild, rebuildUser, setWind, pushParticles, setTheme, resize,
    markField() { fieldDirty = true; dirty = true; },
    showField(on) { field.visible = on; dirty = true; },
    // call every frame; renders only when something changed
    frame() {
      controls.update();
      if (fieldDirty) { fieldTex.needsUpdate = true; fieldDirty = false; }
      if (!dirty) return false;
      renderer.render(scene, camera); dirty = false; return true;
    },
    azimuth() { return controls.getAzimuthalAngle(); },
    camera, controls, renderer,
    stats() { const i = renderer.info; return { calls: i.render.calls, triangles: i.render.triangles, solid: !!solid, user: userGroup.children.length, cam: camera.position.toArray().map((v) => +v.toFixed(1)) }; },
    resetView() { camera.position.set(...HOME); controls.target.set(0, 0, 30); controls.update(); dirty = true; },
    dispose() { renderer.dispose(); renderer.forceContextLoss(); renderer.domElement.remove(); },
  };
}
