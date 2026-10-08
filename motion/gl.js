/* ChronoCell-5D film: the WebGL layer (three.js r128 + bloom).
 * Everything is a pure function of the film time t (seconds): GL.render(t) can be called for any t in any order,
 * which is what frame-by-frame rendering needs. The chromosome is the app's own reference model of chr22 and the
 * deletion movie is the 4D workspace's 22q11.2 preset (motion/build_data.py); the helix and the dust are artwork. */
"use strict";

const GL = (() => {
  const D = window.DATA;
  const W = 1920, H = 1080;
  const canvas = document.getElementById("gl");
  const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, preserveDrawingBuffer: true, powerPreference: "high-performance" });
  renderer.setPixelRatio(1);
  renderer.setSize(W, H, false);
  renderer.setClearColor(0x04050b, 1);
  const scene = new THREE.Scene();
  const camera = new THREE.PerspectiveCamera(34, W / H, 0.01, 200);

  const composer = new THREE.EffectComposer(renderer);
  composer.addPass(new THREE.RenderPass(scene, camera));
  const bloom = new THREE.UnrealBloomPass(new THREE.Vector2(W, H), 0.9, 0.55, 0.12);
  composer.addPass(bloom);

  // ---------------------------------------------------------------------------------------------- helpers
  const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
  const lin = (t, a, b) => clamp((t - a) / (b - a));
  const ease = {
    io3: (x) => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2),
    o3: (x) => 1 - Math.pow(1 - x, 3),
    oExpo: (x) => (x >= 1 ? 1 : 1 - Math.pow(2, -10 * x)),
    ioExpo: (x) => (x <= 0 ? 0 : x >= 1 ? 1 : x < 0.5 ? Math.pow(2, 20 * x - 10) / 2 : (2 - Math.pow(2, -20 * x + 10)) / 2),
    io2: (x) => (x < 0.5 ? 2 * x * x : 1 - Math.pow(-2 * x + 2, 2) / 2),
  };
  const mix = (a, b, p) => a + (b - a) * p;
  function rng(seed) {
    let s = seed >>> 0;
    return () => {
      s = (s + 0x6d2b79f5) >>> 0;
      let q = s;
      q = Math.imul(q ^ (q >>> 15), q | 1);
      q ^= q + Math.imul(q ^ (q >>> 7), q | 61);
      return ((q ^ (q >>> 14)) >>> 0) / 4294967296;
    };
  }
  const hex = (h) => new THREE.Color(h);
  // the app's genomic-position ramp (blue -> violet -> terracotta -> ochre), brightened for a dark background
  const RAMP = [hex("#4553ff"), hex("#8a5cff"), hex("#e0573a"), hex("#ffb347")];
  function ramp(u, out) {
    const x = clamp(u) * (RAMP.length - 1), i = Math.min(Math.floor(x), RAMP.length - 2), f = x - i;
    out.copy(RAMP[i]).lerp(RAMP[i + 1], f);
    return out;
  }

  // ---------------------------------------------------------------------------------------------- point shaders
  const VERT = `
    attribute float asize; attribute float aalpha; attribute vec3 acolor;
    varying vec3 vColor; varying float vAlpha; uniform float uScale;
    void main() {
      vColor = acolor; vAlpha = aalpha;
      vec4 mv = modelViewMatrix * vec4(position, 1.0);
      gl_PointSize = asize * uScale / max(-mv.z, 0.05);
      gl_Position = projectionMatrix * mv;
    }`;
  const FRAG_GLOW = `
    varying vec3 vColor; varying float vAlpha; uniform float uOpacity;
    void main() {
      vec2 d = gl_PointCoord * 2.0 - 1.0; float r = dot(d, d);
      if (r > 1.0) discard;
      float a = pow(1.0 - r, 2.2) * vAlpha * uOpacity;
      gl_FragColor = vec4(vColor * a, a);
    }`;
  const FRAG_BEAD = `   // shaded sphere impostor
    varying vec3 vColor; varying float vAlpha; uniform float uOpacity;
    void main() {
      vec2 d = gl_PointCoord * 2.0 - 1.0; d.y = -d.y; float r = dot(d, d);
      if (r > 1.0) discard;
      vec3 n = vec3(d, sqrt(1.0 - r));
      vec3 L = normalize(vec3(-0.45, 0.55, 0.75));
      float diff = max(dot(n, L), 0.0);
      float spec = pow(max(dot(reflect(-L, n), vec3(0.0, 0.0, 1.0)), 0.0), 24.0);
      vec3 c = vColor * (0.28 + 0.85 * diff) + vec3(0.9) * spec * 0.55 + vColor * pow(1.0 - n.z, 3.0) * 0.6;
      if (vAlpha * uOpacity < 0.02) discard;
      gl_FragColor = vec4(c, vAlpha * uOpacity);
    }`;
  function pointsLayer(n, frag, additive) {
    const g = new THREE.BufferGeometry();
    const pos = new Float32Array(n * 3), col = new Float32Array(n * 3), size = new Float32Array(n), alpha = new Float32Array(n);
    g.setAttribute("position", new THREE.BufferAttribute(pos, 3));
    g.setAttribute("acolor", new THREE.BufferAttribute(col, 3));
    g.setAttribute("asize", new THREE.BufferAttribute(size, 1));
    g.setAttribute("aalpha", new THREE.BufferAttribute(alpha, 1));
    const m = new THREE.ShaderMaterial({
      vertexShader: VERT, fragmentShader: frag, transparent: true, depthTest: !additive, depthWrite: !additive,
      blending: additive ? THREE.AdditiveBlending : THREE.NormalBlending,
      uniforms: { uScale: { value: H * 0.5 }, uOpacity: { value: 1 } },
    });
    const p = new THREE.Points(g, m);
    p.frustumCulled = false;
    return { obj: p, pos, col, size, alpha, mat: m, geo: g,
      dirty() { for (const k of ["position", "acolor", "asize", "aalpha"]) g.attributes[k].needsUpdate = true; } };
  }

  // ---------------------------------------------------------------------------------------------- dust
  const ND = 2600;
  const dust = pointsLayer(ND, FRAG_GLOW, true);
  const dustHome = new Float32Array(ND * 3), dustPh = new Float32Array(ND);
  {
    const r = rng(11);
    const tint = [hex("#6f7cff"), hex("#a98bff"), hex("#ff8a5c"), hex("#9fe7ff")];
    for (let i = 0; i < ND; i++) {
      const u = r() * 2 - 1, th = r() * Math.PI * 2, rad = 2.2 + Math.pow(r(), 0.6) * 9;
      const s = Math.sqrt(1 - u * u);
      dustHome[3 * i] = rad * s * Math.cos(th) * 1.6;
      dustHome[3 * i + 1] = rad * u * 0.9;
      dustHome[3 * i + 2] = rad * s * Math.sin(th) - 2;
      dustPh[i] = r() * 100;
      const c = tint[Math.floor(r() * tint.length)];
      dust.col.set([c.r, c.g, c.b], 3 * i);
      dust.size[i] = 0.06 + r() * 0.12;
    }
  }
  scene.add(dust.obj);

  // ---------------------------------------------------------------------------------------------- double helix
  const NHS = 520, RUNG_EVERY = 4, RUNG_PTS = 5;
  const NR = Math.floor(NHS / RUNG_EVERY) * RUNG_PTS;
  const NH = 2 * NHS + NR;
  const helix = pointsLayer(NH, FRAG_GLOW, true);
  const helixPos = new Float32Array(NH * 3);   // in helix space (x along the axis), before spin
  const helixStart = new Float32Array(NH * 3); // where each particle comes from (inside the dust cloud)
  const helixDelay = new Float32Array(NH), helixAxial = new Float32Array(NH);
  {
    const r = rng(21), L = 5.6, R = 0.46;
    let k = 0;
    const put = (x, y, z, c, s) => {
      helixPos.set([x, y, z], 3 * k);
      helix.col.set([c.r, c.g, c.b], 3 * k);
      helix.size[k] = s;
      const u = r() * 2 - 1, th = r() * Math.PI * 2, rad = 3 + r() * 6, q = Math.sqrt(1 - u * u);
      helixStart.set([rad * q * Math.cos(th) * 1.5, rad * u, rad * q * Math.sin(th) - 1], 3 * k);
      helixDelay[k] = r();
      helixAxial[k] = (x + L / 2) / L;
      k++;
    };
    const cA = hex("#6b7bff"), cB = hex("#c08bff"), cR = hex("#ff9a6a");
    for (let s = 0; s < 2; s++)
      for (let i = 0; i < NHS; i++) {
        const x = -L / 2 + (L * i) / (NHS - 1), a = i * 0.21 + s * Math.PI;
        put(x, R * Math.cos(a), R * Math.sin(a), s ? cB : cA, 0.16);
      }
    for (let i = 0; i < NHS; i += RUNG_EVERY)
      for (let j = 1; j <= RUNG_PTS; j++) {
        const x = -L / 2 + (L * i) / (NHS - 1), a = i * 0.21, f = j / (RUNG_PTS + 1);
        const y = mix(R * Math.cos(a), R * Math.cos(a + Math.PI), f), z = mix(R * Math.sin(a), R * Math.sin(a + Math.PI), f);
        put(x, y, z, cR, 0.09);
      }
  }
  scene.add(helix.obj);

  // ---------------------------------------------------------------------------------------------- chromosome (chr22)
  const F = D.fold;
  const NB = F.xyz.length;
  const fold = new THREE.Group();
  scene.add(fold);
  const beads = pointsLayer(NB, FRAG_BEAD, false);
  const glow = pointsLayer(NB, FRAG_GLOW, true);
  const lineGeo = new THREE.BufferGeometry();
  const linePos = new Float32Array(NB * 3), lineCol = new Float32Array(NB * 3);
  lineGeo.setAttribute("position", new THREE.BufferAttribute(linePos, 3));
  lineGeo.setAttribute("color", new THREE.BufferAttribute(lineCol, 3));
  const lineMat = new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.55, blending: THREE.AdditiveBlending, depthWrite: false });
  const line = new THREE.Line(lineGeo, lineMat);
  line.frustumCulled = false;
  fold.add(line, beads.obj, glow.obj);
  const foldXYZ = new Float32Array(NB * 3), fiberXYZ = new Float32Array(NB * 3);
  {
    const c = new THREE.Color();
    for (let i = 0; i < NB; i++) {
      const p = F.xyz[i];
      foldXYZ.set([p[0], p[1], p[2]], 3 * i);
      const u = F.pos[i];
      // the stretched fibre: the chain laid out along x as a slow coil (artwork; only the end state is the model)
      const x = -3.3 + 6.6 * u, a = u * 38;
      fiberXYZ.set([x, 0.16 * Math.cos(a) + 0.05 * Math.sin(u * 7), 0.16 * Math.sin(a)], 3 * i);
      ramp(u, c);
      if (!F.valid[i]) c.setRGB(0.35, 0.36, 0.42);
      beads.col.set([c.r, c.g, c.b], 3 * i);
      glow.col.set([c.r * 0.9, c.g * 0.9, c.b * 0.9], 3 * i);
      lineCol.set([c.r * 0.8, c.g * 0.8, c.b * 0.8], 3 * i);
      beads.size[i] = F.valid[i] ? 0.06 : 0.035;
      glow.size[i] = 0.2;
    }
  }
  // helix particle -> bead it hands over to (ordered along the axis, so the stream keeps its order)
  const handTo = new Int32Array(NH);
  for (let k = 0; k < NH; k++) handTo[k] = Math.round(helixAxial[k] * (NB - 1));

  // ---------------------------------------------------------------------------------------------- 22q11.2 deletion (4D)
  const S = D.sv;
  const NS = S.disp.length;
  const sv = new THREE.Group();
  scene.add(sv);
  const svBeads = pointsLayer(NS, FRAG_BEAD, false), svGlow = pointsLayer(NS, FRAG_GLOW, true);
  const svLineGeo = new THREE.BufferGeometry();
  const svLinePos = new Float32Array(NS * 3), svLineCol = new Float32Array(NS * 3);
  svLineGeo.setAttribute("position", new THREE.BufferAttribute(svLinePos, 3));
  svLineGeo.setAttribute("color", new THREE.BufferAttribute(svLineCol, 3));
  const svLine = new THREE.Line(svLineGeo, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.5, blending: THREE.AdditiveBlending, depthWrite: false }));
  svLine.frustumCulled = false;
  sv.add(svLine, svBeads.obj, svGlow.obj);
  const calm = hex("#7d89ff"), hot = hex("#ff5a2a"), white = hex("#fff1d6");

  // ---------------------------------------------------------------------------------------------- per-frame state
  const tmp = new THREE.Color();
  function setDust(t, level) {
    for (let i = 0; i < ND; i++) {
      const ph = dustPh[i];
      dust.pos[3 * i] = dustHome[3 * i] + 0.25 * Math.sin(t * 0.21 + ph);
      dust.pos[3 * i + 1] = dustHome[3 * i + 1] + 0.18 * Math.cos(t * 0.17 + ph * 1.3);
      dust.pos[3 * i + 2] = dustHome[3 * i + 2] + 0.2 * Math.sin(t * 0.13 + ph * 0.7);
      dust.alpha[i] = level * (0.35 + 0.65 * (0.5 + 0.5 * Math.sin(t * 1.3 + ph * 3.1)));
    }
    dust.dirty();
  }

  function setHelix(t) {
    // 0.4-4.4 s: particles stream in from the dust and form the helix; 7.6-9.4 s: the helix unwinds into the fibre
    const spin = 0.55 * t, cs = Math.cos(spin), sn = Math.sin(spin);
    const tilt = -0.18, ct = Math.cos(tilt), st = Math.sin(tilt);
    const unwind = ease.io3(lin(t, 7.6, 9.4));
    const fade = 1 - lin(t, 9.0, 9.7);
    let vis = 0;
    for (let k = 0; k < NH; k++) {
      const p = ease.oExpo(lin(t, 0.4 + helixDelay[k] * 2.0, 2.4 + helixDelay[k] * 2.0));
      let x = helixPos[3 * k], y = helixPos[3 * k + 1], z = helixPos[3 * k + 2];
      const y2 = y * cs - z * sn, z2 = y * sn + z * cs;          // spin about the axis
      let hx = x * ct - y2 * st, hy = x * st + y2 * ct, hz = z2;  // tilt the whole helix a little
      const sx = helixStart[3 * k], sy = helixStart[3 * k + 1], sz = helixStart[3 * k + 2];
      const sw = (1 - p) * 1.2;                                   // incoming particles swirl
      let px = mix(sx * Math.cos(sw) - sz * Math.sin(sw), hx, p);
      let py = mix(sy, hy, p);
      let pz = mix(sx * Math.sin(sw) + sz * Math.cos(sw), hz, p);
      if (unwind > 0) {
        const b = handTo[k];
        const q = ease.io3(clamp(unwind * 1.25 - helixAxial[k] * 0.25));
        px = mix(px, fiberXYZ[3 * b], q);
        py = mix(py, fiberXYZ[3 * b + 1], q);
        pz = mix(pz, fiberXYZ[3 * b + 2], q);
      }
      helix.pos[3 * k] = px; helix.pos[3 * k + 1] = py; helix.pos[3 * k + 2] = pz;
      helix.alpha[k] = p * fade * (0.75 + 0.25 * Math.sin(t * 3 + k));
      vis += helix.alpha[k];
    }
    helix.obj.visible = vis > 0.5;
    helix.dirty();
  }

  function setFold(t) {
    // 9.0-9.7 s the fibre appears where the helix ends; 9.6-14.8 s each bead folds into the model, a wave along
    // the chromosome; afterwards the model turns slowly. Opacity by scene.
    let op = 0;
    if (t < 21.5) op = lin(t, 8.9, 9.6) * (1 - lin(t, 20.6, 21.4)) * mix(1, 0.42, ease.io3(lin(t, 15.7, 16.6)));
    else op = lin(t, 84.4, 85.4);
    fold.visible = op > 0.001;
    if (!fold.visible) return;
    for (let i = 0; i < NB; i++) {
      const u = F.pos[i];
      const a = 9.7 + 3.4 * u;
      const p = ease.io3(lin(t, a, a + 1.6));
      const bulge = Math.sin(Math.PI * p);
      const fx = fiberXYZ[3 * i], fy = fiberXYZ[3 * i + 1], fz = fiberXYZ[3 * i + 2];
      const gx = foldXYZ[3 * i], gy = foldXYZ[3 * i + 1], gz = foldXYZ[3 * i + 2];
      const x = mix(fx, gx, p), y = mix(fy, gy, p) + bulge * 0.35 * Math.sin(u * 40), z = mix(fz, gz, p) + bulge * 0.6;
      beads.pos[3 * i] = glow.pos[3 * i] = linePos[3 * i] = x;
      beads.pos[3 * i + 1] = glow.pos[3 * i + 1] = linePos[3 * i + 1] = y;
      beads.pos[3 * i + 2] = glow.pos[3 * i + 2] = linePos[3 * i + 2] = z;
      beads.alpha[i] = op;
      glow.alpha[i] = op * (0.14 + 0.3 * F.epi[i]) * (0.6 + 0.4 * Math.sin(t * 2.2 + u * 60));
    }
    lineMat.opacity = 0.8 * op;
    beads.dirty(); glow.dirty();
    lineGeo.attributes.position.needsUpdate = true;
    // whole-model motion
    const settle = ease.io3(lin(t, 9.6, 15.5));
    fold.rotation.y = t < 21.5 ? 0.18 * (t - 9.6) * settle : 0.2 * (t - 84.4) - 0.6;
    fold.rotation.x = t < 21.5 ? -0.25 * settle : -0.25;
    let sc = 1;
    if (t >= 84) sc = mix(0.82, 1.0, ease.o3(lin(t, 84.4, 86.2)));
    fold.position.set(0, 0, 0);
    if (t > 15.6 && t < 21.5) sc = mix(1, 0.8, ease.io3(lin(t, 15.6, 17.0)));   // title: smaller, behind the name
    fold.scale.setScalar(sc);
  }

  function setSV(t) {
    // the 4D workspace's 22q11.2 deletion preset: 24 relaxation frames after the deletion, played during workspace 02
    const t0 = TIMES.tour + 1 * TIMES.ws, t1 = t0 + TIMES.ws;
    const op = lin(t, t0 + 0.15, t0 + 0.7) * (1 - lin(t, t1 - 0.45, t1));
    sv.visible = op > 0.001;
    if (!sv.visible) return;
    const fp = clamp((t - (t0 + 0.6)) / (TIMES.ws - 1.4)) * (S.frames - 1);
    const f0 = Math.floor(fp), f1 = Math.min(f0 + 1, S.frames - 1), w = fp - f0;
    const A = S.xyz[f0], B = S.xyz[f1];
    const heat = fp / (S.frames - 1);
    for (let i = 0; i < NS; i++) {
      const x = mix(A[i][0], B[i][0], w), y = mix(A[i][1], B[i][1], w), z = mix(A[i][2], B[i][2], w);
      svBeads.pos.set([x, y, z], 3 * i); svGlow.pos.set([x, y, z], 3 * i); svLinePos.set([x, y, z], 3 * i);
      const d = S.disp[i] * heat;
      tmp.copy(calm).lerp(hot, clamp(d * 1.6)).lerp(white, clamp(d * 2 - 1.2));
      svBeads.col.set([tmp.r, tmp.g, tmp.b], 3 * i); svGlow.col.set([tmp.r, tmp.g, tmp.b], 3 * i);
      svLineCol.set([tmp.r * 0.7, tmp.g * 0.7, tmp.b * 0.7], 3 * i);
      svBeads.size[i] = 0.055; svGlow.size[i] = 0.16 + 0.35 * d;
      svBeads.alpha[i] = op; svGlow.alpha[i] = op * (0.08 + 0.6 * d);
    }
    svLine.material.opacity = 0.45 * op;
    svBeads.dirty(); svGlow.dirty();
    svLineGeo.attributes.position.needsUpdate = true; svLineGeo.attributes.color.needsUpdate = true;
    sv.rotation.set(-0.2, 0.25 * (t - t0) - 0.4, 0);
    sv.position.set(0.95, -0.05, 0);
    sv.scale.setScalar(0.95);
  }

  function setCamera(t) {
    let r = 7.2, az = 0, el = 0.0, tx = 0, ty = 0;
    if (t < 7.6) { r = mix(7.6, 6.1, ease.io2(lin(t, 0, 7.6))); }
    else if (t < 16) {
      const q = ease.io3(lin(t, 7.6, 15.8));
      r = mix(6.1, 4.6, q); az = mix(0, 0.55, q); el = mix(0, 0.22, q);
      tx = -1.15 * ease.io3(lin(t, 10.2, 11.8)) * (1 - ease.io3(lin(t, 15.0, 16.0)));
    } else if (t < 21.5) { r = 4.6; az = 0.55 + 0.05 * (t - 16); el = 0.22; }
    else if (t < 84) { r = 5.2; az = 0.02 * (t - 21.5); el = 0.05; }
    else { r = mix(6.0, 5.4, ease.o3(lin(t, 84.4, 92))); az = 0.1; el = 0.12; ty = -0.62; }
    // the target is offset sideways in the camera's own frame, so the model moves across the screen, not in depth
    const cx = Math.cos(az), sx = Math.sin(az);
    const ox = tx * cx, oz = -tx * sx;
    camera.position.set(ox + r * sx * Math.cos(el), ty + r * Math.sin(el), oz + r * cx * Math.cos(el));
    camera.lookAt(ox, ty, oz);
  }

  function render(t) {
    // background dust: bright in the opening, an atmosphere elsewhere
    let level = 0.35;
    if (t < 9) level = 0.25 + 0.75 * (1 - lin(t, 2.5, 6)) * lin(t, 0, 0.6);
    else if (t > 84) level = 0.5;
    setDust(t, level);
    setHelix(t);
    setFold(t);
    setSV(t);
    setCamera(t);
    // bloom: punchy in the 3D scenes, softer behind text
    let strength = 0.75;
    if (t < 16) strength = 0.95;
    if (t > 15.8 && t < 17.2) strength = 0.95 + 1.6 * Math.exp(-Math.pow((t - 16.1) / 0.35, 2)); // title flash
    if (t > 21.5 && t < 84) strength = 0.55;
    bloom.strength = strength;
    composer.render();
  }

  return { render };
})();
