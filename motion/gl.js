/* ChronoCell-5D film: the 3D layer (three.js r128).
 * Everything is a pure function of the film time t, so any frame can be drawn in any order (frame-by-frame rendering).
 * The chromosome is the app's own reference model of chr22, drawn as a lit tube like the 3D workspace draws it; the
 * deletion movie is the 4D workspace's 22q11.2 preset (motion/build_data.py). Two renderers, one per 3D subject, each
 * with a transparent canvas that the film moves into the scene showing it. */
"use strict";

const GL = (() => {
  const D = window.DATA;
  const W = 1920, H = 1080, PR = 1.5;                       // drawn at 1.5x and scaled down: clean tube edges
  const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
  const lin = (t, a, b) => clamp((t - a) / (b - a));
  const mix = (a, b, p) => a + (b - a) * p;
  const io3 = (x) => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2);
  const io2 = (x) => (x < 0.5 ? 2 * x * x : 1 - Math.pow(-2 * x + 2, 2) / 2);
  const hex = (h) => new THREE.Color(h);

  // the app's genomic-position ramp (theme.py: cobalt -> violet -> terracotta -> ochre), a little brighter for video
  const RAMP = ["#3A46F0", "#8350D8", "#E9532A", "#EEAA2A"].map(hex);
  const GREY = hex("#8A8FA0");
  function ramp(u, out) {
    const x = clamp(u) * (RAMP.length - 1), i = Math.min(Math.floor(x), RAMP.length - 2), f = x - i;
    return out.copy(RAMP[i]).lerp(RAMP[i + 1], f);
  }

  function makeView(fov) {
    const canvas = document.createElement("canvas");
    canvas.className = "glc";
    const renderer = new THREE.WebGLRenderer({ canvas, antialias: true, alpha: true, preserveDrawingBuffer: true });
    renderer.setPixelRatio(PR);
    renderer.setSize(W, H, false);
    renderer.setClearColor(0x000000, 0);
    Object.assign(canvas.style, { position: "absolute", left: "0", top: "0", width: W + "px", height: H + "px" });
    const scene = new THREE.Scene();
    scene.add(new THREE.HemisphereLight(0xffffff, 0x2b2e45, 0.78));
    const key = new THREE.DirectionalLight(0xffffff, 0.82);
    key.position.set(-3, 4, 5);
    const rim = new THREE.DirectionalLight(0xaab6ff, 0.5);
    rim.position.set(4, -1.5, -5);
    const fill = new THREE.DirectionalLight(0xfff1e0, 0.22);
    fill.position.set(5, 1, 3);
    scene.add(key, rim, fill);
    const camera = new THREE.PerspectiveCamera(fov, W / H, 0.01, 100);
    return { canvas, renderer, scene, camera };
  }
  function attach(v, host) { if (v.canvas.parentNode !== host) host.appendChild(v.canvas); }
  // camera on an orbit around the origin; (sx, sy) = where on screen the origin lands
  function orbit(cam, r, az, el, sx, sy) {
    cam.position.set(r * Math.sin(az) * Math.cos(el), r * Math.sin(el), r * Math.cos(az) * Math.cos(el));
    cam.lookAt(0, 0, 0);
    cam.setViewOffset(W, H, -(sx - W / 2), -(sy - H / 2), W, H);
  }
  function glowTexture() {
    const c = document.createElement("canvas");
    c.width = c.height = 128;
    const g = c.getContext("2d");
    const gr = g.createRadialGradient(64, 64, 0, 64, 64, 64);
    gr.addColorStop(0, "rgba(255,255,255,1)");
    gr.addColorStop(0.25, "rgba(255,236,190,.85)");
    gr.addColorStop(1, "rgba(255,200,120,0)");
    g.fillStyle = gr;
    g.fillRect(0, 0, 128, 128);
    return new THREE.CanvasTexture(c);
  }

  // ---------------------------------------------------------------------------------------------- chr22 fold
  const F = D.fold, NB = F.xyz.length;
  const fv = makeView(30);
  const curve = new THREE.CatmullRomCurve3(F.xyz.map((p) => new THREE.Vector3(p[0], p[1], p[2])), false, "centripetal");
  const SEG = NB * 3, RAD = 7;
  const foldGeo = new THREE.TubeGeometry(curve, SEG, 0.0165, RAD, false);
  {
    const col = new Float32Array(foldGeo.attributes.position.count * 3), c = new THREE.Color();
    for (let i = 0; i <= SEG; i++) {
      const k = Math.min(NB - 1, Math.round((i / SEG) * (NB - 1)));     // equal bonds: arc length ~ bead index
      if (F.valid[k]) ramp(F.pos[k], c); else c.copy(GREY);
      for (let j = 0; j <= RAD; j++) col.set([c.r, c.g, c.b], 3 * (i * (RAD + 1) + j));
    }
    foldGeo.setAttribute("color", new THREE.BufferAttribute(col, 3));
  }
  const tubeMat = () => new THREE.MeshStandardMaterial({ vertexColors: true, roughness: 0.4, metalness: 0.0 });
  const fold = new THREE.Mesh(foldGeo, tubeMat());
  fold.frustumCulled = false;
  const head = new THREE.Sprite(new THREE.SpriteMaterial({ map: glowTexture(), blending: THREE.AdditiveBlending,
    depthTest: false, transparent: true }));
  const group = new THREE.Group();
  group.add(fold, head);
  fv.scene.add(group);
  const QUADS = RAD * 6;

  /* mode "chapter": the fold draws itself (5.0-9.6 s) on the dark stage, then turns; flies into the camera at the end.
     mode "ws1": the finished model on paper, turning; the Turntable button (25.6 s) speeds it up. */
  function drawFold(t, mode, host, opts) {
    attach(fv, host);
    let prog = 1, r = 3.1, az = 0, el = 0.28, sx = 1250, sy = 520;
    if (mode === "free") {                                   // the finished model, camera given by the caller
      ({ r, az, el, sx, sy } = { r: 6.3, az: 0.15 * t, el: 0.25, sx: 960, sy: 520, ...opts });
    } else if (mode === "chapter") {
      prog = io2(lin(t, 5.0, 9.6));
      az = -0.6 + 0.16 * (t - 4.5) + 0.5 * io3(lin(t, 9.4, 14.5));
      r = mix(3.9, 5.1, io3(lin(t, 4.5, 10.5))) * mix(1, 0.3, io3(lin(t, 14.1, 14.9)));
      el = mix(0.12, 0.3, io3(lin(t, 4.5, 11)));
      sx = mix(1180, 1270, io3(lin(t, 4.5, 9.0)));
      sy = 500;
    } else {
      const k = Math.max(0, t - 25.6);
      az = 0.9 + 0.15 * (t - 24.5) + 1.1 * (k - (1 - Math.exp(-3 * k)) / 3);
      r = 6.3;
      el = 0.25;
      sx = 1300;
      sy = 455;
    }
    const n = Math.max(1, Math.floor(prog * SEG));
    foldGeo.setDrawRange(0, n * QUADS);
    head.visible = prog < 0.999 && prog > 0.001;
    if (head.visible) {
      head.position.copy(curve.getPointAt(Math.min(prog, 1)));
      const s = 0.2 + 0.04 * Math.sin(t * 18);
      head.scale.set(s, s, s);
    }
    group.rotation.set(0, 0, 0);
    orbit(fv.camera, r, az, el, sx, sy);
    fv.renderer.render(fv.scene, fv.camera);
  }

  // ---------------------------------------------------------------------------------------------- 22q11.2 deletion
  const S = D.sv, NS = S.disp.length;
  const sv = makeView(30);
  const svMesh = new THREE.Mesh(new THREE.BufferGeometry(), tubeMat());
  svMesh.frustumCulled = false;
  sv.scene.add(svMesh);
  const calm = hex("#5664F5"), hot = hex("#F0552E"), white = hex("#FFE7B0");
  let svKey = "";
  function drawSV(t, t0, dur, host) {
    attach(sv, host);
    const fp = clamp((t - t0) / dur) * (S.frames - 1);
    const key = fp.toFixed(4);
    if (key !== svKey) {
      svKey = key;
      const f0 = Math.floor(fp), f1 = Math.min(f0 + 1, S.frames - 1), w = fp - f0;
      const A = S.xyz[f0], B = S.xyz[f1];
      const pts = [];
      for (let i = 0; i < NS; i++) pts.push(new THREE.Vector3(mix(A[i][0], B[i][0], w), mix(A[i][1], B[i][1], w), mix(A[i][2], B[i][2], w)));
      const seg = NS * 2, rad = 6;
      const g = new THREE.TubeGeometry(new THREE.CatmullRomCurve3(pts, false, "centripetal"), seg, 0.0155, rad, false);
      const col = new Float32Array(g.attributes.position.count * 3), c = new THREE.Color();
      const heat = fp / (S.frames - 1);
      for (let i = 0; i <= seg; i++) {
        const k = Math.min(NS - 1, Math.round((i / seg) * (NS - 1)));
        const d = S.disp[k] * heat;
        c.copy(calm).lerp(hot, clamp(d * 1.7)).lerp(white, clamp(d * 2.2 - 1.3));
        for (let j = 0; j <= rad; j++) col.set([c.r, c.g, c.b], 3 * (i * (rad + 1) + j));
      }
      g.setAttribute("color", new THREE.BufferAttribute(col, 3));
      svMesh.geometry.dispose();
      svMesh.geometry = g;
    }
    orbit(sv.camera, 5.8, 0.5 + 0.22 * (t - t0), 0.22, 1300, 470);
    sv.renderer.render(sv.scene, sv.camera);
  }

  return { drawFold, drawSV };
})();
