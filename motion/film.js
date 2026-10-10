/* ChronoCell-5D reel: the timeline. Every element's state is a pure function of the film time t, so window.seek(t)
 * can draw any frame in any order (motion/render.py calls it once per frame).
 * Every number on screen comes from window.DATA (motion/build_data.py): the app's model, its contact list, the 4D
 * preset, the Scoreboard rows and the result files. App pictures are real screenshots (motion/capture_app.py). */
"use strict";

// ------------------------------------------------------------------------------------------------- timing (s)
const T = { open: 0, fold: 4.5, map: 14.5, title: 20, ws1: 24.5, ws2: 29, ws3: 34, ws4: 38.5, stack: 43, quantum: 47.5,
  bento: 52.5, score: 58, loops: 67, mol: 72.5, mit: 76.5, fail: 80, words: 84, outro: 90, end: 100 };
const W = 1920, H = 1080, FPS = 60;
const D = window.DATA, F = D.fold, TS = D.tests, HI = D.hi;
const C = { cobalt: "#3340D1", cobaltHi: "#5F6BFF", cobaltDk: "#1F2896", terra: "#E8582C", terraDk: "#A93A16",
  ochre: "#F5B931", ochreDk: "#B07F12", violet: "#6E45D6", violetDk: "#432796", paper: "#F2F2EF", paper2: "#E9E9E4",
  ink: "#1C1E1B", ink2: "#474A45", muted: "#62645F", night: "#111318", pass: "#27B26B", passDk: "#16673D",
  fail: "#F0552E", other: "#8D93A8" };
const STATUS = { pass: C.pass, fail: C.fail, other: C.other };

const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
const lin = (t, a, b) => (b === a ? (t >= b ? 1 : 0) : clamp((t - a) / (b - a)));
const mix = (a, b, p) => a + (b - a) * p;
const E = {
  o2: (x) => 1 - (1 - x) * (1 - x),
  o3: (x) => 1 - Math.pow(1 - x, 3),
  o4: (x) => 1 - Math.pow(1 - x, 4),
  i3: (x) => x * x * x,
  io2: (x) => (x < 0.5 ? 2 * x * x : 1 - Math.pow(-2 * x + 2, 2) / 2),
  io3: (x) => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2),
  oExpo: (x) => (x >= 1 ? 1 : 1 - Math.pow(2, -10 * x)),
  iExpo: (x) => (x <= 0 ? 0 : Math.pow(2, 10 * x - 10)),
  oBack: (x) => { const c1 = 1.7, c3 = c1 + 1; return x <= 0 ? 0 : x >= 1 ? 1 : 1 + c3 * Math.pow(x - 1, 3) + c1 * Math.pow(x - 1, 2); },
};
const hexA = (h, a) => { const n = parseInt(h.slice(1), 16); return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`; };
function rng(seed) {
  let s = seed >>> 0;
  return () => { s = (s + 0x6d2b79f5) >>> 0; let q = s; q = Math.imul(q ^ (q >>> 15), q | 1); q ^= q + Math.imul(q ^ (q >>> 7), q | 61); return ((q ^ (q >>> 14)) >>> 0) / 4294967296; };
}
const fmt = (n) => Math.round(n).toLocaleString("en-US");
const SUBS = "₀₁₂₃₄₅₆₇₈₉";
const subU = (m) => m.replace(/\d/g, (d) => SUBS[+d]);
const subH = (m) => m.replace(/(\d+)/g, "<sub>$1</sub>");

// ------------------------------------------------------------------------------------------------- DOM helpers
function mk(tag, cls, parent, html, style) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  if (style) Object.assign(e.style, style);
  if (parent) parent.appendChild(e);
  return e;
}
const img = (name) => (name.includes("/") ? name : `assets/app/${name}.png`);
function chars(el, text) {
  el.innerHTML = "";
  return [...text].map((c) => { const s = mk("span", "ch", el); s.textContent = c === " " ? " " : c; return s; });
}
function setT(el, tr, op) {
  el.style.transform = tr;
  if (op != null) el.style.opacity = clamp(op).toFixed(4);
}
const show = (el, on) => { el.style.display = on ? "block" : "none"; };
function extrude(col, n = 10, dx = 1, dy = 1.3, soft = "rgba(0,0,0,.28)") {
  const a = [];
  for (let i = 1; i <= n; i++) a.push(`${(i * dx).toFixed(1)}px ${(i * dy).toFixed(1)}px 0 ${col}`);
  a.push(`${(n * dx + 8).toFixed(1)}px ${(n * dy + 18).toFixed(1)}px 30px ${soft}`);
  return a.join(",");
}
/* letters / words drop in (overshoot), optional exit */
function animSpans(spans, lt, o) {
  spans.forEach((s, i) => {
    const a = o.in + i * (o.st ?? 0.03);
    const pin = (o.ein ?? E.oBack)(lin(lt, a, a + (o.din ?? 0.5)));
    const b = o.out != null ? o.out + i * (o.sto ?? 0.01) : 1e9;
    const pout = E.i3(lin(lt, b, b + (o.dout ?? 0.3)));
    const y = (1 - pin) * (o.fy ?? 80) - pout * (o.ty ?? 60);
    const sc = mix(o.fs ?? 1, 1, pin);
    setT(s, `translate3d(0,${y.toFixed(2)}px,0) rotate(${((1 - pin) * (o.fr ?? 0)).toFixed(2)}deg) scale(${sc.toFixed(4)})`,
      clamp(pin * 1.6) * (1 - pout));
  });
}
/* pop: an element scales in with overshoot */
function pop(el, t, at, o = {}) {
  const p = E.oBack(lin(t, at, at + (o.d ?? 0.45)));
  setT(el, `translate3d(${(o.x ?? 0) * (1 - p)}px,${(o.y ?? 30) * (1 - p)}px,0) rotate(${((o.r ?? 0) * (1 - p)).toFixed(2)}deg) scale(${mix(o.s ?? 0.6, 1, p).toFixed(4)})`,
    clamp(p * 2) * (o.op ?? 1));
  return p;
}
const rise = (el, t, at, dy = 24, d = 0.4) => setT(el, `translate3d(0,${((1 - E.o3(lin(t, at, at + d))) * dy).toFixed(2)}px,0)`, lin(t, at, at + d * 0.7));

/* odometer: one rolling column per digit; set(v) shows v (higher columns roll only when the one below passes 9) */
function odo(parent, template, style, cls = "dn") {
  const el = mk("div", "odo " + cls, parent, null, style);
  const dec = template.includes(".") ? template.split(".")[1].length : 0;
  const nd = (template.match(/\d/g) || []).length;
  const parts = [];
  let j = 0;
  for (const ch of template) {
    if (/\d/.test(ch)) {
      const col = mk("span", "col", el), strip = mk("span", "strip", col);
      for (let d = 0; d <= 10; d++) mk("span", "", strip, String(d % 10));
      parts.push({ col, strip, place: nd - 1 - j });
      j++;
    } else parts.push({ sep: mk("span", "sep", el, ch), place: nd - j });
  }
  const digits = parts.filter((p) => p.strip).sort((a, b) => a.place - b.place);
  return {
    el,
    set(v) {
      const N = Math.max(0, v) * Math.pow(10, dec);
      let prev = 0;
      for (const p of digits) {
        let q;
        if (p.place === 0) q = N % 10;
        else q = (Math.floor(N / Math.pow(10, p.place)) % 10) + (prev > 9 ? prev - 9 : 0);
        prev = q;
        p.strip.style.transform = `translateY(${(-q).toFixed(4)}em)`;
      }
      for (const p of parts) {
        const lead = p.place > dec;
        const vis = !lead || N >= Math.pow(10, p.place) * 0.999;
        (p.col || p.sep).style.display = vis ? "inline-block" : "none";
      }
    },
  };
}
const logoSvg = (size, c1, c2, sw = 1.8) => `<svg width="${size}" height="${size}" viewBox="0 0 26 26"><path d="M5 20c3-9 5-13 8-13s5 4 8 13" stroke="${c1}" stroke-width="${sw}"/><path d="M5 6c3 9 5 13 8 13s5-4 8-13" stroke="${c2}" stroke-width="${sw}"/></svg>`;
const CHECK = (size, col, sw = 5) => `<svg width="${size}" height="${size}" viewBox="0 0 40 40"><path d="M9 21 L17 29 L32 12" fill="none" stroke="${col}" stroke-width="${sw}" stroke-linecap="round" stroke-linejoin="round"/></svg>`;
const CURSOR_SVG = `<svg width="46" height="58" viewBox="0 0 23 29"><path d="M2 1.5 L2 23 L7.6 18 L11.4 26.6 L15 25 L11.3 16.6 L18.6 16.4 Z" fill="#fff" stroke="#111" stroke-width="1.5" stroke-linejoin="round"/></svg>`;
const CROWN = `<svg width="64" height="48" viewBox="0 0 64 48"><path d="M4 40 L8 10 L22 24 L32 4 L42 24 L56 10 L60 40 Z" fill="${C.ochre}" stroke="${C.ochreDk}" stroke-width="3" stroke-linejoin="round"/></svg>`;

// ------------------------------------------------------------------------------------------------- registries
const CUES = [];
const cue = (t, k, x = {}) => CUES.push({ t: +t.toFixed(3), k, ...x });
const BURSTS = [];
function burst(t, x, y, o = {}) {
  BURSTS.push({ t, x, y, n: 110, seed: 13 + BURSTS.length * 977, speed: 1150, spread: Math.PI * 2, angle: -Math.PI / 2,
    life: 2.3, colors: ["#ffffff", C.ochre, C.terra, C.cobaltHi, C.pass], ...o });
  if (o.sound !== false) cue(t, "sparkle");
}
const SHAKES = [];
const shake = (t, amp = 14) => SHAKES.push({ t, amp });
const FLASHES = [];
const flash = (t, a = 0.6, d = 0.12) => FLASHES.push({ t, a, d });
const TRACKS = [];
function track(from, to, keys, clicks = [], ring = "#ffffff") {
  TRACKS.push({ from, to, keys, clicks, ring });
  clicks.forEach((c) => cue(c, "click"));
}
function trackPos(tr, t) {
  let x = tr.keys[0][1], y = tr.keys[0][2];
  for (let i = 1; i < tr.keys.length; i++) {
    const [t0, x0, y0] = tr.keys[i - 1], [t1, x1, y1] = tr.keys[i];
    if (t >= t0) { const p = E.io3(lin(t, t0, t1)); x = mix(x0, x1, p); y = mix(y0, y1, p); }
  }
  return [x, y];
}

// ------------------------------------------------------------------------------------------------- stage
const stage = document.getElementById("stage");
const world = mk("div", "layer", stage);
const L = { sc: mk("div", "layer persp", world) };
L.fx = mk("canvas", "c2d", world);
L.fx.width = W; L.fx.height = H;
const fxg = L.fx.getContext("2d");
L.hud = mk("div", "layer", stage);
L.cur = mk("div", "layer", stage);
L.ov = mk("div", "layer", stage);

const SCENES = [];
function scene(o) {
  const root = mk("div", "scene", L.sc);
  const s = { ...o, root, bg: mk("div", "fill", root), host: mk("div", "fill", root), ui: mk("div", "fill persp", root) };
  o.build(s);
  SCENES.push(s);
  if (o.inT && o.inT.sound !== null) cue(o.inT.at, o.inT.sound || "whoosh");
  return s;
}
function trans(s, t) {
  let tf = "", clip = "none", op = 1;
  const i = s.inT, o = s.outT;
  if (i) {
    const d = i.d || 0.45, p = E.io3(lin(t, i.at, i.at + d));
    if (p < 1) {
      if (i.type === "iris") {
        const x = i.x ?? 960, y = i.y ?? 540, R = Math.hypot(Math.max(x, W - x), Math.max(y, H - y)) + 30;
        clip = `circle(${(E.i3(lin(t, i.at, i.at + d)) * 0.25 * R + E.io3(lin(t, i.at, i.at + d)) * 0.75 * R).toFixed(1)}px at ${x}px ${y}px)`;
      } else if (i.type === "push-up") tf += `translate3d(0,${((1 - p) * H).toFixed(1)}px,0) `;
      else if (i.type === "push-left") tf += `translate3d(${((1 - p) * W).toFixed(1)}px,0,0) `;
      else if (i.type === "zoom") { tf += `scale(${mix(1.3, 1, E.oExpo(lin(t, i.at, i.at + d))).toFixed(4)}) `; op *= lin(t, i.at, i.at + 0.2); }
    }
  }
  if (o) {
    const d = o.d || 0.45, p = E.io3(lin(t, o.at, o.at + d));
    if (p > 0) {
      if (o.type === "push-up") tf += `translate3d(0,${(-p * H).toFixed(1)}px,0) `;
      else if (o.type === "push-left") tf += `translate3d(${(-p * W).toFixed(1)}px,0,0) `;
      else if (o.type === "zoom") { tf += `scale(${mix(1, 2.4, E.iExpo(lin(t, o.at, o.at + d))).toFixed(4)}) `; op *= 1 - lin(t, o.at + 0.12, o.at + d); }
      else if (o.type === "shrink") tf += `scale(${mix(1, 0.88, p).toFixed(4)}) `;
    }
  }
  s.root.style.transform = tf || "none";
  s.root.style.clipPath = clip;
  s.root.style.opacity = op.toFixed(4);
}

// ---- backgrounds and furniture
function bgDark(s, glow, gx = 50, gy = 50) {
  mk("div", "fill", s.bg, null, { background: C.night });
  mk("div", "fill dots-dark", s.bg);
  mk("div", "fill", s.bg, null, { background: `radial-gradient(circle at ${gx}% ${gy}%, ${hexA(glow, 0.55)} 0%, ${hexA(glow, 0.2)} 26%, rgba(0,0,0,0) 60%)` });
}
function bgPaper(s, tone = C.paper) {
  mk("div", "fill", s.bg, null, { background: tone });
  mk("div", "fill dots-light", s.bg);
}
function bgBurst(s, color) {
  mk("div", "fill", s.bg, null, { background: color });
  const b = mk("div", "burst", s.bg);
  mk("div", "fill", s.bg, null, { background: "radial-gradient(circle at 50% 45%, rgba(255,255,255,.24) 0%, rgba(255,255,255,0) 42%, rgba(0,0,0,.24) 100%)" });
  return b;
}
function rings(parent, x, y, radii, color) {
  const size = radii[radii.length - 1] * 2 + 40, c = size / 2;
  return mk("div", "abs", parent, `<svg width="${size}" height="${size}" viewBox="0 0 ${size} ${size}">` +
    radii.map((r, i) => `<circle cx="${c}" cy="${c}" r="${r}" fill="none" stroke="${color}" stroke-width="${i === 1 ? 3 : 2}"${i === 1 ? ' stroke-dasharray="46 30"' : ""}/>`).join("") + "</svg>",
  { left: `${x - c}px`, top: `${y - c}px`, width: `${size}px`, height: `${size}px` });
}
function headline(parent, text, x, y, o = {}) {
  const el = mk("div", "abs " + (o.cls || "dc"), parent, null, { left: `${x}px`, top: `${y}px`, fontSize: `${o.size || 140}px`,
    color: o.color || "#fff", textShadow: extrude(o.ext || "rgba(0,0,0,.4)", o.depth ?? 10, 1, 1.3, o.soft) });
  return { el, ch: chars(el, text) };
}
function kick(parent, text, x, y, color) { return mk("div", "abs kick", parent, text, { left: `${x}px`, top: `${y}px`, color }); }
function lead(parent, text, x, y, w, color, size = 28) {
  return mk("div", "abs lead", parent, text, { left: `${x}px`, top: `${y}px`, width: `${w}px`, color, fontSize: `${size}px` });
}
function card(parent, name, x, y, w, h, tag, o = {}) {
  const c = mk("div", "card", parent, null, { left: `${x}px`, top: `${y}px`, width: `${w}px`, height: `${h}px`, ...o.style });
  const im = mk("img", "", c);
  im.src = img(name);
  if (o.pos) im.style.objectPosition = o.pos;
  if (tag) mk("div", "tag", c, tag, o.tagStyle);
  return c;
}
function win(parent, name, x, y, w, h, url) {
  const e = mk("div", "win", parent, null, { left: `${x}px`, top: `${y}px`, width: `${w}px`, height: `${h}px` });
  const bar = mk("div", "bar", e, "<i></i><i></i><i></i>");
  mk("div", "url", bar, url);
  const vp = mk("div", "vp", e);
  const im = mk("img", "", vp);
  im.src = img(name);
  return e;
}
function stamp(parent, text, x, y, rot, o) {
  const el = mk("div", "stamp", parent, text, { left: `${x}px`, top: `${y}px`, color: o.color, borderColor: o.border || o.color,
    background: o.bg || "transparent", boxShadow: o.shadow || "0 18px 40px rgba(0,0,0,.25)" });
  el.dataset.rot = rot;
  return el;
}
function stampAnim(el, t, at) {
  const p = lin(t, at, at + 0.22), q = E.o3(p);
  setT(el, `rotate(${(+el.dataset.rot + (1 - q) * 8).toFixed(2)}deg) scale(${mix(2.6, 1, q).toFixed(4)})`, clamp(p * 3));
}

// ================================================================================================= 01 the workstation
const OPEN_CLICK = 2.0;
track(1.0, 3.5, [[1.0, 1640, 1010], [1.85, 1332, 455], [3.5, 1380, 640]], [OPEN_CLICK]);
scene({ t0: 0, t1: T.fold + 0.35,
  build(s) {
    bgDark(s, C.cobalt, 30, 50);
    const words = [...F.chromosomes, "Micro-C", "Hi-C", "H3K27ac", "CTCF", "cohesin", "TADs", "loops", "compartments", "enhancers"];
    const wrap = mk("div", "abs", s.bg, null, { left: "-520px", top: "-280px", width: "2960px", height: "1640px", transform: "rotate(-11deg)" });
    s.rows = [];
    for (let r = 0; r < 10; r++) {
      const big = r % 2 === 1;
      const row = mk("div", "abs", wrap, null, { left: "0", top: `${r * 160}px`, whiteSpace: "nowrap",
        font: `800 ${big ? 68 : 44}px/1 var(--display)`, fontStretch: big ? "100%" : "85%",
        color: r % 3 === 1 ? hexA(C.cobaltHi, 0.3) : "rgba(255,255,255,.075)", filter: `blur(${big ? 2.4 : 1.1}px)` });
      const seq = [];
      for (let k = 0; k < 44; k++) seq.push(words[(k * 7 + r * 5) % words.length]);
      row.textContent = seq.join("    ·    ");
      s.rows.push(row);
    }
    s.ring = rings(s.bg, 485, 540, [150, 270, 420], "rgba(130,145,255,.18)");
    s.p3 = mk("div", "abs", s.ui, null, { left: "410px", top: "465px", width: "1100px", height: "150px", transformOrigin: "75px 75px" });
    const pill = mk("div", "abs", s.p3, null, { left: "0", top: "0", width: "1100px", height: "150px", borderRadius: "999px",
      background: `linear-gradient(180deg, #5A67FF 0%, ${C.cobalt} 72%)`,
      boxShadow: `0 34px 90px ${hexA(C.cobalt, 0.6)}, inset 0 2px 0 rgba(255,255,255,.35), 0 0 0 1px rgba(255,255,255,.14)` });
    s.disc = mk("div", "abs logo", pill, logoSvg(70, C.cobalt, C.terra, 2.3), { left: "20px", top: "20px", width: "110px", height: "110px",
      borderRadius: "50%", background: "#fff", display: "flex", alignItems: "center", justifyContent: "center" });
    mk("div", "abs", pill, "ChronoCell-5D", { left: "158px", top: "28px", font: "900 66px/1 var(--display)", fontStretch: "100%", letterSpacing: "-0.01em", color: "#fff" });
    mk("div", "abs", pill, "CHROMATIN 3D / 4D WORKSTATION", { left: "162px", top: "102px", font: "600 19px/1 var(--mono)", letterSpacing: ".18em", color: "rgba(255,255,255,.75)" });
    s.btn = mk("div", "abs", pill, null, { left: "776px", top: "33px", width: "300px", height: "84px", borderRadius: "999px",
      background: "#fff", overflow: "hidden", boxShadow: "0 6px 16px rgba(0,0,0,.18)" });
    const bs = { left: "0", top: "0", width: "300px", height: "84px", font: "900 34px/84px var(--display)", fontStretch: "88%", textAlign: "center" };
    s.btnA = mk("div", "abs", s.btn, "Open chr22", { ...bs, color: C.cobalt });
    s.btnB = mk("div", "abs", s.btn, "chr22 ✓", { ...bs, color: C.ink, background: C.ochre });
    burst(OPEN_CLICK + 0.04, 1332, 452, { n: 46, speed: 760, life: 1.4, colors: ["#fff", C.ochre, C.cobaltHi], sound: false });
    cue(0.2, "rise", { d: 0.9 });
    cue(1.08, "hit");
    cue(3.45, "rise", { d: 0.85 });
  },
  update(s, t) {
    s.rows.forEach((r, i) => { r.style.transform = `translate3d(${((i % 2 ? -1 : 1) * t * (38 + i * 6) - 260).toFixed(1)}px,0,0)`; });
    s.ring.style.transform = `rotate(${(t * 14).toFixed(2)}deg)`;
    const p = E.oBack(lin(t, 0.25, 1.05)), q = clamp(p);
    const settle = E.io3(lin(t, 2.45, 3.3));
    const zp = lin(t, 3.55, 4.4), z = Math.exp(3.6 * E.i3(zp));        // fly through the logo disc
    const rx = mix(26, 8, q) * (1 - settle), rz = mix(-24, -6, q) * (1 - settle);
    const tx = (960 - 485) * E.io3(zp), ty = (1 - q) * 160;
    setT(s.p3, `translate3d(${tx.toFixed(1)}px,${ty.toFixed(1)}px,0) rotateX(${rx.toFixed(2)}deg) rotateZ(${rz.toFixed(2)}deg) scale(${(mix(0.55, 1, q) * mix(1, 1.06, settle) * z).toFixed(4)})`, clamp(p * 2));
    const b = E.oBack(lin(t, OPEN_CLICK, OPEN_CLICK + 0.32));
    s.btnA.style.transform = `translateY(${(-b * 84).toFixed(1)}px)`;
    s.btnB.style.transform = `translateY(${((1 - b) * 84).toFixed(1)}px)`;
    s.disc.firstChild.style.opacity = (1 - lin(zp, 0.1, 0.45)).toFixed(3);
  } });

// ================================================================================================= 02 the fold
scene({ t0: T.fold - 0.2, t1: T.map + 0.46, inT: { type: "iris", at: T.fold - 0.2, d: 0.55 }, outT: { type: "zoom", at: T.map, d: 0.45 },
  build(s) {
    bgDark(s, C.cobalt, 66, 46);
    const u = s.ui;
    s.kick = kick(u, "CHR22 · HUMAN · GRCh38", 118, 196, "#8E98FF");
    s.kick.style.textTransform = "none";
    s.row = mk("div", "abs", u, null, { left: "108px", top: "236px", display: "flex", alignItems: "baseline", gap: "18px" });
    s.mb = odo(s.row, F.size_mb.toFixed(2), { fontSize: "210px", color: "#fff", textShadow: extrude(C.cobaltDk, 12) });
    mk("div", "dn", s.row, "Mb", { fontSize: "72px", color: "#fff" });
    s.lab = kick(u, "OF DNA, DRAWN BEAD BY BEAD", 118, 468, "rgba(255,255,255,.82)");
    const cr = mk("div", "abs", u, null, { left: "118px", top: "520px", display: "flex", gap: "14px" });
    s.chips = [`${fmt(F.beads)} beads`, `${F.resolution_kb} kb each`, `R<sub style="font-size:.7em;letter-spacing:0">g</sub>&nbsp;${F.rg_nm} nm`].map((h) =>
      mk("div", "pill", cr, h, { position: "relative", textTransform: "none", letterSpacing: ".02em", fontSize: "22px",
        background: "rgba(255,255,255,.08)", border: "1px solid rgba(255,255,255,.2)", color: "#fff" }));
    const h1 = mk("div", "abs dc", u, null, { left: "108px", top: "640px", fontSize: "124px" });
    const a = mk("div", "", h1, null, { color: "#fff", textShadow: extrude(C.cobaltDk, 10) });
    const b = mk("div", "", h1, null, { color: C.ochre, textShadow: extrude(C.ochreDk, 10) });
    s.l1 = chars(a, "THE SHAPE");
    s.l2 = chars(b, "IS THE SWITCH");
    s.sent = lead(u, "Beads that touch in 3D can switch each other’s genes on or off.", 118, 878, 740, "rgba(255,255,255,.8)", 25);
    const ide = mk("div", "abs", u, null, { left: "930px", top: "952px", width: "870px", height: "14px", borderRadius: "7px", background: "rgba(255,255,255,.12)", overflow: "hidden" });
    s.ideFill = mk("div", "abs", ide, null, { left: "0", top: "0", width: "870px", height: "14px", background: "linear-gradient(90deg,#3A46F0,#8350D8,#E9532A,#EEAA2A)" });
    s.ideLab = [mk("div", "abs mono", u, "pter · 0 Mb", { left: "930px", top: "978px", font: "500 16px/1 var(--mono)", color: "rgba(255,255,255,.5)" }),
      mk("div", "abs mono", u, `${F.size_mb} Mb · qter`, { left: "1660px", top: "978px", font: "500 16px/1 var(--mono)", color: "rgba(255,255,255,.5)" })];
    s.ideMark = mk("div", "abs", u, "", { top: "916px", font: "600 19px/1 var(--mono)", color: "#fff", whiteSpace: "nowrap" });
    s.note = mk("div", "abs mono", u, "THE APP’S SYNTHETIC REFERENCE MODEL OF CHR22 · SEED 7", { right: "40px", top: "62px", font: "500 15px/1 var(--mono)", letterSpacing: ".12em", color: "rgba(255,255,255,.45)" });
    cue(5.0, "rise", { d: 4.6 });
    cue(9.62, "hit");
    flash(9.62, 0.18, 0.1);
    s.chips.forEach((_, i) => cue(5.7 + i * 0.16, "pop"));
  },
  update(s, t) {
    GL.drawFold(t, "chapter", s.host);
    const prog = E.io2(lin(t, 5.0, 9.6));
    s.mb.set(prog * F.size_mb);
    rise(s.kick, t, 4.75);
    pop(s.row, t, 4.85, { s: 0.7, y: 40 });
    rise(s.lab, t, 5.25);
    s.chips.forEach((c, i) => pop(c, t, 5.7 + i * 0.16, { s: 0.5, y: 20 }));
    animSpans(s.l1, t, { in: 9.9, st: 0.035, din: 0.5, fy: -140, fs: 1.5 });
    animSpans(s.l2, t, { in: 10.25, st: 0.035, din: 0.5, fy: -140, fs: 1.5 });
    rise(s.sent, t, 10.9);
    const io = lin(t, 4.9, 5.3);
    s.ideFill.style.clipPath = `inset(0 ${((1 - prog) * 100).toFixed(2)}% 0 0)`;
    s.ideFill.parentNode.style.opacity = io;
    s.ideLab.forEach((e) => (e.style.opacity = io));
    s.ideMark.textContent = `chr22 : ${(prog * F.size_mb).toFixed(2)} Mb`;
    s.ideMark.style.left = `${(930 + prog * 870 - 80).toFixed(1)}px`;
    s.ideMark.style.opacity = io;
    s.note.style.opacity = (0.9 * lin(t, 5.4, 5.9)).toFixed(3);
  } });

// ================================================================================================= 03 the contacts
const MAP = F.map, NBIN = MAP.bins;
const MX = 960, MY = 118, MS = 860, PX = MS / NBIN;
let HOV = [40, 52];
{
  let best = -1;
  for (let i = 24; i < 100; i++) for (let j = i + 8; j < Math.min(NBIN, i + 21); j++) {
    const c = MAP.count[i * NBIN + j];
    if (c > best) { best = c; HOV = [i, j]; }
  }
}
const HOVX = MX + (HOV[1] + 0.5) * PX, HOVY = MY + (HOV[0] + 0.5) * PX;
track(17.9, 20.2, [[17.9, 1880, 1070], [18.55, HOVX + 2, HOVY + 2], [20.2, HOVX + 20, HOVY + 14]], []);
const CMAP = (() => {
  const st = [[0, [248, 248, 246]], [0.25, [235, 201, 181]], [0.55, [217, 116, 74]], [0.8, [168, 56, 15]], [1, [74, 23, 5]]];
  const out = [];
  for (let k = 0; k < 256; k++) {
    const u = k / 255;
    let i = 0;
    while (i < st.length - 2 && u > st[i + 1][0]) i++;
    const f = (u - st[i][0]) / (st[i + 1][0] - st[i][0]);
    const c = st[i][1].map((v, n) => Math.round(mix(v, st[i + 1][1][n], clamp(f))));
    out.push(`rgb(${c[0]},${c[1]},${c[2]})`);
  }
  return out;
})();
scene({ t0: T.map, t1: T.title + 0.55, inT: { type: "zoom", at: T.map, d: 0.45, sound: null },
  build(s) {
    bgPaper(s);
    const u = s.ui;
    s.frame = mk("div", "abs", u, null, { left: `${MX - 26}px`, top: `${MY - 26}px`, width: `${MS + 52}px`, height: `${MS + 52}px`, borderRadius: "26px", background: "#fff", boxShadow: "0 34px 80px rgba(0,0,0,.16), 0 0 0 1px rgba(0,0,0,.05)" });
    s.cv = mk("canvas", "c2d", u);
    s.cv.width = W; s.cv.height = H;
    s.g = s.cv.getContext("2d");
    s.px = [];
    const r = rng(99), cx = MX + MS / 2, cy = MY + MS / 2;
    for (let i = 0; i < NBIN; i++) for (let j = i; j < NBIN; j++) {
      const lev = MAP.level[i * NBIN + j];
      if (lev <= 2) continue;
      const delay = 15.0 + ((j - i) / (NBIN - 1)) * 1.9 + r() * 0.35;
      const a = r() * Math.PI * 2, rad = 520 + r() * 760;
      s.px.push([i, j, lev, delay, cx + Math.cos(a) * rad, cy + Math.sin(a) * rad * 0.75]);
    }
    s.cache = document.createElement("canvas");
    s.cache.width = W; s.cache.height = H;
    const g = s.cache.getContext("2d");
    for (const [i, j, lev] of s.px) {
      g.fillStyle = CMAP[lev];
      g.fillRect(MX + j * PX, MY + i * PX, PX + 0.6, PX + 0.6);
      g.fillRect(MX + i * PX, MY + j * PX, PX + 0.6, PX + 0.6);
    }
    s.kick = kick(u, "SIMULATED MICRO-C · CHR22", 118, 196, C.terra);
    s.odo = odo(u, fmt(F.contacts), { position: "absolute", left: "108px", top: "236px", fontSize: "176px", color: C.ink, textShadow: extrude(C.terra, 11, 1, 1.3, "rgba(0,0,0,.12)") });
    s.lab = kick(u, "CONTACT PAIRS", 118, 432, C.ink);
    s.lab2 = mk("div", "abs mono", u, `${fmt(MAP.total)} counts · ${NBIN} × ${NBIN} squares of ${MAP.bin_kb} kb`, { left: "118px", top: "472px", font: "500 21px/1 var(--mono)", color: C.muted });
    const h1 = mk("div", "abs dc", u, null, { left: "108px", top: "560px", fontSize: "124px" });
    s.l1 = chars(mk("div", "", h1, null, { color: C.ink, textShadow: extrude("#C9C9C2", 9, 1, 1.3, "rgba(0,0,0,.1)") }), "WHO TOUCHES");
    s.l2 = chars(mk("div", "", h1, null, { color: C.terra, textShadow: extrude(C.terraDk, 9, 1, 1.3, "rgba(0,0,0,.12)") }), "WHOM");
    s.sent = lead(u, "Experiments such as Hi-C and Micro-C measure which pieces of DNA touch. Darker squares touch more often.", 118, 800, 740, C.ink2, 25);
    s.tip = mk("div", "abs", u, null, { left: `${HOVX > 1480 ? HOVX - 330 : HOVX + 34}px`, top: `${HOVY - 92}px`, padding: "14px 18px", borderRadius: "14px", background: C.ink, color: "#fff", boxShadow: "0 16px 34px rgba(0,0,0,.3)", whiteSpace: "nowrap" });
    const mb = (k) => ((k * MAP.bin_kb) / 1000).toFixed(1);
    mk("div", "mono", s.tip, `${mb(HOV[0])} Mb × ${mb(HOV[1])} Mb`, { font: "500 18px/1 var(--mono)", color: "rgba(255,255,255,.7)" });
    mk("div", "dn", s.tip, `${fmt(MAP.count[HOV[0] * NBIN + HOV[1]])} counts`, { fontSize: "34px", marginTop: "8px" });
    flash(T.map + 0.02, 0.55, 0.1);
    cue(T.map, "whoosh");
    cue(15.0, "roll", { d: 2.6 });
    cue(18.6, "pop");
  },
  update(s, t) {
    pop(s.frame, t, 14.6, { s: 0.9, y: 30 });
    const g = s.g;
    g.clearRect(0, 0, W, H);
    if (t > 17.95) g.drawImage(s.cache, 0, 0);
    else {
      for (const [i, j, lev, delay, sx, sy] of s.px) {
        const dt = t - delay;
        if (dt < 0) continue;
        const q = E.o3(clamp(dt / 0.6)), sz = PX * mix(0.35, 1, q) + 0.6;
        g.fillStyle = CMAP[lev];
        const cx = MX + MS / 2, cy = MY + MS / 2;
        g.fillRect(mix(sx, MX + j * PX, q), mix(sy, MY + i * PX, q), sz, sz);
        g.fillRect(mix(cx + (sy - cy), MX + i * PX, q), mix(cy + (sx - cx), MY + j * PX, q), sz, sz);
      }
    }
    const hp = lin(t, 18.55, 18.8);
    if (hp > 0) {
      g.strokeStyle = C.cobalt;
      g.lineWidth = 3;
      const k = 3 + 2 * Math.sin((t - 18.55) * 9);
      g.globalAlpha = hp;
      g.strokeRect(HOVX - PX / 2 - k, HOVY - PX / 2 - k, PX + 2 * k, PX + 2 * k);
      g.globalAlpha = 1;
    }
    pop(s.tip, t, 18.62, { s: 0.6, y: 16 });
    s.odo.set(F.contacts * E.o3(lin(t, 15.0, 17.6)));
    rise(s.kick, t, 14.75);
    pop(s.odo.el, t, 14.85, { s: 0.7, y: 40 });
    rise(s.lab, t, 15.3);
    rise(s.lab2, t, 15.45);
    animSpans(s.l1, t, { in: 16.2, st: 0.035, din: 0.5, fy: -140, fs: 1.5 });
    animSpans(s.l2, t, { in: 16.55, st: 0.05, din: 0.5, fy: -140, fs: 1.5 });
    rise(s.sent, t, 17.1);
  } });

// ================================================================================================= title
scene({ t0: T.title - 0.2, t1: T.ws1 + 0.46, inT: { type: "iris", at: T.title - 0.2, d: 0.5, x: HOVX, y: HOVY, sound: null }, outT: { type: "push-up", at: T.ws1 },
  build(s) {
    s.burst = bgBurst(s, C.cobalt);
    const u = s.ui;
    s.logo = mk("div", "abs logo", u, logoSvg(190, "#fff", C.ochre, 2.2), { left: "865px", top: "130px", width: "190px", height: "190px" });
    s.name = chars(mk("div", "abs dn", u, null, { left: "0", width: "1920px", top: "400px", textAlign: "center", fontSize: "196px", color: "#fff", textShadow: extrude(C.cobaltDk, 13) }), "ChronoCell-5D");
    s.bar = mk("div", "abs", u, null, { left: "560px", top: "640px", width: "800px", height: "14px", borderRadius: "7px", background: C.ochre });
    s.tag = mk("div", "pill", u, "Chromatin 3D / 4D workstation", { left: "0", top: "700px", background: C.ink, color: "#fff", fontSize: "30px", padding: "20px 34px" });
    s.small = mk("div", "abs mono", u, `8 WORKSPACES · ${TS.n} TESTS ON HELD-OUT REAL DATA`, { left: "0", width: "1920px", textAlign: "center", top: "830px", font: "600 22px/1 var(--mono)", letterSpacing: ".2em", color: "rgba(255,255,255,.82)" });
    cue(T.title + 0.25, "boom");
    flash(T.title + 0.25, 0.55, 0.14);
    shake(T.title + 0.3, 16);
    burst(T.title + 0.42, 960, 520, { n: 170, speed: 1500, life: 2.6, colors: ["#fff", C.ochre, C.terra, "#9AA4FF", C.ink] });
  },
  update(s, t) {
    s.burst.style.transform = `rotate(${(t * 7).toFixed(2)}deg)`;
    pop(s.logo, t, T.title + 0.2, { s: 0.2, r: -40, y: 0, d: 0.5 });
    animSpans(s.name, t, { in: T.title + 0.3, st: 0.035, din: 0.55, fy: -320, fr: 10 });
    const bp = E.o4(lin(t, T.title + 0.95, T.title + 1.4));
    setT(s.bar, `scaleX(${bp.toFixed(4)})`, 1);
    if (!s.tag.dataset.w) s.tag.dataset.w = 1;
    s.tag.style.left = `${(960 - s.tag.offsetWidth / 2).toFixed(1)}px`;
    pop(s.tag, t, T.title + 1.25, { s: 0.6, y: 30 });
    rise(s.small, t, T.title + 1.6);
  } });

// ================================================================================================= 04 3D structure (ws 01)
const TURN = 25.6;
track(25.0, 27.5, [[25.0, 1560, 1070], [25.5, 1040, 889], [27.5, 1180, 1010]], [TURN], C.cobalt);
scene({ t0: T.ws1, t1: T.ws2 + 0.46, inT: { type: "push-up", at: T.ws1 }, outT: { type: "shrink", at: T.ws2 },
  build(s) {
    bgPaper(s);
    const u = s.ui;
    s.kick = kick(u, "WORKSPACE 01", 112, 150, C.cobalt);
    s.h = headline(u, "3D STRUCTURE", 102, 188, { size: 150, color: C.ink, ext: C.cobalt, depth: 10, soft: "rgba(0,0,0,.12)" });
    s.lead = lead(u, "See how one chromosome is folded inside the nucleus, and measure it.", 112, 345, 640, C.ink2);
    const M = [["Radius of gyration", F.rg_nm.toFixed(1), "nm"], ["Scaling exponent ν", F.nu.toFixed(3), ""],
      ["Max 3D span", fmt(F.span_nm), "nm"], ["Contour length", F.contour_um.toFixed(2), "µm"]];
    s.m = M.map(([l, v, un], k) => {
      const e = mk("div", "metric", u, null, { left: `${112 + (k % 2) * 362}px`, top: `${480 + Math.floor(k / 2) * 190}px` });
      mk("div", "l", e, l);
      const row = mk("div", "v", e);
      const o = odo(row, v, null, "");
      if (un) mk("span", "u", row, un);
      return { e, o, v: parseFloat(v.replace(/,/g, "")) };
    });
    s.fig = mk("div", "abs mono", u, "FIG. 1 — RECONSTRUCTED FOLD · CHR22", { left: "1010px", top: "112px", font: "500 18px/1 var(--mono)", letterSpacing: ".1em", color: C.muted });
    const B = [["Turntable", 960, 160], ["Pause", 1132, 110], ["Iso", 1470, 76], ["Front", 1556, 98], ["Top", 1664, 76], ["Side", 1750, 90]];
    s.btns = B.map(([n, x, w]) => mk("div", "btn", u, n, { left: `${x}px`, top: "862px", width: `${w}px`, textAlign: "center", padding: "0" }));
    s.m.forEach((_, k) => cue(T.ws1 + 0.75 + k * 0.12, "pop"));
    cue(T.ws1 + 0.85, "roll", { d: 1.3 });
  },
  update(s, t) {
    GL.drawFold(t, "ws1", s.host);
    rise(s.kick, t, T.ws1 + 0.3);
    animSpans(s.h.ch, t, { in: T.ws1 + 0.35, st: 0.028, din: 0.5, fy: -150, fs: 1.4 });
    rise(s.lead, t, T.ws1 + 0.65);
    s.m.forEach((m, k) => {
      pop(m.e, t, T.ws1 + 0.75 + k * 0.12, { s: 0.5, r: k % 2 ? 6 : -6 });
      m.o.set(m.v * E.o3(lin(t, T.ws1 + 0.85, T.ws1 + 2.15)));
    });
    rise(s.fig, t, T.ws1 + 0.6);
    s.btns.forEach((b, k) => pop(b, t, T.ws1 + 0.5 + k * 0.05, { s: 0.7, y: 20 }));
    const on = t >= TURN;
    s.btns[0].style.background = on ? "#E4E6F7" : "#f8f8f6";
    s.btns[0].style.color = on ? C.cobalt : C.ink;
    s.btns[0].style.borderColor = on ? C.cobalt : "#d8d8d3";
  } });

// ================================================================================================= 05 4D dynamics (ws 02)
const PLAY = 30.05, SV0 = 30.2, SVD = 3.4;
track(29.5, 31.3, [[29.5, 1120, 1070], [30.0, 697, 907], [31.3, 760, 1010]], [PLAY]);
scene({ t0: T.ws2, t1: T.ws3 + 0.46, inT: { type: "iris", at: T.ws2, d: 0.5, x: 1300, y: 455 }, outT: { type: "push-left", at: T.ws3 },
  build(s) {
    bgDark(s, C.terra, 64, 44);
    const u = s.ui, S = D.sv;
    s.kick = kick(u, "WORKSPACE 02", 112, 150, "#FF9A72");
    s.h = headline(u, "4D DYNAMICS", 102, 188, { size: 150, color: "#fff", ext: C.terraDk, depth: 10 });
    s.lead = lead(u, "Watch the fold change over time, between conditions, or after a DNA rearrangement.", 112, 345, 640, "rgba(255,255,255,.82)");
    s.pill = mk("div", "pill", u, `${S.title} · ${S.deleted_mb} Mb removed`, { left: "112px", top: "470px", background: C.terra, color: "#fff", textTransform: "none", letterSpacing: ".02em", fontSize: "23px" });
    s.card = mk("div", "abs", u, null, { left: "112px", top: "560px", width: "430px", height: "190px", borderRadius: "20px", background: "rgba(255,255,255,.07)", border: "1px solid rgba(255,255,255,.14)", padding: "24px 28px" });
    mk("div", "mono", s.card, "LARGEST DISPLACEMENT", { font: "600 16px/1 var(--mono)", letterSpacing: ".14em", color: "rgba(255,255,255,.6)" });
    const row = mk("div", "", s.card, null, { display: "flex", alignItems: "baseline", gap: "12px", marginTop: "20px" });
    s.odo = odo(row, String(Math.round(S.max_disp_nm)), { fontSize: "100px", color: "#fff" });
    mk("div", "mono", row, "nm", { font: "600 30px/1 var(--mono)", color: "rgba(255,255,255,.7)" });
    s.note = lead(u, `Colour = how far each bead moved. The app’s 22q11.2 preset, ${S.frames} frames.`, 112, 780, 620, "rgba(255,255,255,.6)", 22);
    s.play = mk("div", "btn", u, "Play", { left: "640px", top: "880px", width: "114px", textAlign: "center", padding: "0", background: "rgba(255,255,255,.1)", color: "#fff", borderColor: "rgba(255,255,255,.2)" });
    s.pause = mk("div", "btn", u, "Pause", { left: "764px", top: "880px", width: "120px", textAlign: "center", padding: "0", background: "rgba(255,255,255,.1)", color: "#fff", borderColor: "rgba(255,255,255,.2)" });
    s.track = mk("div", "abs", u, null, { left: "920px", top: "904px", width: "880px", height: "6px", borderRadius: "3px", background: "rgba(255,255,255,.18)" });
    s.fill = mk("div", "abs", s.track, null, { left: "0", top: "0", width: "880px", height: "6px", borderRadius: "3px", background: C.terra, transformOrigin: "0 50%" });
    for (let k = 0; k < S.frames; k++) {
      mk("div", "abs", u, null, { left: `${920 + (880 * k) / (S.frames - 1) - 1}px`, top: "918px", width: "2px", height: "8px", background: "rgba(255,255,255,.3)" });
      if (k % 4 === 0 || k === S.frames - 1) mk("div", "abs mono", u, `${k * 15}`, { left: `${920 + (880 * k) / (S.frames - 1) - 20}px`, width: "40px", textAlign: "center", top: "934px", font: "500 14px/1 var(--mono)", color: "rgba(255,255,255,.45)" });
    }
    s.thumb = mk("div", "abs", u, null, { left: "907px", top: "894px", width: "26px", height: "26px", borderRadius: "50%", background: "#fff", boxShadow: `0 0 0 5px ${hexA(C.terra, 0.45)}` });
    s.val = mk("div", "abs mono", u, "", { left: "1500px", width: "300px", textAlign: "right", top: "862px", font: "600 20px/1 var(--mono)", color: "#fff" });
    cue(SV0, "rise", { d: SVD });
    cue(SV0 + SVD, "hit");
  },
  update(s, t) {
    GL.drawSV(t, SV0, SVD, s.host);
    rise(s.kick, t, T.ws2 + 0.3);
    animSpans(s.h.ch, t, { in: T.ws2 + 0.35, st: 0.028, din: 0.5, fy: -150, fs: 1.4 });
    rise(s.lead, t, T.ws2 + 0.6);
    pop(s.pill, t, T.ws2 + 0.75, { s: 0.6 });
    pop(s.card, t, T.ws2 + 0.9, { s: 0.7 });
    rise(s.note, t, T.ws2 + 1.1);
    const fp = clamp((t - SV0) / SVD), fr = Math.round(fp * (D.sv.frames - 1));
    s.odo.set(D.sv.max_disp_nm * E.o2(fp));
    [s.play, s.pause].forEach((b, k) => pop(b, t, T.ws2 + 0.5 + k * 0.08, { s: 0.7, y: 20 }));
    const on = t >= PLAY;
    s.play.style.background = on ? C.terra : "rgba(255,255,255,.1)";
    s.fill.style.transform = `scaleX(${fp.toFixed(4)})`;
    s.thumb.style.transform = `translateX(${(fp * 880).toFixed(1)}px)`;
    s.val.textContent = `${fr * 15} sweeps`;
    const pi = lin(t, T.ws2 + 0.5, T.ws2 + 0.9);
    [s.track, s.thumb, s.val].forEach((e) => (e.style.opacity = pi));
  } });

// ================================================================================================= 06 compare (ws 03)
track(35.6, 38.3, [[35.6, 900, 1070], [36.0, 190, 690], [37.9, 840, 690], [38.3, 900, 760]], []);
scene({ t0: T.ws3, t1: T.ws4 + 0.46, inT: { type: "push-left", at: T.ws3 }, outT: { type: "push-left", at: T.ws4 },
  build(s) {
    mk("div", "fill", s.bg, null, { background: C.ochre });
    mk("div", "fill halftone", s.bg, null, { opacity: 0.55 });
    const u = s.ui;
    s.kick = kick(u, "WORKSPACE 03", 112, 130, "rgba(28,30,27,.7)");
    s.h = headline(u, "COMPARE", 100, 168, { size: 160, color: C.ink, ext: C.ochreDk, depth: 11, soft: "rgba(0,0,0,.15)" });
    s.lead = lead(u, "Put two biological states side by side; rotating one rotates the other.", 112, 350, 560, C.ink);
    s.win = win(u, "ws03_top", 820, 120, 1020, 620, "localhost:8501 · ChronoCell-5D · 03 Compare");
    s.chart = card(u, "ws03_plot0", 112, 540, 760, 254, "Where they differ");
    s.scan = mk("div", "abs", u, null, { left: "0", top: "560px", width: "4px", height: "222px", borderRadius: "2px", background: C.cobalt, boxShadow: `0 0 0 4px ${hexA(C.cobalt, 0.2)}` });
    s.chip = mk("div", "pill", u, "Distance apart, along the region", { left: "112px", top: "820px", background: C.ink, color: "#fff", fontSize: "20px" });
    cue(T.ws3 + 1.2, "pop");
  },
  update(s, t) {
    rise(s.kick, t, T.ws3 + 0.25);
    animSpans(s.h.ch, t, { in: T.ws3 + 0.3, st: 0.04, din: 0.5, fy: -170, fs: 1.4 });
    rise(s.lead, t, T.ws3 + 0.6);
    const wp = E.oExpo(lin(t, T.ws3 + 0.3, T.ws3 + 1.0));
    setT(s.win, `translate3d(${((1 - wp) * 900).toFixed(1)}px,0,0) perspective(2000px) rotateY(${mix(-38, -12, wp) + 4 * lin(t, T.ws3, T.ws4)}deg) rotateX(4deg)`, wp);
    pop(s.chart, t, T.ws3 + 1.2, { s: 0.6, r: -6 });
    const sp = E.io3(lin(t, 36.0, 37.9));
    s.scan.style.left = `${(150 + sp * 690).toFixed(1)}px`;
    s.scan.style.opacity = lin(t, 35.95, 36.15) * (1 - lin(t, 38.0, 38.3));
    pop(s.chip, t, 36.1, { s: 0.6 });
  } });

// ================================================================================================= 07 drug lab (ws 04)
const DRAG0 = 39.8, DRAG1 = 41.6;
track(39.0, 42.8, [[39.0, 700, 1070], [39.6, 162, 657], [DRAG0, 162, 657], [DRAG1, 802, 657], [42.8, 900, 960]], [39.75], C.cobalt);
scene({ t0: T.ws4, t1: T.stack + 0.46, inT: { type: "push-left", at: T.ws4 }, outT: { type: "zoom", at: T.stack },
  build(s) {
    bgPaper(s, C.paper2);
    const u = s.ui;
    s.kick = kick(u, "WORKSPACE 04", 112, 130, C.terra);
    s.h = headline(u, "DRUG LAB", 100, 168, { size: 176, color: C.ink, ext: C.terra, depth: 11, soft: "rgba(0,0,0,.12)" });
    s.lead = lead(u, "Apply a virtual epigenetic drug and see how far it pushes the fold back toward healthy.", 112, 350, 640, C.ink2);
    s.sl = mk("div", "abs", u, null, { left: "112px", top: "520px", width: "760px", height: "200px", borderRadius: "22px", background: "#fff", boxShadow: "0 20px 50px rgba(0,0,0,.12)", padding: "26px 30px" });
    mk("div", "mono", s.sl, "DOSE · % OF MAXIMUM", { font: "600 17px/1 var(--mono)", letterSpacing: ".14em", color: C.muted });
    const vr = mk("div", "abs", s.sl, null, { right: "30px", top: "16px", display: "flex", alignItems: "baseline", gap: "6px", color: C.ink });
    s.odo = odo(vr, "100", { fontSize: "64px" });
    mk("div", "dn", vr, "%", { fontSize: "40px" });
    s.trk = mk("div", "abs", s.sl, null, { left: "50px", top: "130px", width: "640px", height: "10px", borderRadius: "5px", background: "#E3E3DE" });
    s.trkF = mk("div", "abs", s.trk, null, { left: "0", top: "0", width: "640px", height: "10px", borderRadius: "5px", background: C.cobalt, transformOrigin: "0 50%" });
    s.thumb = mk("div", "abs", s.sl, null, { left: "33px", top: "118px", width: "34px", height: "34px", borderRadius: "50%", background: "#fff", border: `4px solid ${C.cobalt}`, boxShadow: "0 4px 12px rgba(0,0,0,.2)" });
    s.chart = card(u, "ws04_plot1", 960, 130, 860, 454, "Dose–response");
    s.chartImg = s.chart.querySelector("img");
    s.fold = card(u, "ws04_plot0", 1250, 615, 560, 306, "Treated fold", { tagStyle: { background: C.terra } });
    for (let k = 1; k <= 10; k++) cue(mix(DRAG0, DRAG1, k / 10), "tick");
    cue(DRAG1 + 0.1, "pop");
  },
  update(s, t) {
    rise(s.kick, t, T.ws4 + 0.25);
    animSpans(s.h.ch, t, { in: T.ws4 + 0.3, st: 0.04, din: 0.5, fy: -170, fs: 1.4 });
    rise(s.lead, t, T.ws4 + 0.6);
    pop(s.sl, t, T.ws4 + 0.7, { s: 0.7 });
    pop(s.chart, t, T.ws4 + 0.85, { s: 0.7, r: 4 });
    const dose = E.io3(lin(t, DRAG0, DRAG1));
    s.odo.set(dose * 100);
    s.trkF.style.transform = `scaleX(${dose.toFixed(4)})`;
    s.thumb.style.transform = `translateX(${(dose * 640).toFixed(1)}px)`;
    s.chartImg.style.clipPath = `inset(0 ${((1 - dose) * 95).toFixed(2)}% 0 0)`;
    pop(s.fold, t, DRAG1 + 0.1, { s: 0.5, r: -8 });
  } });

// ================================================================================================= 08 genes + guide (ws 05, 08)
const PICK_GENES = 44.4, PICK_GUIDE = 45.8;
track(43.6, 47.3, [[43.6, 1900, 900], [44.3, 1560, 420], [45.3, 1560, 420], [45.7, 1560, 520], [47.3, 1700, 940]], [PICK_GENES, PICK_GUIDE], C.cobalt);
scene({ t0: T.stack, t1: T.quantum + 0.5, inT: { type: "zoom", at: T.stack }, outT: { type: "shrink", at: T.quantum },
  build(s) {
    bgPaper(s);
    const u = s.ui;
    s.kick = kick(u, "WORKSPACES 05 · 08", 112, 84, C.cobalt);
    s.h = headline(u, "GENES + GUIDE", 100, 120, { size: 120, color: C.ink, ext: C.cobalt, depth: 9, soft: "rgba(0,0,0,.12)" });
    s.iso = mk("div", "abs", u, null, { left: "250px", top: "330px", width: "900px", height: "506px", transformStyle: "preserve-3d", transform: "rotateX(55deg) rotateZ(-38deg)" });
    const names = ["ws06_pages", "ws06_quantum", "ws06_top", "ws05_fold", "ws05_top"];
    s.cards = names.map((n, k) => card(s.iso, n, 0, 0, 900, 506, null, { style: { borderRadius: "18px" } }));
    s.btns = [["05 Genes", 380], ["08 Guide", 480], ["Plain words", 580]].map(([n, y]) =>
      mk("div", "abs dc", u, n, { left: "1380px", top: `${y}px`, width: "380px", height: "80px", borderRadius: "999px", background: "#fff", color: C.ink, fontSize: "46px", lineHeight: "80px", textAlign: "center", boxShadow: "0 12px 30px rgba(0,0,0,.12)" }));
    s.hero = card(u, "ws06_top", 360, 210, 1200, 675, null);
    s.cap = mk("div", "pill", u, "What everything means, in plain words", { left: "360px", top: "915px", background: C.ink, color: "#fff", fontSize: "21px" });
    cue(PICK_GUIDE + 0.05, "whoosh");
  },
  update(s, t) {
    rise(s.kick, t, T.stack + 0.2);
    animSpans(s.h.ch, t, { in: T.stack + 0.25, st: 0.025, din: 0.5, fy: -130, fs: 1.4 });
    const g = E.oBack(lin(t, PICK_GENES, PICK_GENES + 0.45));
    s.cards.forEach((c, k) => {
      const a = E.oBack(lin(t, T.stack + 0.3 + k * 0.1, T.stack + 0.8 + k * 0.1));
      let z = k * 64 * a + (1 - clamp(a)) * -300, x = 0;
      if (k >= 3) { z += 150 * g * (k - 2); x += 80 * g; }
      c.style.transform = `translate3d(${x.toFixed(1)}px,0,${z.toFixed(1)}px)`;
      c.style.opacity = (clamp(a * 2) * (k === 2 ? 1 - lin(t, PICK_GUIDE, PICK_GUIDE + 0.08) : 1)).toFixed(3);
    });
    const sw = 1 + 0.03 * Math.sin((t - T.stack) * 1.6);
    s.iso.style.transform = `rotateX(${55 + 3 * Math.sin((t - T.stack) * 0.9)}deg) rotateZ(${-38 + 4 * lin(t, T.stack, T.quantum)}deg) scale(${sw})`;
    s.btns.forEach((b, k) => {
      pop(b, t, T.stack + 0.5 + k * 0.1, { s: 0.6, x: 60, y: 0 });
      const sel = (k === 0 && t >= PICK_GENES && t < PICK_GUIDE) || (k === 1 && t >= PICK_GUIDE);
      b.style.background = sel ? C.cobalt : "#fff";
      b.style.color = sel ? "#fff" : C.ink;
    });
    const hp = E.oExpo(lin(t, PICK_GUIDE + 0.05, PICK_GUIDE + 0.75));
    setT(s.hero, `translate3d(${((1 - hp) * -260).toFixed(1)}px,${((1 - hp) * 120).toFixed(1)}px,0) perspective(2000px) rotateX(${((1 - hp) * 50).toFixed(2)}deg) rotateZ(${((1 - hp) * -30).toFixed(2)}deg) scale(${mix(0.5, 1, hp).toFixed(4)})`, lin(t, PICK_GUIDE + 0.05, PICK_GUIDE + 0.2));
    pop(s.cap, t, PICK_GUIDE + 0.7, { s: 0.6 });
  } });

// ================================================================================================= 09 quantum lab (ws 06)
const QTABS = ["TAD boundaries", "Molecule (VQE)", "Noise &amp; mitigation", "Lattice fold", "Quantum walk", "Drug combination"];
const qx = (k) => 112 + k * 292 + 135;
const QCLICK = [48.5, 50.0, 51.3], QTAB = [3, 2, 4];
track(48.0, 52.3, [[48.0, 1300, 1070], [48.45, qx(3), 404], [49.5, qx(3), 404], [49.95, qx(2), 404], [50.8, qx(2), 404], [51.25, qx(4), 404], [52.3, 1500, 640]], QCLICK);
scene({ t0: T.quantum, t1: T.bento + 0.46, inT: { type: "iris", at: T.quantum, d: 0.5, x: 1560, y: 520 }, outT: { type: "push-up", at: T.bento },
  build(s) {
    bgDark(s, C.cobalt, 50, 70);
    const u = s.ui;
    s.ring = rings(s.bg, 960, 760, [260, 420, 600], "rgba(130,145,255,.14)");
    s.kick = kick(u, "WORKSPACE 06 · SIMULATED QUANTUM COMPUTER", 112, 96, "#8E98FF");
    s.h = headline(u, "QUANTUM LAB", 100, 132, { size: 150, color: "#fff", ext: C.cobaltDk, depth: 10 });
    s.lead = lead(u, "Try ChronoCell’s problems on a simulated quantum computer, next to the classical answer.", 112, 290, 1300, "rgba(255,255,255,.8)");
    s.tabs = QTABS.map((n, k) => mk("div", "tab", u, n, { left: `${112 + k * 292}px`, top: "375px", background: "rgba(255,255,255,.08)", color: "rgba(255,255,255,.85)", border: "1px solid rgba(255,255,255,.16)" }));
    s.panels = [
      [card(u, "q_lattice_plot0", 112, 480, 800, 364, "QAOA · simulated quantum"), card(u, "q_lattice_plot1", 1008, 480, 800, 364, "Exact optimum · classical", { tagStyle: { background: C.ink } })],
      [card(u, "q_noise_plot0", 112, 500, 1696, 382, "Noisy vs mitigated", { tagStyle: { background: C.terra } })],
      [card(u, "q_walk_plot0", 112, 470, 1696, 248, "Quantum walk"), card(u, "q_walk_plot1", 112, 738, 1696, 248, "Classical walk", { tagStyle: { background: C.ink } })],
    ];
  },
  update(s, t) {
    s.ring.style.transform = `rotate(${((t - T.quantum) * -10).toFixed(2)}deg)`;
    rise(s.kick, t, T.quantum + 0.3);
    animSpans(s.h.ch, t, { in: T.quantum + 0.3, st: 0.03, din: 0.5, fy: -150, fs: 1.4 });
    rise(s.lead, t, T.quantum + 0.55);
    let act = -1;
    QCLICK.forEach((c, i) => { if (t >= c) act = i; });
    s.tabs.forEach((b, k) => {
      pop(b, t, T.quantum + 0.45 + k * 0.06, { s: 0.6, y: 20 });
      const on = act >= 0 && QTAB[act] === k;
      b.style.background = on ? "#fff" : "rgba(255,255,255,.08)";
      b.style.color = on ? C.cobalt : "rgba(255,255,255,.85)";
    });
    s.panels.forEach((cards, i) => {
      const a = QCLICK[i] + 0.05, b = i < 2 ? QCLICK[i + 1] : 1e9;
      cards.forEach((c, k) => {
        const pin = E.oExpo(lin(t, a + k * 0.07, a + 0.5 + k * 0.07)), pout = E.i3(lin(t, b, b + 0.25));
        const vis = t >= a && t < b + 0.25;
        c.style.display = vis ? "block" : "none";
        if (vis) setT(c, `perspective(1800px) translate3d(0,${((1 - pin) * 140 - pout * 80).toFixed(1)}px,0) rotateX(${((1 - pin) * -60 + pout * 50).toFixed(2)}deg)`, clamp(pin * 1.5) * (1 - pout));
      });
      if (i < 2) cue(b + 0.02, "swish");
    });
  } });

// ================================================================================================= 10 everything (bento)
const TILES = [
  { t: "3D fold", s: "Whole chr22, bead by bead", bg: C.cobalt, fg: "#fff" },
  { t: "4D time-lapse", s: "Rearrangements, played frame by frame", bg: C.ink, fg: "#fff" },
  { t: "Drug lab", s: "A virtual epigenetic drug, dose by dose", bg: C.ochre, fg: C.ink },
  { t: "Genes", s: "Open, active chromatin vs buried", bg: "#fff", fg: C.ink },
  { t: "Quantum lab", s: "VQE, QAOA and quantum walks, on a simulated chip", bg: C.terra, fg: "#fff" },
  { t: "Scoreboard", s: "Every claim, tested once", bg: C.violet, fg: "#fff" },
];
scene({ t0: T.bento, t1: T.score + 0.5, inT: { type: "push-up", at: T.bento }, outT: { type: "shrink", at: T.score },
  build(s) {
    bgPaper(s, C.paper2);
    const u = s.ui;
    s.tiles = TILES.map((d, k) => {
      const e = mk("div", "tile", u, null, { left: `${80 + (k % 3) * 596}px`, top: `${96 + Math.floor(k / 3) * 448}px`, background: d.bg, color: d.fg });
      mk("div", "tt", e, d.t);
      mk("div", "ts", e, d.s);
      cue(T.bento + 0.4 + k * 0.11, "pop");
      return e;
    });
    const [t1, t2, t3, t4, t5, t6] = s.tiles;
    const disc = (p, name, x, y, sz) => { const c = mk("div", "abs", p, null, { left: `${x}px`, top: `${y}px`, width: `${sz}px`, height: `${sz}px`, borderRadius: "50%", overflow: "hidden", background: "#fff", boxShadow: "0 10px 30px rgba(0,0,0,.25)" }); const im = mk("img", "", c, null, { width: "100%", height: "100%", objectFit: "cover" }); im.src = img(name); return c; };
    s.d1 = disc(t1, "fold3d", 330, 160, 210);
    const r1 = mk("div", "abs", t1, null, { left: "36px", top: "268px" });
    s.o1 = odo(r1, fmt(F.beads), { fontSize: "74px" });
    mk("div", "mono", r1, "BEADS", { font: "700 20px/1 var(--mono)", letterSpacing: ".12em", marginTop: "10px" });
    s.ticks = [];
    for (let k = 0; k < D.sv.frames; k++) s.ticks.push(mk("div", "abs", t2, null, { left: `${36 + k * 20.6}px`, top: "300px", width: "14px", height: "64px", borderRadius: "4px", background: "rgba(255,255,255,.14)" }));
    const r2 = mk("div", "abs", t2, null, { left: "36px", top: "214px", display: "flex", alignItems: "baseline", gap: "10px" });
    s.o2 = odo(r2, String(D.sv.frames), { fontSize: "70px" });
    mk("div", "mono", r2, "FRAMES", { font: "700 20px/1 var(--mono)", letterSpacing: ".12em" });
    const tr3 = mk("div", "abs", t3, null, { left: "36px", top: "318px", width: "496px", height: "12px", borderRadius: "6px", background: "rgba(28,30,27,.18)" });
    s.f3 = mk("div", "abs", tr3, null, { left: "0", top: "0", width: "496px", height: "12px", borderRadius: "6px", background: C.ink, transformOrigin: "0 50%" });
    s.th3 = mk("div", "abs", t3, null, { left: "22px", top: "306px", width: "36px", height: "36px", borderRadius: "50%", background: "#fff", border: `5px solid ${C.ink}` });
    const r3 = mk("div", "abs", t3, null, { left: "36px", top: "214px", display: "flex", alignItems: "baseline", gap: "6px" });
    s.o3 = odo(r3, "100", { fontSize: "64px" });
    mk("div", "dn", r3, "% dose", { fontSize: "34px" });
    s.d4 = disc(t4, "genes3d", 310, 140, 220);
    s.p4 = [["Open", C.terra], ["Buried", C.cobalt]].map(([n, c], k) => mk("div", "pill", t4, n, { left: `${36 + k * 150}px`, top: "300px", border: `3px solid ${c}`, color: c, fontSize: "20px", padding: "12px 20px" }));
    const m5 = mk("div", "abs", t5, null, { left: "300px", top: "160px", width: "232px", height: "106px", borderRadius: "14px", overflow: "hidden", background: "#fff", boxShadow: "0 10px 30px rgba(0,0,0,.25)" });
    const im5 = mk("img", "", m5, null, { width: "100%", height: "100%", objectFit: "cover" });
    im5.src = img("q_lattice_plot0");
    s.q5 = [];
    for (let k = 0; k < 8; k++) s.q5.push(mk("div", "abs", t5, null, { left: `${36 + k * 62}px`, top: "318px", width: "40px", height: "40px", borderRadius: "50%", border: "4px solid #fff" }));
    const r6 = mk("div", "abs", t6, null, { left: "36px", top: "200px", display: "flex", alignItems: "baseline", gap: "10px" });
    s.o6 = odo(r6, String(TS.n), { fontSize: "92px" });
    mk("div", "mono", r6, "TESTS", { font: "700 20px/1 var(--mono)", letterSpacing: ".12em" });
    const bar = mk("div", "abs", t6, null, { left: "36px", top: "320px", width: "496px", height: "34px", borderRadius: "8px", overflow: "hidden", display: "flex", background: "rgba(255,255,255,.15)" });
    s.seg = [["pass", TS.pass], ["fail", TS.fail], ["other", TS.other]].map(([k, n]) => mk("div", "", bar, null, { height: "34px", width: `${(496 * n) / TS.n}px`, background: STATUS[k], transformOrigin: "0 50%" }));
    s.seglab = mk("div", "abs mono", t6, `${TS.pass} passed · ${TS.fail} failed · ${TS.other} mixed`, { left: "36px", top: "368px", font: "600 17px/1 var(--mono)", color: "rgba(255,255,255,.85)" });
    cue(T.bento + 0.8, "roll", { d: 1.3 });
  },
  update(s, t) {
    const lt = t - T.bento;
    s.tiles.forEach((e, k) => pop(e, t, T.bento + 0.4 + k * 0.11, { s: 0.4, r: k % 2 ? 7 : -7, y: 60, d: 0.55 }));
    const roll = E.o3(lin(t, T.bento + 0.8, T.bento + 2.1));
    s.o1.set(F.beads * roll);
    s.d1.style.transform = `rotate(${(lt * 25).toFixed(2)}deg)`;
    s.o2.set(D.sv.frames * roll);
    const on = Math.floor(((lt - 0.8) * 12) % (D.sv.frames + 6));
    s.ticks.forEach((e, k) => (e.style.background = lt > 0.8 && k <= on ? C.terra : "rgba(255,255,255,.14)"));
    const dose = 0.5 - 0.5 * Math.cos(Math.max(0, lt - 0.8) * 2.2);
    s.f3.style.transform = `scaleX(${dose.toFixed(4)})`;
    s.th3.style.transform = `translateX(${(dose * 496).toFixed(1)}px)`;
    s.o3.set(dose * 100);
    const which = Math.floor(Math.max(0, lt - 0.8) / 0.8) % 2;
    s.p4.forEach((p, k) => { const c = k ? C.cobalt : C.terra; p.style.background = which === k ? c : "transparent"; p.style.color = which === k ? "#fff" : c; });
    s.q5.forEach((q, k) => { const a = 0.5 + 0.5 * Math.sin(lt * 5 - k * 0.7); q.style.background = `rgba(255,255,255,${(a * 0.9).toFixed(3)})`; q.style.transform = `scale(${(0.8 + 0.25 * a).toFixed(3)})`; });
    s.o6.set(TS.n * roll);
    s.seg.forEach((g, k) => (g.style.transform = `scaleX(${E.o3(lin(t, T.bento + 1.0 + k * 0.25, T.bento + 1.6 + k * 0.25)).toFixed(4)})`));
    s.seglab.style.opacity = lin(t, T.bento + 1.8, T.bento + 2.1);
  } });

// ================================================================================================= 11 the scoreboard (ws 07)
const NC = [1250, 560];
const AREAS = Object.keys(TS.areas);
const ASHORT = { "Structure & imaging": "STRUCTURE", "Analysis & perturbations": "ANALYSIS", "Quantum lab": "QUANTUM", "Drug lab": "DRUG LAB" };
// the tests sort into four rows: passed; failed, then fixed by a retest (new method, new data); failed and still open; mixed
const rowOf = (r) => (r.status === "fail" ? (r.fixed_by ? "fixed" : "open") : r.status);
const NET = (() => {
  const nodes = [], hubs = [];
  let a0 = -Math.PI / 2;
  const rowCount = { pass: 0, fixed: 0, open: 0, other: 0 };
  AREAS.forEach((area) => {
    const rows = TS.rows.filter((r) => r.area === area);
    const span = (2 * Math.PI * rows.length) / TS.n;
    hubs.push({ area, ang: a0 + span / 2 });
    rows.forEach((r, k) => {
      const ang = a0 + (span * (k + 0.5)) / rows.length;
      const row = rowOf(r), slot = rowCount[row]++;
      nodes.push({ ...r, ang, rad: k % 2 ? 405 : 330, hub: hubs.length - 1, row, slot });
    });
    a0 += span;
  });
  nodes.forEach((n, i) => (n.order = i));
  return { nodes, hubs };
})();
const ROWY = { pass: 452, fixed: 578, open: 704, other: 830 };
const STAMP_SB = 64.0;
scene({ t0: T.score, t1: T.loops + 0.46, inT: { type: "iris", at: T.score, d: 0.5 }, outT: { type: "push-left", at: T.loops },
  build(s) {
    bgDark(s, C.cobalt, 64, 52);
    const u = s.ui;
    s.cv = mk("canvas", "c2d", u);
    s.cv.width = W; s.cv.height = H;
    s.g = s.cv.getContext("2d");
    const h = mk("div", "abs dx", u, null, { left: "100px", top: "96px", fontSize: "78px", color: "#fff", textShadow: extrude(C.cobaltDk, 9) });
    s.l1 = chars(mk("div", "", h), "EVERY CLAIM");
    s.l2 = chars(mk("div", "", h, null, { color: C.ochre, textShadow: extrude(C.ochreDk, 9) }), "IS A TEST");
    s.sub = lead(u, "Written down before its data were read. Run once. Kept on the record, pass or fail.", 112, 290, 700, "rgba(255,255,255,.8)", 25);
    s.cnt = [["pass", TS.pass, "passed", C.pass], ["fixed", TS.fixed, "failed · then passed as a retest", C.fail],
      ["open", TS.open, "failed · still open", C.fail], ["other", TS.other, "mixed · blocked · baseline", C.other]].map(([k, n, l, col]) => {
      const e = mk("div", "abs", u, null, { left: "112px", top: `${ROWY[k] - 50}px`, display: "flex", alignItems: "center", gap: "18px" });
      const o = odo(e, String(n), { fontSize: "96px", color: col });
      mk("div", "kick", e, k === "fixed" ? `failed · <span style="color:${C.pass}">then passed as a retest</span>` : l,
        { color: col, width: "250px", lineHeight: "1.3" });
      return { e, o, n };
    });
    s.stamp = stamp(u, "KEPT ON THE RECORD", 980, 600, -6, { color: C.fail, bg: "rgba(17,19,24,.88)" });
    s.win = win(u, "ws08_top", 260, 150, 1400, 834, "localhost:8501 · ChronoCell-5D · 07 Scoreboard");
    s.wtag = mk("div", "pill", u, "Workspace 07 · Scoreboard", { left: "260px", top: "90px", background: C.ochre, color: C.ink, fontSize: "21px" });
    cue(T.score + 0.45, "hit");
    shake(T.score + 0.5, 8);
    cue(60.1, "rise", { d: 0.9 });
    cue(61.6, "whoosh");
    cue(62.3, "roll", { d: 1.1 });
    cue(STAMP_SB, "stamp");
    shake(STAMP_SB + 0.05, 18);
    cue(65.1, "whoosh");
  },
  update(s, t) {
    const g = s.g;
    g.clearRect(0, 0, W, H);
    const rot = 0.05 * (t - T.score);
    const fly = (n) => E.io3(lin(t, 61.6 + n.order * 0.012, 62.4 + n.order * 0.012));
    const netFade = 1 - lin(t, 61.6, 62.1);
    const pos = (ang, rad, p) => [NC[0] + Math.cos(ang + rot) * rad * p, NC[1] + Math.sin(ang + rot) * rad * p];
    const hubP = NET.hubs.map((hb, k) => { const p = E.oBack(lin(t, T.score + 0.45 + k * 0.05, T.score + 1.0 + k * 0.05)); return { p, xy: pos(hb.ang, 175, p) }; });
    // lines
    if (netFade > 0) {
      g.lineWidth = 2;
      hubP.forEach((h) => { g.strokeStyle = `rgba(255,255,255,${(0.3 * clamp(h.p) * netFade).toFixed(3)})`; g.beginPath(); g.moveTo(NC[0], NC[1]); g.lineTo(h.xy[0], h.xy[1]); g.stroke(); });
    }
    const leaf = NET.nodes.map((n) => {
      const p = E.oBack(lin(t, T.score + 0.6 + n.order * 0.018, T.score + 1.25 + n.order * 0.018));
      const [x0, y0] = pos(n.ang, n.rad, p);
      const f = fly(n), tx = 660 + n.slot * 68, ty = ROWY[n.row];          // 15 in a row stay behind the window (x < 1660)
      const x = mix(x0, tx, f), y = mix(y0, ty, f) - Math.sin(Math.PI * f) * 90;
      return { n, p, x, y, f };
    });
    if (netFade > 0) leaf.forEach(({ n, p, x, y }) => {
      const hb = hubP[n.hub].xy;
      g.strokeStyle = `rgba(255,255,255,${(0.18 * clamp(p) * netFade).toFixed(3)})`;
      g.beginPath(); g.moveTo(hb[0], hb[1]); g.lineTo(x, y); g.stroke();
    });
    // centre and hubs
    const cp = E.oBack(lin(t, T.score + 0.35, T.score + 0.85));
    if (netFade > 0 && cp > 0) {
      g.globalAlpha = netFade;
      g.fillStyle = C.cobalt;
      g.beginPath(); g.arc(NC[0], NC[1], 58 * clamp(cp, 0, 1.2), 0, 6.2832); g.fill();
      g.strokeStyle = "rgba(255,255,255,.5)"; g.lineWidth = 3; g.stroke();
      g.fillStyle = "#fff"; g.font = "700 20px 'IBM Plex Mono'"; g.textAlign = "center"; g.textBaseline = "middle";
      g.fillText(`${TS.n}`, NC[0], NC[1] - 9);
      g.font = "600 12px 'IBM Plex Mono'"; g.fillText("TESTS", NC[0], NC[1] + 13);
      hubP.forEach((h, k) => {
        if (h.p <= 0) return;
        g.fillStyle = "#262a36"; g.beginPath(); g.arc(h.xy[0], h.xy[1], 30 * clamp(h.p, 0, 1.2), 0, 6.2832); g.fill();
        g.strokeStyle = "rgba(255,255,255,.45)"; g.lineWidth = 2; g.stroke();
        g.fillStyle = "rgba(255,255,255,.8)"; g.font = "700 13px 'IBM Plex Mono'";
        g.fillText(ASHORT[NET.hubs[k].area] || NET.hubs[k].area, h.xy[0], h.xy[1] + 50);
      });
      g.globalAlpha = 1;
    }
    // leaves
    g.textAlign = "center"; g.textBaseline = "middle";
    leaf.forEach(({ n, p, x, y, f }, i) => {
      if (p <= 0) return;
      const ct = 60.1 + (i / NET.nodes.length) * 0.9, cl = lin(t, ct, ct + 0.12);
      const col = STATUS[n.status];
      const r = 24 * clamp(p, 0, 1.25) * (1 + 0.35 * Math.exp(-Math.pow((t - ct - 0.06) / 0.08, 2)));
      g.globalAlpha = clamp(p * 2);
      g.fillStyle = cl > 0.5 ? col : "#ECEEF7";
      g.beginPath(); g.arc(x, y, r, 0, 6.2832); g.fill();
      if (n.row === "fixed" && f > 0.6) {                    // a fail that a later retest passed: a green ring
        g.strokeStyle = C.pass; g.lineWidth = 5; g.globalAlpha = lin(f, 0.6, 1);
        g.beginPath(); g.arc(x, y, r + 5, 0, 6.2832); g.stroke(); g.globalAlpha = clamp(p * 2);
      }
      if (cl > 0 && t < ct + 0.5) { g.strokeStyle = col; g.lineWidth = 3; g.globalAlpha = 1 - lin(t, ct, ct + 0.5); g.beginPath(); g.arc(x, y, r + 6 + 26 * lin(t, ct, ct + 0.5), 0, 6.2832); g.stroke(); g.globalAlpha = 1; }
      g.fillStyle = cl > 0.5 ? "#fff" : C.ink;
      g.font = `700 ${n.id.length > 3 ? 12 : 14}px 'IBM Plex Mono'`;
      g.fillText(n.id, x, y + 1);
    });
    g.globalAlpha = 1;
    animSpans(s.l1, t, { in: T.score + 0.3, st: 0.03, din: 0.5, fy: -120, fs: 1.4 });
    animSpans(s.l2, t, { in: T.score + 0.55, st: 0.03, din: 0.5, fy: -120, fs: 1.4 });
    rise(s.sub, t, T.score + 0.9);
    s.cnt.forEach((c, k) => { pop(c.e, t, 62.2 + k * 0.12, { s: 0.6, x: -40, y: 0 }); c.o.set(c.n * E.o3(lin(t, 62.3, 63.4))); });
    stampAnim(s.stamp, t, STAMP_SB);
    const wp = E.oExpo(lin(t, 65.1, 65.8));
    setT(s.win, `perspective(2000px) translate3d(0,${((1 - wp) * 900).toFixed(1)}px,0) rotateX(${((1 - wp) * 28).toFixed(2)}deg)`, lin(t, 65.1, 65.25));
    pop(s.wtag, t, 65.6, { s: 0.6 });
  } });

// ================================================================================================= 12 DNA loops (gate 6b)
const G0 = 67.8, SPD = 0.33;
const LOOPS = [{ key: "hmec", x: 112, title: "HMEC · breast epithelium" }, { key: "hap1", x: 1000, title: "HAP-1 · near-haploid line" }];
const TOOLS = [["mustache", "Mustache", C.ochre], ["chromosight", "chromosight", C.other], ["chronocell", "ChronoCell", C.cobaltHi]];
const SLOT = [395, 515, 635];
scene({ t0: T.loops, t1: T.mol + 0.5, inT: { type: "push-left", at: T.loops }, outT: { type: "shrink", at: T.mol },
  build(s) {
    bgDark(s, C.cobalt, 20, 18);
    const u = s.ui;
    mk("div", "abs dc", s.bg, "Gate 6b", { right: "50px", top: "760px", fontSize: "300px", color: "rgba(255,255,255,.05)", textTransform: "none" });
    s.icon = mk("div", "abs", u, `<svg width="92" height="92" viewBox="0 0 92 92"><rect width="92" height="92" rx="20" fill="${C.cobalt}"/><path d="M30 58c-9-9-9-24 2-31s27-3 29 10" fill="none" stroke="#fff" stroke-width="7" stroke-linecap="round"/><path d="M52 31l10 6-2-12" fill="none" stroke="#fff" stroke-width="7" stroke-linecap="round" stroke-linejoin="round"/></svg>`, { left: "112px", top: "100px" });
    s.h = headline(u, "TOP LOOP CALLERS", 226, 96, { size: 112, color: "#fff", ext: C.cobaltDk, depth: 8 });
    s.sub = mk("div", "abs mono", u, "F1 against ENCODE reference loops · cell types it had never seen · higher is better", { left: "114px", top: "232px", font: "500 21px/1 var(--mono)", letterSpacing: ".04em", color: "rgba(255,255,255,.62)" });
    s.panels = LOOPS.map((P) => {
      const v = HI.loops[P.key];
      const ttl = mk("div", "abs kick", u, P.title, { left: `${P.x}px`, top: "318px", color: "rgba(255,255,255,.85)", letterSpacing: ".12em" });
      SLOT.forEach((y, k) => mk("div", "abs mono", u, `#${k + 1}`, { left: `${P.x}px`, top: `${y + 26}px`, font: "700 28px/1 var(--mono)", color: "rgba(255,255,255,.45)" }));
      const rows = TOOLS.map(([key, name, col], k) => {
        const row = mk("div", "abs", u, null, { left: `${P.x + 64}px`, top: "0", width: "800px", height: "96px" });
        mk("div", "abs", row, `<span style="display:inline-block;width:14px;height:14px;border-radius:50%;background:${col};margin-right:12px"></span>${name}`,
          { left: "0", top: "16px", height: "64px", width: "240px", borderRadius: "14px", background: "rgba(255,255,255,.08)", font: "800 28px/64px var(--display)", fontStretch: "90%", padding: "0 18px", color: "#fff", whiteSpace: "nowrap" });
        const bar = mk("div", "abs", row, null, { left: "256px", top: "16px", height: "64px", width: `${(v[key] * 540).toFixed(1)}px`, borderRadius: "12px", background: col, transformOrigin: "0 50%",
          boxShadow: key === "chronocell" ? `0 0 34px ${hexA(C.cobaltHi, 0.6)}` : "none" });
        const val = mk("div", "abs dn", row, "", { top: "26px", fontSize: "44px", color: "#fff" });
        return { key, row, bar, val, v: v[key], slot0: k };
      });
      const crown = mk("div", "abs", u, CROWN, { left: "0", top: "0" });
      return { P, ttl, rows, crown, v };
    });
    s.panels.forEach((pn) => {
      const tf = G0 + pn.v.chronocell / SPD;
      const x = pn.P.x + 64 + 256 + pn.v.chronocell * 540;
      burst(tf + 0.05, x, SLOT[0] + 48, { n: 60, speed: 800, life: 1.6, colors: [C.ochre, "#fff", C.cobaltHi] });
    });
    cue(G0, "roll", { d: 2.2 });
    cue(G0 + 0.47 / SPD, "whoosh");
  },
  update(s, t) {
    pop(s.icon, t, T.loops + 0.3, { s: 0.3, r: -30, y: 0 });
    animSpans(s.h.ch, t, { in: T.loops + 0.35, st: 0.022, din: 0.5, fy: -120, fs: 1.4 });
    rise(s.sub, t, T.loops + 0.6);
    const sw = E.io3(lin(t, G0 + 0.47 / SPD, G0 + 0.47 / SPD + 0.4));
    s.panels.forEach((pn, pi) => {
      rise(pn.ttl, t, T.loops + 0.55 + pi * 0.1);
      pn.rows.forEach((r) => {
        const cur = Math.min(r.v, SPD * Math.max(0, t - G0));
        const slot = r.key === "chronocell" ? mix(2, 0, sw) : r.key === "mustache" ? mix(0, 1, sw) : mix(1, 2, sw);
        const y = mix(SLOT[Math.floor(slot)], SLOT[Math.min(2, Math.floor(slot) + 1)], slot - Math.floor(slot));
        pop(r.row, t, T.loops + 0.6 + r.slot0 * 0.08 + pi * 0.1, { s: 0.8, x: -60, y: 0 });
        r.row.style.top = `${y.toFixed(1)}px`;
        r.bar.style.transform = `scaleX(${(cur / r.v).toFixed(4)})`;
        r.val.textContent = cur.toFixed(2);
        r.val.style.left = `${(256 + cur * 540 + 20).toFixed(1)}px`;
        r.val.style.opacity = lin(t, G0, G0 + 0.2);
      });
      const tf = G0 + pn.v.chronocell / SPD, x = pn.P.x + 64 + 256 + pn.v.chronocell * 540;
      pn.crown.style.left = `${x + 30}px`; pn.crown.style.top = `${SLOT[0] - 38}px`;
      pop(pn.crown, t, tf, { s: 0.2, r: -30, y: -40 });
    });
  } });

// ================================================================================================= 13 molecules (gate Q6d)
const Q6 = HI.q6d, CHIP0 = 73.0, CHIPD = 0.11, STAMP_Q6 = 75.0;
scene({ t0: T.mol, t1: T.mit + 0.46, inT: { type: "iris", at: T.mol, d: 0.5, x: 960, y: 540, sound: "boom" }, outT: { type: "push-up", at: T.mit },
  build(s) {
    s.burst = bgBurst(s, C.terra);
    const u = s.ui;
    s.disc = mk("div", "abs", u, CHECK(230, C.terra, 6), { left: "112px", top: "150px", width: "290px", height: "290px", borderRadius: "50%", background: "#fff", display: "flex", alignItems: "center", justifyContent: "center", boxShadow: `14px 18px 0 ${C.terraDk}, 0 40px 60px rgba(0,0,0,.25)` });
    const r = mk("div", "abs", u, null, { left: "450px", top: "150px", display: "flex", alignItems: "baseline" });
    s.odo = odo(r, String(Q6.cases), { fontSize: "250px", color: "#fff", textShadow: extrude(C.terraDk, 13) });
    mk("div", "dn", r, `/${Q6.cases}`, { fontSize: "250px", color: "rgba(255,255,255,.55)", textShadow: extrude(hexA(C.terraDk, 0.6), 13) });
    s.row = r;
    s.lab = kick(u, "new molecules within chemical accuracy", 118, 470, "#fff");
    s.sub = lead(u, "Bonds stretched toward breaking, solved on a simulated quantum computer (ADAPT-VQE). × = how far each bond is stretched.", 118, 520, 900, "rgba(255,255,255,.9)", 25);
    s.chips = Q6.rows.map((m, k) => {
      const e = mk("div", "chip", u, `<span class="f">${subH(m.m)}</span><span class="s">${m.s}×</span><span class="ok">✓</span>`,
        { left: `${1180 + (k % 2) * 320}px`, top: `${120 + Math.floor(k / 2) * 112}px`, width: "300px" });
      cue(CHIP0 + k * CHIPD, "tick");
      return e;
    });
    s.stamp = stamp(u, `WORST ${Q6.worst} mHa · LIMIT ${HI.chem_accuracy} mHa`, 130, 690, -5, { color: C.ink, bg: C.ochre, border: C.ink });
    burst(74.55, 640, 300, { n: 150, speed: 1300, colors: ["#fff", C.ochre, C.ink, C.cobaltHi] });
    cue(STAMP_Q6, "stamp");
    shake(STAMP_Q6 + 0.05, 16);
  },
  update(s, t) {
    s.burst.style.transform = `rotate(${(t * -8).toFixed(2)}deg)`;
    pop(s.disc, t, T.mol + 0.3, { s: 0.2, r: -60, y: 0, d: 0.55 });
    s.disc.style.transform += ` rotate(${(3 * Math.sin((t - T.mol) * 3)).toFixed(2)}deg)`;
    pop(s.row, t, T.mol + 0.35, { s: 0.6, y: 40 });
    s.odo.set(t < CHIP0 ? 0 : Math.min(Q6.cases, (t - CHIP0) / CHIPD + 1));
    rise(s.lab, t, T.mol + 0.6);
    rise(s.sub, t, T.mol + 0.75);
    s.chips.forEach((c, k) => pop(c, t, CHIP0 + k * CHIPD, { s: 0.3, r: k % 2 ? 8 : -8, y: 40, d: 0.4 }));
    stampAnim(s.stamp, t, STAMP_Q6);
  } });

// ================================================================================================= 14 noisy chips (gate Q9)
const Q9 = HI.q9, FALL0 = 77.8, STAMP_Q9 = 79.1;
const LY = (e) => 860 - ((Math.log10(Math.max(e, 0.1)) + 1) / (Math.log10(30) + 1)) * 600;
scene({ t0: T.mit, t1: T.fail + 0.5, inT: { type: "push-up", at: T.mit }, outT: { type: "shrink", at: T.fail },
  build(s) {
    bgDark(s, C.pass, 72, 55);
    const u = s.ui;
    const h = mk("div", "abs dc", u, null, { left: "100px", top: "100px", fontSize: "124px" });
    s.l1 = chars(mk("div", "", h, null, { color: "#fff", textShadow: extrude("#000", 9) }), "A NOISY CHIP,");
    s.l2 = chars(mk("div", "", h, null, { color: C.pass, textShadow: extrude(C.passDk, 9) }), "REPAIRED");
    const r = mk("div", "abs", u, null, { left: "108px", top: "420px", display: "flex", alignItems: "baseline", gap: "16px" });
    s.odo = odo(r, Q9.noisy.toFixed(2).padStart(5, "0"), { fontSize: "180px" });
    mk("div", "mono", r, "mHa", { font: "600 44px/1 var(--mono)", color: "rgba(255,255,255,.75)" });
    s.row = r;
    s.lab = kick(u, `median energy error · ${Q9.cases} new molecules`, 118, 640, "rgba(255,255,255,.75)");
    s.sub = lead(u, "Simulated at today’s best noise level, then corrected with symmetry checks and zero-noise extrapolation.", 118, 690, 700, "rgba(255,255,255,.75)", 24);
    const svgW = 940, svgH = 760;
    s.svg = mk("div", "abs", u, `<svg width="${svgW}" height="${svgH}" viewBox="900 200 ${svgW} ${svgH}"></svg>`, { left: "900px", top: "200px" }).firstChild;
    const NS = "http://www.w3.org/2000/svg";
    const el = (tag, at) => { const e = document.createElementNS(NS, tag); for (const k in at) e.setAttribute(k, at[k]); s.svg.appendChild(e); return e; };
    [0.1, 1, 10].forEach((v) => {
      el("line", { x1: 960, x2: 1820, y1: LY(v), y2: LY(v), stroke: "rgba(255,255,255,.12)", "stroke-width": 1 });
      const tx = el("text", { x: 948, y: LY(v) + 6, fill: "rgba(255,255,255,.5)", "font-family": "IBM Plex Mono", "font-size": 17, "text-anchor": "end" });
      tx.textContent = String(v);
    });
    s.lim = el("line", { x1: 960, x2: 1820, y1: LY(HI.chem_accuracy), y2: LY(HI.chem_accuracy), stroke: C.ochre, "stroke-width": 3, "stroke-dasharray": "12 9" });
    s.limT = el("text", { x: 1820, y: LY(HI.chem_accuracy) - 12, fill: C.ochre, "font-family": "IBM Plex Mono", "font-size": 18, "font-weight": 600, "text-anchor": "end" });
    s.limT.textContent = `CHEMICAL ACCURACY · ${HI.chem_accuracy} mHa`;
    s.dots = Q9.rows.map((m, k) => {
      const x = 990 + k * 52;
      const trail = el("line", { x1: x, x2: x, y1: LY(m.noisy), y2: LY(m.noisy), stroke: "rgba(255,255,255,.25)", "stroke-width": 2, "stroke-dasharray": "4 5" });
      const c = el("circle", { cx: x, cy: LY(m.noisy), r: 13, fill: C.fail });
      const lab = el("text", { x, y: 892, fill: "rgba(255,255,255,.6)", "font-family": "IBM Plex Mono", "font-size": 14, "text-anchor": "end", transform: `rotate(-55 ${x} 892)` });
      lab.textContent = `${subU(m.m)} ${m.s}×`;
      return { m, x, trail, c, lab, k };
    });
    s.stamp = stamp(u, `${Q9.reduction}× SMALLER`, 1150, 210, -6, { color: C.pass, bg: "rgba(17,19,24,.88)" });
    cue(FALL0, "roll", { d: 1.2 });
    cue(STAMP_Q9, "stamp");
    shake(STAMP_Q9 + 0.05, 14);
  },
  update(s, t) {
    animSpans(s.l1, t, { in: T.mit + 0.3, st: 0.03, din: 0.5, fy: -140, fs: 1.4 });
    animSpans(s.l2, t, { in: T.mit + 0.55, st: 0.035, din: 0.5, fy: -140, fs: 1.4 });
    const v = mix(Q9.noisy, Q9.median, E.io3(lin(t, FALL0, FALL0 + 1.2)));
    s.odo.set(v);
    s.row.style.color = v < HI.chem_accuracy ? C.pass : C.fail;
    pop(s.row, t, T.mit + 0.5, { s: 0.6 });
    rise(s.lab, t, T.mit + 0.7);
    rise(s.sub, t, T.mit + 0.85);
    s.dots.forEach((d) => {
      const a = pop(d.c, t, T.mit + 0.6 + d.k * 0.04, { s: 0.2, y: 0 });
      d.c.style.transform = "none";
      d.c.setAttribute("r", (13 * clamp(a, 0, 1.3)).toFixed(2));
      const q = E.oBack(lin(t, FALL0 + d.k * 0.06, FALL0 + 0.55 + d.k * 0.06));
      const y = mix(LY(d.m.noisy), LY(d.m.fixed), q);
      d.c.setAttribute("cy", y.toFixed(1));
      d.c.setAttribute("fill", y > LY(HI.chem_accuracy) ? C.pass : C.fail);
      d.trail.setAttribute("y2", y.toFixed(1));
      d.trail.style.opacity = lin(t, FALL0 + d.k * 0.06, FALL0 + 0.2 + d.k * 0.06);
      d.lab.style.opacity = lin(t, T.mit + 0.6, T.mit + 1.0);
    });
    s.lim.style.opacity = lin(t, T.mit + 0.5, T.mit + 0.8);
    s.limT.style.opacity = lin(t, T.mit + 0.6, T.mit + 0.9);
    stampAnim(s.stamp, t, STAMP_Q9);
  } });

// ================================================================================================= 15 the honest fail (gate Q7b)
const Q7 = HI.q7b, GROW0 = 80.7, STAMP_F = 82.4;
scene({ t0: T.fail, t1: T.words + 0.02, inT: { type: "iris", at: T.fail, d: 0.5, sound: "low" },
  build(s) {
    bgDark(s, C.terra, 40, 30);
    const u = s.ui;
    const h = mk("div", "abs dx", u, null, { left: "100px", top: "104px", fontSize: "88px" });
    s.l1 = chars(mk("div", "", h, null, { color: "#fff", textShadow: extrude("#000", 8) }), "AND WHERE IT LOST,");
    s.l2 = chars(mk("div", "", h, null, { color: C.terra, textShadow: extrude(C.terraDk, 8) }), "IT SAYS SO.");
    s.sub = lead(u, `On ${Q7.complexes} protein–drug complexes it had never seen, the quantum route placed ${Q7.qaoa} % of the drugs correctly; plain random search placed ${Q7.random} %. Gate Q7b.`, 112, 316, 1300, "rgba(255,255,255,.8)", 26);
    s.bars = [["Quantum docking (QAOA)", Q7.qaoa, C.terra, 470], ["Random search, same score", Q7.random, "#fff", 680]].map(([l, v, col, y]) => {
      const lab = kick(u, l, 114, y, "rgba(255,255,255,.8)");
      const bar = mk("div", "abs", u, null, { left: "112px", top: `${y + 40}px`, width: `${v * 14}px`, height: "100px", borderRadius: "16px", background: col, transformOrigin: "0 50%" });
      const r = mk("div", "abs", u, null, { top: `${y + 46}px`, display: "flex", alignItems: "baseline", gap: "6px", color: "#fff" });
      const o = odo(r, String(v), { fontSize: "86px" });
      mk("div", "dn", r, "%", { fontSize: "54px" });
      return { lab, bar, r, o, v };
    });
    s.stamp = stamp(u, "KEPT ON THE RECORD", 1090, 560, -6, { color: "#fff", bg: C.terra, border: "#fff" });
    cue(GROW0, "roll", { d: 1.2 });
    cue(STAMP_F, "stamp");
    shake(STAMP_F + 0.05, 18);
  },
  update(s, t) {
    animSpans(s.l1, t, { in: T.fail + 0.3, st: 0.022, din: 0.5, fy: -120, fs: 1.4 });
    animSpans(s.l2, t, { in: T.fail + 0.6, st: 0.03, din: 0.5, fy: -120, fs: 1.4 });
    rise(s.sub, t, T.fail + 0.75);
    s.bars.forEach((b, k) => {
      rise(b.lab, t, GROW0 - 0.2 + k * 0.15);
      const p = E.o3(lin(t, GROW0 + k * 0.2, GROW0 + 1.1 + k * 0.2));
      b.bar.style.transform = `scaleX(${p.toFixed(4)})`;
      b.o.set(b.v * p);
      b.r.style.left = `${(112 + b.v * 14 * p + 24).toFixed(1)}px`;
      b.r.style.opacity = lin(t, GROW0 + k * 0.2, GROW0 + 0.2 + k * 0.2);
    });
    stampAnim(s.stamp, t, STAMP_F);
  } });

// ================================================================================================= 16 in one word
const WORDS = [["FOLD.", C.cobalt, "#fff", C.cobaltDk, "See chr22 in 3D", C.ink], ["TIME.", C.terra, "#fff", C.terraDk, "Watch it change in 4D", C.ink],
  ["DOSE.", C.ochre, C.ink, C.ochreDk, "Try a virtual epigenetic drug", C.ink], ["GENES.", C.paper, C.ink, "#C4C4BC", "Open or buried", C.cobalt],
  ["QUBITS.", C.ink, "#fff", "#000", "A simulated quantum chip", C.cobalt], ["PROOF.", C.violet, "#fff", C.violetDk, `${TS.n} tests · run once`, C.ink]];
scene({ t0: T.words, t1: T.outro + 0.5,
  build(s) {
    s.w = WORDS.map(([word, bg, fg, ext, line, pillBg], k) => {
      const g = mk("div", "fill", s.ui);
      mk("div", "fill", g, null, { background: bg });
      const deco = bg === C.paper ? mk("div", "fill halftone", g, null, { opacity: 0.6 }) : mk("div", "burst", g);
      if (bg === C.ink) deco.style.opacity = 0.5;
      const el = mk("div", "abs dx", g, null, { left: "0", width: "1920px", top: "300px", textAlign: "center", fontSize: "290px", color: fg, textShadow: extrude(ext, 15, 1.1, 1.4) });
      const ch = chars(el, word);
      let bar = null;
      if (bg === C.paper) bar = mk("div", "abs", g, null, { left: "560px", top: "590px", width: "800px", height: "22px", borderRadius: "11px", background: C.cobalt, transformOrigin: "0 50%" });
      const pill = mk("div", "pill", g, `<span style="font-size:1.1em">▲</span> ${line}`, { left: "0", top: "690px", background: pillBg, color: "#fff", fontSize: "30px", padding: "20px 36px" });
      const at = T.words + k;
      cue(at, "hit");
      shake(at + 0.04, 12);
      burst(at + 0.08, 960, 420, { n: 70, speed: 1100, life: 1.5, colors: bg === C.ochre ? [C.ink, "#fff", C.terra] : ["#fff", C.ochre, C.terra, C.cobaltHi], sound: false });
      return { g, ch, bar, pill, at };
    });
    flash(T.words, 0.5, 0.08);
  },
  update(s, t) {
    s.w.forEach((w, k) => {
      const on = t >= w.at && (k === s.w.length - 1 || t < w.at + 1);
      w.g.style.display = on ? "block" : "none";
      if (!on) return;
      const lt = t - w.at;
      w.ch.forEach((c, i) => {
        const p = lin(lt, i * 0.022, i * 0.022 + 0.17), q = E.o3(p);
        setT(c, `scale(${(mix(2.3, 1, q) * (1 + 0.04 * lt)).toFixed(4)}) rotate(${((1 - q) * -6).toFixed(2)}deg)`, clamp(p * 2.5));
      });
      if (w.bar) w.bar.style.transform = `scaleX(${E.o4(lin(lt, 0.15, 0.45)).toFixed(4)})`;
      w.pill.style.left = `${(960 - w.pill.offsetWidth / 2).toFixed(1)}px`;
      pop(w.pill, lt, 0.2, { s: 0.6, y: 30, d: 0.35 });
      const deco = w.g.children[1];
      if (deco.classList.contains("burst")) deco.style.transform = `rotate(${(t * 9).toFixed(2)}deg)`;
    });
  } });

// ================================================================================================= outro
scene({ t0: T.outro, t1: T.end + 0.1, inT: { type: "iris", at: T.outro, d: 0.45, sound: "boom" },
  build(s) {
    s.burst = bgBurst(s, C.cobalt);
    const u = s.ui;
    s.logo = mk("div", "abs logo", u, logoSvg(150, "#fff", C.ochre, 2.2), { left: "885px", top: "92px", width: "150px", height: "150px" });
    s.name = chars(mk("div", "abs dn", u, null, { left: "0", width: "1920px", top: "290px", textAlign: "center", fontSize: "180px", color: "#fff", textShadow: extrude(C.cobaltDk, 13) }), "ChronoCell-5D");
    s.tag = mk("div", "pill", u, "Chromatin 3D / 4D workstation", { left: "0", top: "510px", background: C.ink, color: "#fff", fontSize: "28px", padding: "18px 32px" });
    s.card = mk("div", "abs", u, null, { left: "310px", top: "650px", width: "1300px", height: "190px", borderRadius: "28px", background: C.ink, boxShadow: "0 30px 70px rgba(0,0,0,.35)", textAlign: "center", paddingTop: "34px" });
    mk("div", "dc", s.card, `${TS.n} tests · <span style="color:${C.pass}">${TS.pass} passed</span> · <span style="color:${C.fail}">${TS.fail} failed</span> · <span style="color:${C.other}">${TS.other} mixed</span>`, { fontSize: "66px", color: "#fff" });
    mk("div", "mono", s.card, `${TS.fixed} fails later passed as a retest · ${TS.open} still open · read from the app`, { font: "600 21px/1 var(--mono)", letterSpacing: ".12em", textTransform: "uppercase", color: "rgba(255,255,255,.7)", marginTop: "26px" });
    s.line = mk("div", "abs mono", u, "— A CHRONOCELL-5D MOTION REEL —", { left: "0", width: "1920px", textAlign: "center", top: "905px", font: "600 24px/1 var(--mono)", letterSpacing: ".3em", color: "rgba(255,255,255,.88)" });
    burst(T.outro + 0.35, 960, 400, { n: 190, speed: 1600, life: 2.8, colors: ["#fff", C.ochre, C.terra, "#9AA4FF", C.ink] });
    burst(T.outro + 2.45, 960, 745, { n: 90, speed: 1100, life: 2.2, colors: ["#fff", C.pass, C.ochre] });
    flash(T.outro + 0.3, 0.5, 0.12);
    shake(T.outro + 0.32, 16);
    cue(T.outro + 2.4, "hit");
  },
  update(s, t) {
    s.burst.style.transform = `rotate(${(t * 6).toFixed(2)}deg)`;
    pop(s.logo, t, T.outro + 0.25, { s: 0.2, r: -40, y: 0, d: 0.5 });
    animSpans(s.name, t, { in: T.outro + 0.35, st: 0.035, din: 0.55, fy: -300, fr: 10 });
    s.tag.style.left = `${(960 - s.tag.offsetWidth / 2).toFixed(1)}px`;
    pop(s.tag, t, T.outro + 1.3, { s: 0.6 });
    pop(s.card, t, T.outro + 2.4, { s: 0.6, y: 60, d: 0.55 });
    rise(s.line, t, T.outro + 3.2);
  } });

// ================================================================================================= HUD and overlays
const CHAP = [[T.open, "01", "The workstation", "chromatin 3D / 4D, in one app"], [T.fold, "02", "The fold", "human chr22, reconstructed in 3D"],
  [T.map, "03", "The contacts", "what Micro-C measures"], [T.title, null], [T.ws1, "04", "3D structure", "workspace 01"],
  [T.ws2, "05", "4D dynamics", "workspace 02"], [T.ws3, "06", "Compare", "workspace 03"], [T.ws4, "07", "Drug lab", "workspace 04"],
  [T.stack, "08", "Genes + Guide", "workspaces 05 · 08"], [T.quantum, "09", "Quantum lab", "workspace 06 · simulated"],
  [T.bento, "10", "Everything", "eight workspaces, one app"], [T.score, "11", "The scoreboard", "workspace 07 · every claim tested"],
  [T.loops, "12", "DNA loops", "gate 6b · new cell types"], [T.mol, "13", "Molecules", "gate Q6d · simulated chip"],
  [T.mit, "14", "Noisy chips", "gate Q9 · error mitigation"], [T.fail, "15", "The honest fail", "gate Q7b · quantum docking"],
  [T.words, "16", "In one word", "what each workspace is for"], [T.outro, null], [T.end + 1, null]];
CHAP.forEach((c) => {
  if (!c[1]) return;
  const b = mk("div", "badge", L.hud);
  mk("div", "n", b, c[1]);
  const tx = mk("div", "", b);
  mk("div", "t", tx, c[2]);
  mk("div", "s", tx, c[3]);
  c.el = b;
});
const tcEl = mk("div", "tc", L.hud);
function hud(t) {
  CHAP.forEach((c, i) => {
    if (!c.el) return;
    const next = CHAP[i + 1][0], on = t >= c[0] && t < next;
    c.el.style.display = on ? "flex" : "none";
    if (!on) return;
    const a = E.oBack(lin(t, c[0] + 0.3, c[0] + 0.75)), b = E.i3(lin(t, next - 0.3, next - 0.05));
    setT(c.el, `translate3d(${((1 - clamp(a)) * -40).toFixed(1)}px,${(b * 40).toFixed(1)}px,0) scale(${mix(0.85, 1, clamp(a, 0, 1.1)).toFixed(4)})`, clamp(a * 1.4) * (1 - b));
  });
  const f = Math.floor(t * FPS), ss = Math.floor(t), p2 = (n) => String(n).padStart(2, "0");
  tcEl.textContent = `CHRONOCELL-5D  ·  00:${p2(Math.floor(ss / 60))}:${p2(ss % 60)}:${p2(f % FPS)}`;
  tcEl.style.opacity = (1 - lin(t, T.end - 1.5, T.end - 0.8)).toFixed(3);
}
const CUR = mk("div", "cursor", L.cur, CURSOR_SVG);
function cursor(t) {
  const tr = TRACKS.find((c) => t >= c.from && t <= c.to);
  if (!tr) { CUR.style.opacity = 0; return; }
  const [x, y] = trackPos(tr, t);
  let sc = 1;
  for (const c of tr.clicks) { const d = t - c; if (d > -0.08 && d < 0.2) sc = Math.min(sc, 1 - 0.2 * Math.sin(Math.PI * clamp((d + 0.08) / 0.28))); }
  setT(CUR, `translate3d(${(x - 6).toFixed(1)}px,${(y - 4).toFixed(1)}px,0) scale(${sc.toFixed(4)})`, lin(t, tr.from, tr.from + 0.15) * (1 - lin(t, tr.to - 0.15, tr.to)));
}
function drawBurst(g, b, dt) {
  const r = rng(b.seed);
  for (let i = 0; i < b.n; i++) {
    const ang = b.angle + (r() - 0.5) * b.spread, sp = b.speed * (0.3 + 0.7 * r());
    const vx = Math.cos(ang) * sp, vy = Math.sin(ang) * sp;
    const k = 2.4, gk = (1 - Math.exp(-k * dt)) / k;
    const x = b.x + vx * gk, y = b.y + vy * gk + (950 * (dt - gk)) / k;
    const w = 10 + r() * 13, h = 6 + r() * 8, rot = r() * 6.283 + dt * (r() - 0.5) * 16, flip = Math.cos(dt * (5 + r() * 11) + r() * 6.283);
    const col = b.colors[Math.floor(r() * b.colors.length)], shape = r();
    g.globalAlpha = 1 - lin(dt, b.life * 0.55, b.life);
    g.fillStyle = col;
    g.save();
    g.translate(x, y);
    g.rotate(rot);
    g.scale(1, Math.abs(flip) < 0.15 ? 0.15 : flip);
    if (shape < 0.72) g.fillRect(-w / 2, -h / 2, w, h);
    else { g.beginPath(); g.arc(0, 0, h * 0.62, 0, 6.2832); g.fill(); }
    g.restore();
  }
  g.globalAlpha = 1;
}
const flashEl = mk("div", "flash", L.ov), blackEl = mk("div", "black", L.ov);
function overlays(t) {
  fxg.clearRect(0, 0, W, H);
  for (const b of BURSTS) { const dt = t - b.t; if (dt >= 0 && dt < b.life) drawBurst(fxg, b, dt); }
  for (const tr of TRACKS) for (const c of tr.clicks) {
    const d = t - c;
    if (d < 0 || d >= 0.6) continue;
    const [x, y] = trackPos(tr, c);
    [0, 0.09].forEach((o) => {
      const e = d - o;
      if (e < 0) return;
      fxg.strokeStyle = tr.ring;
      fxg.lineWidth = 4 * (1 - e / 0.6);
      fxg.globalAlpha = 1 - e / 0.51;
      fxg.beginPath(); fxg.arc(x, y, 12 + 80 * E.o3(clamp(e / 0.5)), 0, 6.2832); fxg.stroke();
    });
    fxg.globalAlpha = 1;
  }
  let fl = 0;
  for (const f of FLASHES) { const d = t - f.t; if (d >= 0 && d < f.d * 6) fl += f.a * Math.exp(-d / f.d); }
  flashEl.style.opacity = Math.min(1, fl).toFixed(3);
  let sx = 0, sy = 0;
  for (const s of SHAKES) {
    const d = t - s.t;
    if (d < 0 || d > 0.6) continue;
    const a = s.amp * Math.exp(-d * 8);
    sx += a * Math.sin(d * 95 + s.t * 7);
    sy += a * Math.cos(d * 77 + s.t * 3);
  }
  world.style.transform = Math.abs(sx) + Math.abs(sy) > 0.05 ? `translate3d(${sx.toFixed(2)}px,${sy.toFixed(2)}px,0) scale(1.02)` : "none";
  blackEl.style.opacity = Math.max(1 - lin(t, 0, 0.35), lin(t, T.end - 1.0, T.end - 0.05)).toFixed(3);
}

// ================================================================================================= seek + boot
function seek(t) {
  for (const s of SCENES) {
    const on = t >= s.t0 && t < s.t1;
    show(s.root, on);
    if (on) { trans(s, t); s.update(s, t, t - s.t0); }
  }
  hud(t);
  cursor(t);
  overlays(t);
  return t;
}
window.seek = seek;
window.TIMES = T;
window.DURATION = T.end;
window.CUES = CUES.sort((a, b) => a.t - b.t);
window.READY = (async () => {
  await document.fonts.ready;
  await Promise.all(["900 100px Archivo", "600 20px 'IBM Plex Mono'", "700 13px 'IBM Plex Mono'", "500 20px 'Inter Tight'"].map((f) => document.fonts.load(f)));
  const imgs = [...document.querySelectorAll("img")];
  await Promise.all(imgs.map((i) => (i.complete ? Promise.resolve() : new Promise((ok) => { i.onload = ok; i.onerror = ok; })).then(() => i.decode?.().catch(() => {}))));
  seek(0);
  return true;
})();

// preview: index.html?t=12.5 shows one frame; ?play plays in real time; ?bare hides everything but backgrounds and 3D
(async () => {
  const q = new URLSearchParams(location.search);
  await window.READY;
  if (q.has("bare")) { SCENES.forEach((s) => (s.ui.style.visibility = "hidden")); [L.fx, L.hud, L.cur, L.ov].forEach((e) => (e.style.display = "none")); }
  if (q.has("t")) seek(parseFloat(q.get("t")));
  if (q.has("play")) {
    const t0 = performance.now() - 1000 * parseFloat(q.get("play") || "0");
    const loop = () => { seek(((performance.now() - t0) / 1000) % T.end); requestAnimationFrame(loop); };
    loop();
  }
})();
