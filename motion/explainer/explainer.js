/* ChronoCell-5D explainer: what the app is, the app running live, and how its accuracy is tested.
 * Same machinery as the reel (../film.js): every frame is a pure function of time t; window.seek(t) draws it (and waits
 * for the screen-recording videos to land on the right frame). Screen recordings: clips/*.mp4 (record_app.py).
 * Terminal sessions: real output of the test suite and of recheck.py (record_terminal.py). Numbers: ../data/data.js. */
"use strict";

const T = { title: 0, prob: 9, tour: 41, clips: 44, test: 167, pytest: 187, recheck: 209, results: 252, outro: 274, end: 288 };
const W = 1920, H = 1080, FPS = 30;
const D = window.DATA, F = D.fold, TS = D.tests, HI = D.hi, EX = window.EX;
const C = { cobalt: "#3340D1", cobaltHi: "#5F6BFF", cobaltDk: "#1F2896", terra: "#E8582C", terraDk: "#A93A16", ochre: "#F5B931",
  ochreDk: "#B07F12", violet: "#6E45D6", paper: "#F2F2EF", paper2: "#E9E9E4", ink: "#1C1E1B", ink2: "#474A45", muted: "#62645F",
  night: "#111318", pass: "#27B26B", fail: "#F0552E", other: "#8D93A8" };

// ------------------------------------------------------------------------------------------------- helpers
const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
const lin = (t, a, b) => (b === a ? (t >= b ? 1 : 0) : clamp((t - a) / (b - a)));
const mix = (a, b, p) => a + (b - a) * p;
const E = {
  o3: (x) => 1 - Math.pow(1 - x, 3),
  i3: (x) => x * x * x,
  io3: (x) => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2),
  oExpo: (x) => (x >= 1 ? 1 : 1 - Math.pow(2, -10 * x)),
  oBack: (x) => { const c1 = 1.6, c3 = c1 + 1; return x <= 0 ? 0 : x >= 1 ? 1 : 1 + c3 * Math.pow(x - 1, 3) + c1 * Math.pow(x - 1, 2); },
};
const hexA = (h, a) => { const n = parseInt(h.slice(1), 16); return `rgba(${n >> 16},${(n >> 8) & 255},${n & 255},${a})`; };
const fmt = (n) => Math.round(n).toLocaleString("en-US");
function mk(tag, cls, parent, html, style) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  if (style) Object.assign(e.style, style);
  if (parent) parent.appendChild(e);
  return e;
}
function setT(el, tr, op) { el.style.transform = tr; if (op != null) el.style.opacity = clamp(op).toFixed(4); }
const show = (el, on, how = "block") => { el.style.display = on ? how : "none"; };
function extrude(col, n = 9, dx = 1, dy = 1.3, soft = "rgba(0,0,0,.25)") {
  const a = [];
  for (let i = 1; i <= n; i++) a.push(`${(i * dx).toFixed(1)}px ${(i * dy).toFixed(1)}px 0 ${col}`);
  a.push(`${n * dx + 8}px ${n * dy + 18}px 30px ${soft}`);
  return a.join(",");
}
function pop(el, t, at, o = {}) {
  const p = E.oBack(lin(t, at, at + (o.d ?? 0.45)));
  setT(el, `translate3d(${(o.x ?? 0) * (1 - p)}px,${(o.y ?? 30) * (1 - p)}px,0) scale(${mix(o.s ?? 0.7, 1, p).toFixed(4)})`, clamp(p * 2));
  return p;
}
const rise = (el, t, at, dy = 24, d = 0.45) => setT(el, `translate3d(0,${((1 - E.o3(lin(t, at, at + d))) * dy).toFixed(2)}px,0)`, lin(t, at, at + d * 0.7));
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
  return { el, set(v) {
    const N = Math.max(0, v) * Math.pow(10, dec);
    let prev = 0;
    for (const p of digits) {
      const q = p.place === 0 ? N % 10 : (Math.floor(N / Math.pow(10, p.place)) % 10) + (prev > 9 ? prev - 9 : 0);
      prev = q;
      p.strip.style.transform = `translateY(${(-q).toFixed(4)}em)`;
    }
    for (const p of parts) (p.col || p.sep).style.display = p.place <= dec || N >= Math.pow(10, p.place) * 0.999 ? "inline-block" : "none";
  } };
}
const logoSvg = (size, c1, c2, sw = 2) => `<svg width="${size}" height="${size}" viewBox="0 0 26 26"><path d="M5 20c3-9 5-13 8-13s5 4 8 13" stroke="${c1}" stroke-width="${sw}" fill="none" stroke-linecap="round"/><path d="M5 6c3 9 5 13 8 13s5-4 8-13" stroke="${c2}" stroke-width="${sw}" fill="none" stroke-linecap="round"/></svg>`;

// ------------------------------------------------------------------------------------------------- registries
const CUES = [];
const cue = (t, k, x = {}) => CUES.push({ t: +t.toFixed(3), k, ...x });
const CAPS = [];                                            // [from, to, html]
const cap = (a, b, html) => CAPS.push([a, b, html]);

// ------------------------------------------------------------------------------------------------- stage
const stage = document.getElementById("stage");
const L = { sc: mk("div", "layer persp", stage), hud: mk("div", "layer", stage), ov: mk("div", "layer", stage) };
const SCENES = [];
function scene(o) {
  const root = mk("div", "scene", L.sc);
  const s = { ...o, root, bg: mk("div", "fill", root), host: mk("div", "fill", root), ui: mk("div", "fill persp", root) };
  o.build(s);
  SCENES.push(s);
  if (o.inT) cue(o.inT.at, o.inT.sound || "whoosh");
  return s;
}
function trans(s, t) {
  let tf = "", clip = "none", op = 1;
  const i = s.inT, o = s.outT;
  if (i) {
    const d = i.d || 0.5, p = E.io3(lin(t, i.at, i.at + d));
    if (p < 1) {
      if (i.type === "iris") {
        const x = i.x ?? 960, y = i.y ?? 540, R = Math.hypot(Math.max(x, W - x), Math.max(y, H - y)) + 30;
        clip = `circle(${(p * R).toFixed(1)}px at ${x}px ${y}px)`;
      } else if (i.type === "push-up") tf += `translate3d(0,${((1 - p) * H).toFixed(1)}px,0) `;
      else if (i.type === "push-left") tf += `translate3d(${((1 - p) * W).toFixed(1)}px,0,0) `;
      else if (i.type === "fade") op *= p;
    }
  }
  if (o) {
    const d = o.d || 0.5, p = E.io3(lin(t, o.at, o.at + d));
    if (p > 0) {
      if (o.type === "push-up") tf += `translate3d(0,${(-p * H).toFixed(1)}px,0) `;
      else if (o.type === "push-left") tf += `translate3d(${(-p * W).toFixed(1)}px,0,0) `;
      else if (o.type === "fade") op *= 1 - p;
    }
  }
  s.root.style.transform = tf || "none";
  s.root.style.clipPath = clip;
  s.root.style.opacity = op.toFixed(4);
}
function bgDark(s, glow = C.cobalt, gx = 50, gy = 50) {
  mk("div", "fill", s.bg, null, { background: C.night });
  mk("div", "fill dots-dark", s.bg);
  mk("div", "fill", s.bg, null, { background: `radial-gradient(circle at ${gx}% ${gy}%, ${hexA(glow, 0.45)} 0%, ${hexA(glow, 0.16)} 28%, rgba(0,0,0,0) 62%)` });
}
function bgPaper(s, tone = C.paper) { mk("div", "fill", s.bg, null, { background: tone }); mk("div", "fill dots-light", s.bg); }
function bgBurst(s, color) {
  mk("div", "fill", s.bg, null, { background: color });
  const b = mk("div", "burst", s.bg);
  mk("div", "fill", s.bg, null, { background: "radial-gradient(circle at 50% 45%, rgba(255,255,255,.2) 0%, rgba(255,255,255,0) 42%, rgba(0,0,0,.22) 100%)" });
  return b;
}
function head(parent, text, x, y, size, color, ext, cls = "dc") {
  return mk("div", "abs " + cls, parent, text, { left: `${x}px`, top: `${y}px`, fontSize: `${size}px`, color, textShadow: extrude(ext, 8, 1, 1.3, "rgba(0,0,0,.18)") });
}

// ================================================================================================= title
scene({ t0: T.title, t1: T.prob + 0.55, outT: { type: "push-up", at: T.prob },
  build(s) {
    s.burst = bgBurst(s, C.cobalt);
    const u = s.ui;
    s.logo = mk("div", "abs logo", u, logoSvg(170, "#fff", C.ochre, 2.2), { left: "875px", top: "150px" });
    s.name = mk("div", "abs dn", u, "ChronoCell-5D", { left: "0", width: "1920px", top: "380px", textAlign: "center", fontSize: "180px", color: "#fff", textShadow: extrude(C.cobaltDk, 12) });
    s.sub = mk("div", "abs", u, "What it is · the app, live · how it is tested", { left: "0", width: "1920px", top: "610px", textAlign: "center", font: "600 40px/1.2 var(--sans)", color: "#fff" });
    s.pill = mk("div", "pill", u, "A five-minute tour", { left: "780px", top: "700px", background: C.ink, color: "#fff", fontSize: "26px" });
    cue(0.35, "boom");
  },
  update(s, t) {
    s.burst.style.transform = `rotate(${(t * 5).toFixed(2)}deg)`;
    pop(s.logo, t, 0.3, { s: 0.2, y: 0, d: 0.6 });
    pop(s.name, t, 0.6, { s: 0.8, y: 50, d: 0.6 });
    rise(s.sub, t, 1.2);
    s.pill.style.left = `${(960 - s.pill.offsetWidth / 2).toFixed(1)}px`;
    pop(s.pill, t, 1.6, { s: 0.6 });
  } });

// ================================================================================================= the problem
const MAP = F.map, NB = MAP.bins;
const CMAP = (() => {
  const st = [[0, [248, 248, 246]], [0.25, [235, 201, 181]], [0.55, [217, 116, 74]], [0.8, [168, 56, 15]], [1, [74, 23, 5]]];
  return Array.from({ length: 256 }, (_, k) => {
    const u = k / 255;
    let i = 0;
    while (i < st.length - 2 && u > st[i + 1][0]) i++;
    const f = clamp((u - st[i][0]) / (st[i + 1][0] - st[i][0]));
    const c = st[i][1].map((v, n) => Math.round(mix(v, st[i + 1][1][n], f)));
    return `rgb(${c[0]},${c[1]},${c[2]})`;
  });
})();
const P = { a: T.prob, b: T.prob + 8, c: T.prob + 18, d: T.prob + 26, e: T.tour };
cap(P.a + 0.6, P.b - 0.1, "Every human cell holds about <b>two metres of DNA</b>, packed into a nucleus a hundred times thinner than a hair.");
cap(P.b, P.b + 4.6, "<b>How</b> it is folded decides which genes can be read.");
cap(P.b + 4.6, P.c - 0.1, "Experiments such as Hi-C and Micro-C can’t photograph the fold. They measure <b>which pieces of DNA touch</b>.");
cap(P.c, P.d - 0.1, "ChronoCell-5D turns those contacts into a <b>3D model</b>, and says how sure it is about each distance.");
cap(P.d, P.e - 0.1, "3D: the shape. 4th dimension: change over time. 5th: <b>what the changes mean</b>.");
scene({ t0: T.prob, t1: T.tour + 0.55, inT: { type: "push-up", at: T.prob }, outT: { type: "push-left", at: T.tour },
  build(s) {
    bgDark(s, C.cobalt, 66, 46);
    const u = s.ui;
    s.k1 = kick(u, "Inside every human cell", 120, 200, "#8E98FF");
    s.big = mk("div", "abs dc", u, "2 METRES", { left: "110px", top: "240px", fontSize: "190px", color: "#fff", textShadow: extrude(C.cobaltDk, 11) });
    s.small = mk("div", "abs lead", u, "of DNA, folded into a nucleus about 10 micrometres across", { left: "120px", top: "430px", width: "640px", color: "rgba(255,255,255,.82)", fontSize: "30px" });
    // contact map card
    s.card = mk("div", "abs", u, null, { left: "1000px", top: "130px", width: "760px", height: "760px", borderRadius: "26px", background: "#fff", boxShadow: "0 30px 80px rgba(0,0,0,.4)" });
    s.cv = mk("canvas", "", s.card, null, { position: "absolute", left: "30px", top: "30px", width: "700px", height: "700px" });
    s.cv.width = 700; s.cv.height = 700;
    s.g = s.cv.getContext("2d");
    s.px = [];
    let seed = 99;
    const r = () => { seed = (seed * 16807) % 2147483647; return seed / 2147483647; };
    for (let i = 0; i < NB; i++) for (let j = i; j < NB; j++) {
      const lev = MAP.level[i * NB + j];
      if (lev > 2) s.px.push([i, j, lev, P.b + 4.9 + ((j - i) / (NB - 1)) * 1.6 + r() * 0.3, r() * 6.283, 300 + r() * 500]);
    }
    s.k2 = kick(u, "Simulated Micro-C · chr22", 120, 200, "#FF9A72");
    s.cnt = odo(u, fmt(F.contacts), { position: "absolute", left: "110px", top: "240px", fontSize: "170px", color: "#fff", textShadow: extrude(C.terraDk, 11) });
    s.cntl = mk("div", "abs lead", u, "pairs of DNA pieces seen touching, in the app’s reference model of chr22", { left: "120px", top: "420px", width: "640px", color: "rgba(255,255,255,.82)", fontSize: "30px" });
    s.arrow = mk("div", "abs", u, `<svg width="160" height="60" viewBox="0 0 160 60"><path d="M6 30 H140 M118 10 L144 30 L118 50" fill="none" stroke="#fff" stroke-width="8" stroke-linecap="round" stroke-linejoin="round"/></svg>`, { left: "590px", top: "430px" });
    s.cards = [["3D", "the shape", "Rebuilds the fold of a whole chromosome.", C.cobalt], ["4D", "over time", "Watches it change: time courses, drugs, DNA rearrangements.", C.terra],
      ["5D", "what it means", "Explains the changes: genes, loops, domains, plain words.", C.violet]].map(([n, a, b, col], k) => {
      const e = mk("div", "abs", u, null, { left: `${150 + k * 560}px`, top: "250px", width: "520px", height: "460px", borderRadius: "30px", background: col, padding: "40px 40px", color: "#fff", boxShadow: "0 30px 70px rgba(0,0,0,.35)" });
      mk("div", "dn", e, n, { fontSize: "150px", textShadow: extrude("rgba(0,0,0,.25)", 8) });
      mk("div", "dc", e, a, { fontSize: "64px", marginTop: "20px" });
      mk("div", "lead", e, b, { fontSize: "27px", marginTop: "16px", color: "rgba(255,255,255,.88)" });
      return e;
    });
    cue(P.a + 0.6, "rise", { d: 4.4 });
    cue(P.b + 4.9, "roll", { d: 2.0 });
    cue(P.c, "whoosh");
    cue(P.d, "whoosh");
    [0, 1, 2].forEach((k) => cue(P.d + 0.3 + k * 0.35, "pop"));
  },
  update(s, t) {
    // a. the fold draws itself (the reel's chapter timeline, started 0.4 s in)
    const foldOn = t < P.b + 4.6 || (t >= P.c && t < P.d);
    s.host.style.opacity = t < P.b + 4.6 ? 1 - lin(t, P.b + 4.0, P.b + 4.6) : lin(t, P.c, P.c + 0.6) * (1 - lin(t, P.d - 0.5, P.d));
    if (foldOn) {
      if (t < P.c) GL.drawFold(t - P.a + 4.6, "chapter", s.host);
      else GL.drawFold(t - P.c, "free", s.host, { sx: 1380, sy: 500, r: 6.6, az: 0.25 * (t - P.c) - 0.4 });
    }
    const A = 1 - lin(t, P.b + 1.3, P.b + 1.9);              // clear the text as the camera flies into the fold
    rise(s.k1, t, P.a + 0.4); s.k1.style.opacity *= A;
    pop(s.big, t, P.a + 0.7, { s: 0.6 }); s.big.style.opacity *= A;
    rise(s.small, t, P.a + 1.4); s.small.style.opacity *= A;
    // b. the contact map builds itself
    const g = s.g, px = 700 / NB;
    g.clearRect(0, 0, 700, 700);
    if (t >= P.b + 4.6 && t < P.d) {
      for (const [i, j, lev, d0, a, rad] of s.px) {
        const dt = t - d0;
        if (dt < 0) continue;
        const q = E.o3(clamp(dt / 0.55)), sz = px * mix(0.4, 1, q) + 0.6;
        g.fillStyle = CMAP[lev];
        g.fillRect(mix(350 + Math.cos(a) * rad, j * px, q), mix(350 + Math.sin(a) * rad, i * px, q), sz, sz);
        g.fillRect(mix(350 + Math.sin(a) * rad, i * px, q), mix(350 + Math.cos(a) * rad, j * px, q), sz, sz);
      }
    }
    const cardIn = E.oBack(lin(t, P.b + 4.4, P.b + 4.9)), toLeft = E.io3(lin(t, P.c, P.c + 0.8)), cardOut = lin(t, P.d - 0.5, P.d);
    setT(s.card, `translate3d(${(-850 * toLeft).toFixed(1)}px,${(80 * toLeft).toFixed(1)}px,0) scale(${(mix(0.8, 1, clamp(cardIn)) * mix(1, 0.62, toLeft)).toFixed(4)})`, clamp(cardIn * 2) * (1 - cardOut));
    s.card.style.transformOrigin = "0 0";
    const B = lin(t, P.b + 4.6, P.b + 5.0) * (1 - lin(t, P.c - 0.4, P.c));
    s.k2.style.opacity = B; s.cnt.el.style.opacity = B; s.cntl.style.opacity = B;
    s.cnt.set(F.contacts * E.o3(lin(t, P.b + 4.9, P.b + 7.0)));
    pop(s.arrow, t, P.c + 0.7, { s: 0.3, x: -40, y: 0 });
    s.arrow.style.opacity = Math.min(+s.arrow.style.opacity, 1 - cardOut);
    s.arrow.style.left = "640px"; s.arrow.style.top = "390px";
    s.cards.forEach((e, k) => { pop(e, t, P.d + 0.3 + k * 0.35, { s: 0.5, y: 80, d: 0.55 }); });
  } });
function kick(parent, text, x, y, color) { return mk("div", "abs kick", parent, text, { left: `${x}px`, top: `${y}px`, color }); }

// ================================================================================================= the app, live
scene({ t0: T.tour, t1: T.clips + 0.55, inT: { type: "push-left", at: T.tour }, outT: { type: "push-left", at: T.clips },
  build(s) {
    bgPaper(s);
    const u = s.ui;
    s.live = mk("div", "live", u, `<span class="dot"></span>LIVE SCREEN RECORDINGS`, { left: "120px", top: "330px" });
    s.h = head(u, "THE APP, LIVE", 110, 390, 200, C.ink, C.cobalt);
    s.l = mk("div", "abs lead", u, "The real app, running on this computer: every click, every 3D view and every number on screen is computed as you watch.", { left: "120px", top: "610px", width: "1300px", color: C.ink2, fontSize: "34px" });
  },
  update(s, t) {
    pop(s.live, t, T.tour + 0.4, { s: 0.6 });
    pop(s.h, t, T.tour + 0.5, { s: 0.8, y: 60 });
    rise(s.l, t, T.tour + 0.9);
  } });

// ---- one scene per recording
const SW = { x: 260, y: 96, w: 1400, h: 832 };               // the browser window (video area 1400 x 788)
const VW = SW.w, VH = SW.h - 44;
const CLIPS = [
  { c: "c01_launch", d: 8, n: "01", label: "Eight workspaces", url: "the app’s home",
    caps: [[0, 8, "One app, eight workspaces: <b>3D, 4D, compare, drug lab, genes, guide, quantum lab, scoreboard</b>."]] },
  { c: "c02_structure", d: 18, n: "02", label: "3D structure", url: "01 3D structure",
    zoom: [[0, 1, 0.5, 0.5], [3, 1.28, 0.32, 0.52], [14, 1.28, 0.32, 0.52], [16, 1, 0.5, 0.5]],
    caps: [[0, 6, "<b>3D structure</b>: the fold of a whole chromosome, here human chr22 in 5,082 beads of 10 kb."],
      [6, 12, "Drag to turn it. Colour runs along the DNA, from one end to the other."], [12, 18, "<b>Turntable</b> spins it; the panel measures it: size, scaling, overlaps."]] },
  { c: "c03_dynamics", d: 12, n: "03", label: "4D dynamics", url: "02 4D dynamics",
    caps: [[0, 6, "<b>4D dynamics</b>: watch the fold change over time."], [6, 12, "Here a DNA deletion (22q11.2, 2.1 Mb): press Play and the fold relaxes, frame by frame."]] },
  { c: "c04_compare", d: 12, n: "04", label: "Compare", url: "03 Compare",
    caps: [[0, 6, "<b>Compare</b>: two biological states side by side."], [6, 12, "Rotate one view, and the other follows."]] },
  { c: "c05_druglab", d: 11, n: "05", label: "Drug lab", url: "04 Drug lab",
    caps: [[0, 11, "<b>Drug lab</b>: apply a virtual epigenetic drug, then play the dose up and watch the fold respond."]] },
  { c: "c06_genes", d: 12, n: "06", label: "Genes", url: "05 Genes",
    caps: [[0, 12, "<b>Genes</b>: which genes sit in open, active chromatin, and which are buried and likely silenced."]] },
  { c: "c07_guide", d: 8, to: 8, n: "07", label: "Guide", url: "06 Guide",
    caps: [[0, 8, "<b>Guide</b>: what everything means, in plain words, with pictures."]] },
  { c: "c08_quantum", d: 24, n: "08", label: "Quantum lab", url: "07 Quantum lab",
    zoom: [[0, 1, 0.5, 0.5], [5.6, 1, 0.5, 0.5], [6.4, 1.25, 0.62, 0.86], [12.8, 1.25, 0.62, 0.86], [13.6, 1, 0.5, 0.5]],
    marks: [[6.6, 13.0, 0.03, 0.83, 0.97, 0.99, "VQE energy vs the exact answer"], [19.6, 24.7, 0.02, 0.73, 0.45, 0.78, "the same circuit, with hardware noise"]],
    caps: [[0, 5.6, "<b>Quantum lab</b>: ChronoCell’s problems on a simulated quantum computer."],
      [5.6, 13, "Live: the energy of a HeH⁺ molecule by VQE, checked against the exact answer: <b>error 0.000 mHa</b> (limit 1.6)."],
      [13, 18.6, "Now add realistic hardware noise and run it again…"],
      [18.6, 24, "…the same circuit is off by <b>763 mHa</b>. That is why error mitigation matters."]] },
  { c: "c09_scoreboard", d: 18, n: "09", label: "Scoreboard", url: "08 Scoreboard",
    caps: [[0, 7, "<b>Scoreboard</b>: every accuracy test on held-out real data."], [7, 18, "What passed, what failed, and by how much, read straight from the result files."]] },
];
{
  let t0 = T.clips;
  CLIPS.forEach((k, i) => {
    k.t0 = t0; t0 += k.d; k.t1 = t0; k.meta = EX.clips[k.c];
    k.from = k.from ?? 0; k.to = Math.min(k.to ?? k.meta.secs, k.meta.secs); k.rate = (k.to - k.from) / k.d;
  });
}
CLIPS.forEach((k, i) => {
  k.caps.forEach(([a, b, h]) => cap(k.t0 + a, k.t0 + b - 0.05, h));
  k.meta.clicks.forEach((c) => { if (c.t >= k.from && c.t <= k.to) cue(k.t0 + (c.t - k.from) / k.rate, "click"); });
  scene({ t0: k.t0, t1: k.t1 + 0.55, inT: { type: "push-left", at: k.t0 }, outT: { type: "push-left", at: k.t1 }, clip: k,
    build(s) {
      bgPaper(s, i % 2 ? C.paper : C.paper2);
      const u = s.ui;
      s.win = mk("div", "screen", u, null, { left: `${SW.x}px`, top: `${SW.y}px`, width: `${SW.w}px`, height: `${SW.h}px` });
      const bar = mk("div", "bar", s.win, "<i></i><i></i><i></i>");
      mk("div", "url", bar, `localhost:8501 · ChronoCell-5D · ${k.url}`);
      const vp = mk("div", "vp", s.win);
      s.v = mk("video", "", vp);
      Object.assign(s.v, { src: `clips/${k.c}.mp4`, muted: true, preload: "auto", playsInline: true });
      s.marks = (k.marks || []).map((m) => {
        const e = mk("div", "abs", vp, null, { border: `4px solid ${C.cobalt}`, borderRadius: "10px", boxShadow: `0 0 0 6px ${hexA(C.cobalt, 0.18)}` });
        const l = mk("div", "abs", e, m[6], { left: "-4px", top: "-38px", padding: "6px 12px", borderRadius: "8px", background: C.cobalt, color: "#fff", font: "600 17px/1 var(--mono)", whiteSpace: "nowrap" });
        return { e, m, l };
      });
      s.badge = mk("div", "live", u, `<span class="dot"></span>LIVE${k.rate > 1.1 ? ` · ${k.rate.toFixed(1)}× SPEED` : ""}`, { left: `${SW.x + SW.w - 260}px`, top: `${SW.y + 58}px` });
      s.badge.style.left = "";
      s.badge.style.right = `${W - SW.x - SW.w + 16}px`;
    },
    async update(s, t) {
      const lt = t - k.t0;
      const vt = clamp(k.from + lt * k.rate, 0, k.meta.secs - 0.04);
      // zoom: keyframes [time in the slot, scale, focus x, focus y]
      let z = 1, fx = 0.5, fy = 0.5;
      const Z = k.zoom || [[0, 1, 0.5, 0.5]];
      for (let j = 0; j < Z.length; j++) {
        if (lt >= Z[j][0]) { [, z, fx, fy] = Z[j]; if (j + 1 < Z.length) { const p = E.io3(lin(lt, Z[j][0], Z[j + 1][0])); z = mix(Z[j][1], Z[j + 1][1], p); fx = mix(Z[j][2], Z[j + 1][2], p); fy = mix(Z[j][3], Z[j + 1][3], p); } }
      }
      z *= 1 + 0.025 * lin(lt, 0, k.d);                      // a slow push-in
      const tx = clamp(VW / 2 - fx * VW * z, VW - VW * z, 0), ty = clamp(VH / 2 - fy * VH * z, VH - VH * z, 0);
      s.v.style.transform = `translate3d(${tx.toFixed(1)}px,${ty.toFixed(1)}px,0) scale(${z.toFixed(4)})`;
      s.marks.forEach(({ e, m }) => {
        const on = vt >= m[0] && vt <= m[1];
        e.style.display = on ? "block" : "none";
        if (!on) return;
        const x0 = tx + m[2] * VW * z, y0 = ty + m[3] * VH * z, x1 = tx + m[4] * VW * z, y1 = ty + m[5] * VH * z;
        Object.assign(e.style, { left: `${x0 - 8}px`, top: `${y0 - 6}px`, width: `${x1 - x0 + 16}px`, height: `${y1 - y0 + 12}px` });
        e.style.opacity = lin(vt, m[0], m[0] + 0.3);
      });
      pop(s.badge, t, k.t0 + 0.5, { s: 0.6, y: -10 });
      await seekVideo(s.v, vt);
    } });
});
async function seekVideo(v, t) {
  if (v.readyState >= 2 && Math.abs(v.currentTime - t) < 0.004) return;
  await new Promise((res) => {
    const done = () => { v.removeEventListener("seeked", done); res(); };
    v.addEventListener("seeked", done);
    v.currentTime = t;
    setTimeout(done, 3000);
  });
  if (v.requestVideoFrameCallback) await new Promise((res) => { v.requestVideoFrameCallback(() => res()); setTimeout(res, 150); });
}

// ================================================================================================= how we know it works
const STEP = [["Practice data", "Tune the method on data set aside for practice.", "db"], ["Freeze the rule", "Write the pass / fail rule down (validation/frozen.py) before the test data are read.", "lock"],
  ["Run once", "Test on held-out real data, one time.", "play"], ["Keep the result", "Pass or fail, it goes on the Scoreboard. A retest needs a new method and new data.", "flag"]];
const ICON = {
  db: `<svg width="70" height="70" viewBox="0 0 70 70"><ellipse cx="35" cy="16" rx="24" ry="9" fill="none" stroke="${C.cobalt}" stroke-width="5"/><path d="M11 16v38c0 5 11 9 24 9s24-4 24-9V16M11 35c0 5 11 9 24 9s24-4 24-9" fill="none" stroke="${C.cobalt}" stroke-width="5"/></svg>`,
  lock: `<svg width="70" height="70" viewBox="0 0 70 70"><rect x="12" y="30" width="46" height="34" rx="7" fill="${C.ochre}"/><path d="M22 30V21a13 13 0 0 1 26 0v9" fill="none" stroke="${C.ink}" stroke-width="6"/></svg>`,
  play: `<svg width="70" height="70" viewBox="0 0 70 70"><circle cx="35" cy="35" r="30" fill="${C.terra}"/><path d="M28 21 L50 35 L28 49 Z" fill="#fff"/></svg>`,
  flag: `<svg width="70" height="70" viewBox="0 0 70 70"><path d="M16 64V8" stroke="${C.ink}" stroke-width="6" stroke-linecap="round"/><path d="M18 10h38l-9 13 9 13H18z" fill="${C.pass}"/></svg>`,
};
cap(T.test + 0.4, T.test + 7, "How do we know it works? <b>Every accuracy claim is a test</b>, with its rule written down first.");
cap(T.test + 7, T.test + 14, "Each test runs <b>once</b>, on data the method has never seen, so nobody can tune it to the answer.");
cap(T.test + 14, T.pytest - 0.05, "Failures stay on the record. So far: <b>37 tests · 14 passed · 16 failed · 7 mixed or blocked</b>.");
scene({ t0: T.test, t1: T.pytest + 0.55, inT: { type: "push-left", at: T.test }, outT: { type: "push-up", at: T.pytest },
  build(s) {
    bgPaper(s);
    const u = s.ui;
    s.h = head(u, "HOW WE KNOW IT WORKS", 110, 110, 120, C.ink, C.terra);
    s.cards = STEP.map(([tt, tx, ic], k) => {
      const e = mk("div", "step", u, null, { left: `${110 + k * 445}px`, top: "330px" });
      mk("div", "num", e, `STEP ${k + 1}`);
      mk("div", "ic", e, ICON[ic]);
      mk("div", "tt", e, tt);
      mk("div", "tx", e, tx);
      return e;
    });
    s.arrows = [0, 1, 2].map((k) => mk("div", "abs", u, `<svg width="40" height="40" viewBox="0 0 40 40"><path d="M6 20 H30 M20 9 L32 20 L20 31" fill="none" stroke="${C.ink}" stroke-width="5" stroke-linecap="round" stroke-linejoin="round"/></svg>`, { left: `${512 + k * 445}px`, top: "475px" }));
    const row = mk("div", "abs", u, null, { left: "110px", top: "720px", display: "flex", gap: "60px", alignItems: "baseline" });
    s.nums = [["37", "tests", C.ink], [String(TS.pass), "passed", C.pass], [String(TS.fail), "failed, kept", C.fail], [String(TS.other), "mixed or blocked", C.other]].map(([v, l, col]) => {
      const e = mk("div", "", row, null, { display: "flex", alignItems: "baseline", gap: "14px" });
      const o = odo(e, v, { fontSize: "110px", color: col });
      mk("div", "kick", e, l, { color: col });
      return { e, o, v: +v };
    });
    s.row = row;
    STEP.forEach((_, k) => cue(T.test + 1.0 + k * 1.4, "pop"));
    cue(T.test + 14, "roll", { d: 1.2 });
  },
  update(s, t) {
    pop(s.h, t, T.test + 0.4, { s: 0.8 });
    s.cards.forEach((e, k) => pop(e, t, T.test + 1.0 + k * 1.4, { s: 0.6, y: 60, d: 0.5 }));
    s.arrows.forEach((e, k) => pop(e, t, T.test + 1.7 + k * 1.4, { s: 0.3, x: -20, y: 0 }));
    s.row.style.opacity = lin(t, T.test + 13.8, T.test + 14.2);
    s.nums.forEach((n) => n.o.set(n.v * E.o3(lin(t, T.test + 14, T.test + 15.4))));
  } });

// ---- terminal replay
function terminal(parent, x, y, w, h, title) {
  const e = mk("div", "term", parent, null, { left: `${x}px`, top: `${y}px`, width: `${w}px`, height: `${h}px` });
  const bar = mk("div", "bar", e, "<i></i><i></i><i></i>");
  mk("div", "ttl", bar, title);
  const body = mk("div", "body", e);
  const rows = Math.floor((h - 44 - 36) / 26);
  const lines = Array.from({ length: rows }, () => mk("div", "ln", body));
  return { e, lines, rows };
}
const esc = (s) => s.replace(/&/g, "&amp;").replace(/</g, "&lt;").replace(/>/g, "&gt;");
function paint(s) {
  let h = esc(s);
  if (/^={2,}|^== /.test(s)) return `<span class="hd">${h}</span>`;
  h = h.replace(/\bPASSED\b/g, '<span class="ok">PASSED</span>').replace(/\bFAILED\b/g, '<span class="bad">FAILED</span>')
    .replace(/\bSKIPPED\b/g, '<span class="warn">SKIPPED</span>').replace(/\b(\d+ passed)\b/g, '<span class="ok">$1</span>')
    .replace(/\b(\d+ skipped)\b/g, '<span class="warn">$1</span>').replace(/(-&gt; PASS|  PASS |\bok\b|within\b|beats both)/g, '<span class="ok">$1</span>')
    .replace(/(-&gt; FAIL|  FAIL |MISMATCH|OUTSIDE)/g, '<span class="bad">$1</span>').replace(/  (OTHER) /g, '  <span class="warn">$1</span> ');
  return h;
}
function replay(term, data, lt, speed, typeDur = 1.4, hl = null) {
  const cmd = data.cmd, typed = Math.floor(clamp(lt / typeDur) * cmd.length);
  const tt = (lt - typeDur - 0.2) * speed;                  // time in the recorded session
  const out = [`<span class="pr">$</span> ${esc(cmd.slice(0, typed))}${lt < typeDur + 0.2 ? '<span class="caret"></span>' : ""}`];
  let n = 0;
  if (tt >= 0) for (const [t0] of data.lines) { if (t0 <= tt) n++; else break; }
  const shown = data.lines.slice(0, n).map(([, s]) => s);
  const all = out.concat(shown.map((s) => (hl && hl.test(s) ? `<span class="hl">${paint(s)}</span>` : paint(s))));
  const view = all.slice(Math.max(0, all.length - term.rows));
  term.lines.forEach((el, i) => { const h = view[i] ?? ""; if (el._h !== h) { el.innerHTML = h; el._h = h; } });
  return n;
}

// ---- the test suite
const PY = EX.terms.pytest, PY_SPEED = PY.secs / 17.5;
const PY_SUM = PY.lines.map(([, s]) => s).reverse().find((s) => /passed/.test(s)) || "";
cap(T.pytest + 0.3, T.pytest + 9, "First, the software itself: the project’s <b>308 automated checks</b>, run for real on this computer.");
cap(T.pytest + 9, T.recheck - 0.05, `Result: <b>${(PY_SUM.match(/\d+ passed[^=]*/) || [""])[0].trim()}</b>, shown ${Math.round(PY_SPEED)}× faster than it ran.`);
scene({ t0: T.pytest, t1: T.recheck + 0.55, inT: { type: "push-up", at: T.pytest }, outT: { type: "push-left", at: T.recheck },
  build(s) {
    bgDark(s, C.pass, 50, 40);
    s.term = terminal(s.ui, 110, 80, 1700, 820, "Terminal — ChronoCell-5D — tests (real run, sped up)");
    s.badge = mk("div", "abs dc", s.ui, (PY_SUM.match(/\d+ passed/) || ["passed"])[0].toUpperCase(), { left: "0", width: "1920px", textAlign: "center", top: "390px", fontSize: "150px", color: "#fff", textShadow: extrude("#16673D", 10) });
    cue(T.pytest + 1.6, "roll", { d: 17 });
    cue(T.pytest + 1.6 + PY.secs / PY_SPEED, "hit");
  },
  update(s, t) {
    replay(s.term, PY, t - T.pytest, PY_SPEED, 1.4, /passed|PASSED.*100%/);
    const end = T.pytest + 1.6 + PY.secs / PY_SPEED;
    pop(s.badge, t, end + 0.1, { s: 0.4, y: 0 });
    s.term.e.style.filter = t > end ? `blur(${(4 * lin(t, end, end + 0.4)).toFixed(2)}px)` : "none";
  } });

// ---- the re-check
const RC = EX.terms.recheck, RC_SPEED = RC.secs / 34.0;
const SIDE = [
  [/^== Gate 6b/, "Gate 6b · DNA loops", `F1 ${HI.loops.hmec.chronocell.toFixed(2)} · ${HI.loops.hap1.chronocell.toFixed(2)}`, "beats both published tools on two new cell types", true, /verdict: PASS/],
  [/^== Gate Q6d/, "Gate Q6d · molecules", `${HI.q6d.within} / ${HI.q6d.cases} within`, `chemical accuracy (1.6 mHa); worst ${HI.q6d.worst} mHa`, true, /14 of 14 within/],
  [/^== Gate Q9/, "Gate Q9 · noisy chip", `${HI.q9.reduction}× smaller error`, `${HI.q9.within} / ${HI.q9.cases} molecules rescued by error mitigation`, true, /16 of 16 within/],
  [/^== Gate Q7b/, "Gate Q7b · docking", `${HI.q7b.qaoa} % vs ${HI.q7b.random} %`, "quantum route vs plain random search: a fail, kept", false, /quantum route \d+\/\d+/],
  [/^== Gate Q8/, "Gate Q8 · drug safety", `${HI.q8.passed} / ${HI.q8.of} endpoints`, "17 were needed: a fail, kept", false, /endpoints met/],
];
cap(T.recheck + 0.3, T.recheck + 10, "Then the accuracy. Each test ran once, so we don’t re-run it: <b>we re-check its arithmetic</b>.");
cap(T.recheck + 10, T.recheck + 35.5, "The script recomputes every score from the saved numbers and <b>re-applies the frozen rule</b>.");
cap(T.recheck + 35.5, T.results - 0.05, "Every recomputed verdict <b>matches the Scoreboard</b>, the passes and the fails.");
scene({ t0: T.recheck, t1: T.results + 0.55, inT: { type: "push-left", at: T.recheck }, outT: { type: "push-up", at: T.results },
  build(s) {
    bgDark(s, C.cobalt, 75, 40);
    s.term = terminal(s.ui, 60, 80, 1240, 820, "Terminal — ChronoCell-5D — re-check (real output, paced for reading)");
    s.when = SIDE.map(([sec, , , , , trig]) => {
      const h = RC.lines.findIndex(([, l]) => sec.test(l));
      const i = h < 0 ? -1 : RC.lines.findIndex(([, l], n) => n > h && trig.test(l));
      return i < 0 ? 1e9 : RC.lines[i][0];
    });
    s.side = SIDE.map(([, g, v, sub, ok], k) => {
      const e = mk("div", "side", s.ui, null, { left: "1340px", top: `${80 + k * 152}px` });
      mk("div", "g", e, g);
      mk("div", "v", e, v, { color: ok ? C.ink : C.fail });
      mk("div", "s", e, sub);
      mk("div", "vd", e, ok ? "PASS" : "FAIL", { background: ok ? C.pass : C.fail });
      return e;
    });
    s.all = mk("div", "stamp", s.ui, "ALL VERDICTS MATCH", { left: "360px", top: "430px", color: "#fff", background: C.pass, borderColor: "#fff", fontSize: "72px" });
    s.when.forEach((w) => cue(T.recheck + 1.6 + w / RC_SPEED, "pop"));
    s.done = T.recheck + 1.6 + RC.secs / RC_SPEED;
    cue(s.done + 0.2, "stamp");
  },
  update(s, t) {
    replay(s.term, RC, t - T.recheck, RC_SPEED, 1.4, /^ {3}\S.*(F1 0\.\d+|\d+ of \d+ within|quantum route \d+|endpoints met)|^Recomputed scores/);
    s.side.forEach((e, k) => pop(e, t, T.recheck + 1.6 + s.when[k] / RC_SPEED, { s: 0.6, x: 60, y: 0 }));
    const p = lin(t, s.done + 0.2, s.done + 0.42);
    setT(s.all, `rotate(${(-5 + (1 - E.o3(p)) * 8).toFixed(2)}deg) scale(${mix(2.4, 1, E.o3(p)).toFixed(4)})`, clamp(p * 3));
  } });

// ================================================================================================= results
cap(T.results + 0.3, T.results + 8, "Where ChronoCell-5D stands today, test by test, on the Scoreboard.");
cap(T.results + 8, T.outro - 0.05, "Strong on DNA loops and simulated quantum chemistry; honest about docking, drug safety and calibration.");
scene({ t0: T.results, t1: T.outro + 0.55, inT: { type: "push-up", at: T.results }, outT: { type: "fade", at: T.outro },
  build(s) {
    bgBurst(s, C.cobalt);
    const u = s.ui;
    s.h = head(u, "THE SCOREBOARD", 110, 90, 110, "#fff", C.cobaltDk);
    s.big = [[TS.pass, "passed", C.pass], [TS.fail, "failed, kept on the record", C.fail], [TS.other, "mixed, blocked or baseline", "#C9CCDA"]].map(([n, l, col], k) => {
      const e = mk("div", "abs", u, null, { left: "110px", top: `${270 + k * 190}px`, display: "flex", alignItems: "center", gap: "26px" });
      const o = odo(e, String(n), { fontSize: "150px", color: "#fff", textShadow: extrude("rgba(0,0,0,.25)", 8) });
      mk("div", "", e, null, { width: "26px", height: "26px", borderRadius: "50%", background: col });
      mk("div", "kick", e, l, { color: "#fff" });
      return { e, o, n };
    });
    const rows = [
      [true, "DNA loops", `F1 ${HI.loops.hmec.chronocell.toFixed(2)} / ${HI.loops.hap1.chronocell.toFixed(2)} vs best tool ${Math.max(HI.loops.hmec.mustache, HI.loops.hmec.chromosight).toFixed(2)} / ${Math.max(HI.loops.hap1.mustache, HI.loops.hap1.chromosight).toFixed(2)}`],
      [true, "Molecules", `${HI.q6d.within} / ${HI.q6d.cases} within chemical accuracy`],
      [true, "Noisy chip", `${HI.q9.within} / ${HI.q9.cases} after mitigation, ${HI.q9.reduction}× smaller error`],
      [false, "Docking", `${HI.q7b.qaoa} % vs ${HI.q7b.random} % for random search`],
      [false, "Drug safety", `${HI.q8.passed} / ${HI.q8.of} endpoints, 17 needed`]];
    s.rows = rows.map(([ok, a, b], k) => {
      const e = mk("div", "abs", u, null, { left: "1000px", top: `${260 + k * 128}px`, width: "820px", height: "108px", borderRadius: "22px", background: "#fff", color: C.ink, padding: "18px 26px", boxShadow: "0 16px 40px rgba(0,0,0,.25)" });
      mk("div", "abs", e, ok ? "✓" : "✕", { right: "24px", top: "28px", width: "52px", height: "52px", borderRadius: "50%", background: ok ? C.pass : C.fail, color: "#fff", font: "900 30px/52px var(--sans)", textAlign: "center" });
      mk("div", "dc", e, a, { fontSize: "40px" });
      mk("div", "lead", e, b, { fontSize: "23px", marginTop: "6px", color: C.ink2 });
      cue(T.results + 1.6 + k * 0.3, "pop");
      return e;
    });
    cue(T.results + 0.8, "roll", { d: 1.2 });
  },
  update(s, t) {
    pop(s.h, t, T.results + 0.3, { s: 0.8 });
    s.big.forEach((b, k) => { pop(b.e, t, T.results + 0.6 + k * 0.15, { s: 0.7, x: -40, y: 0 }); b.o.set(b.n * E.o3(lin(t, T.results + 0.8, T.results + 2.0))); });
    s.rows.forEach((e, k) => pop(e, t, T.results + 1.6 + k * 0.3, { s: 0.6, x: 60, y: 0 }));
  } });

// ================================================================================================= outro
cap(T.outro + 0.6, T.end - 1.2, "Try it yourself: three commands, and the app opens in your browser.");
scene({ t0: T.outro, t1: T.end + 0.1, inT: { type: "fade", at: T.outro, d: 0.6 },
  build(s) {
    bgDark(s, C.cobalt, 50, 35);
    const u = s.ui;
    s.logo = mk("div", "abs logo", u, logoSvg(110, "#fff", C.ochre, 2.2), { left: "905px", top: "90px" });
    s.name = mk("div", "abs dn", u, "ChronoCell-5D", { left: "0", width: "1920px", top: "220px", textAlign: "center", fontSize: "130px", color: "#fff", textShadow: extrude(C.cobaltDk, 10) });
    s.term = mk("div", "term", u, null, { left: "360px", top: "430px", width: "1200px", height: "300px" });
    const bar = mk("div", "bar", s.term, "<i></i><i></i><i></i>");
    mk("div", "ttl", bar, "Terminal");
    mk("div", "abs cmd", s.term, '<span class="p">$</span> git clone https://github.com/Sh1voham/ChronoCell-5D.git<br><span class="p">$</span> pip install -r requirements.txt<br><span class="p">$</span> streamlit run app.py',
      { left: "36px", top: "74px", whiteSpace: "nowrap", fontSize: "29px" });
    s.line = mk("div", "abs mono", u, "EVERY CLAIM A TEST · RUN ONCE · FAILURES KEPT", { left: "0", width: "1920px", textAlign: "center", top: "780px", font: "600 26px/1 var(--mono)", letterSpacing: ".22em", color: "rgba(255,255,255,.85)" });
    cue(T.outro + 0.4, "boom");
  },
  update(s, t) {
    pop(s.logo, t, T.outro + 0.4, { s: 0.2, y: 0 });
    pop(s.name, t, T.outro + 0.6, { s: 0.8, y: 40 });
    pop(s.term, t, T.outro + 1.2, { s: 0.8, y: 60 });
    rise(s.line, t, T.outro + 2.0);
  } });

// ================================================================================================= HUD
const CHAP = [[T.title, null], [T.prob, "01", "What it is"], [T.tour, "02", "The app, live"], [T.test, "03", "How it is tested"],
  [T.results, "04", "Where it stands"], [T.outro, null], [T.end + 1, null]];
CHAP.forEach((c) => {
  if (!c[1]) return;
  c.el = mk("div", "chap", L.hud);
  mk("div", "n", c.el, c[1]);
  mk("div", "t", c.el, c[2]);
});
const capEl = mk("div", "cap", L.hud);
const prog = mk("div", "prog", L.hud);
const progBar = mk("i", "", prog);
[T.prob, T.tour, T.test, T.results, T.outro].forEach((x) => mk("s", "", prog, null, { left: `${(x / T.end) * 1920}px` }));
const black = mk("div", "black", L.ov);
let capKey = null;
function hud(t) {
  CHAP.forEach((c, i) => {
    if (!c.el) return;
    const next = CHAP[i + 1][0], on = t >= c[0] && t < next;
    c.el.style.display = on ? "flex" : "none";
    if (on) setT(c.el, `translate3d(${((1 - E.o3(lin(t, c[0] + 0.3, c[0] + 0.8))) * -40).toFixed(1)}px,0,0)`, lin(t, c[0] + 0.3, c[0] + 0.7) * (1 - lin(t, next - 0.3, next)));
  });
  const c = CAPS.find(([a, b]) => t >= a && t < b);
  if (!c) capEl.style.opacity = 0;
  else {
    if (capKey !== c) { capEl.innerHTML = c[2]; capKey = c; }
    setT(capEl, `translate3d(-50%,${((1 - E.o3(lin(t, c[0], c[0] + 0.35))) * 18).toFixed(1)}px,0)`, lin(t, c[0], c[0] + 0.25) * (1 - lin(t, c[1] - 0.2, c[1])));
  }
  progBar.style.transform = `scaleX(${(t / T.end).toFixed(5)})`;
  prog.style.opacity = (lin(t, 1.0, 2.0) * (1 - lin(t, T.end - 2, T.end - 1))).toFixed(3);
  black.style.opacity = Math.max(1 - lin(t, 0, 0.35), lin(t, T.end - 1.2, T.end - 0.05)).toFixed(3);
}

// ================================================================================================= seek + boot
async function seek(t) {
  const waits = [];
  for (const s of SCENES) {
    const on = t >= s.t0 && t < s.t1;
    show(s.root, on);
    if (on) { trans(s, t); const r = s.update(s, t); if (r && r.then) waits.push(r); }
  }
  hud(t);
  await Promise.all(waits);
  return t;
}
window.seek = seek;
window.DURATION = T.end;
window.CUES = CUES.sort((a, b) => a.t - b.t);
window.SECTIONS = [[0, 9, 0.55, false, true], [9, 41, 0.45, false, true], [41, 167, 0.6, true, true], [167, 187, 0.45, false, true],
  [187, 252, 0.6, true, true], [252, 274, 0.75, true, true], [274, 288, 0.4, false, true]];
window.READY = (async () => {
  await document.fonts.ready;
  await Promise.all(["900 100px Archivo", "600 20px 'IBM Plex Mono'", "500 20px 'Inter Tight'"].map((f) => document.fonts.load(f)));
  const vids = [...document.querySelectorAll("video")];
  await Promise.all(vids.map((v) => (v.readyState >= 2 ? null : new Promise((ok) => { v.addEventListener("loadeddata", ok, { once: true }); v.load(); setTimeout(ok, 15000); }))));
  await seek(0);
  return true;
})();
(async () => {
  const q = new URLSearchParams(location.search);
  await window.READY;
  if (q.has("t")) await seek(parseFloat(q.get("t")));
  if (q.has("play")) {
    const t0 = performance.now() - 1000 * parseFloat(q.get("play") || "0");
    const loop = async () => { await seek(((performance.now() - t0) / 1000) % T.end); requestAnimationFrame(loop); };
    loop();
  }
})();
