/* ChronoCell-5D film: the timeline. Every element's state is a pure function of the film time t, so
 * window.seek(t) can draw any frame in any order (motion/render.py calls it once per frame).
 * Facts on screen come from window.DATA (motion/build_data.py): the Scoreboard rows and the result files. */
"use strict";

// ------------------------------------------------------------------------------------------------- timing (s)
const TIMES = { open: 0, fold: 7.5, title: 16, tour: 21, ws: 3.4, mosaic: 48.2, claim: 52.2, count: 56.6, loops: 62.6,
  mol: 66.6, mit: 71.4, fail: 74.8, burst: 77.8, cut: 0.5, outro: 84.4, end: 92 };
const FPS = 60;

const D = window.DATA;
const clamp = (x, a = 0, b = 1) => Math.min(b, Math.max(a, x));
const lin = (t, a, b) => (b === a ? (t >= b ? 1 : 0) : clamp((t - a) / (b - a)));
const mix = (a, b, p) => a + (b - a) * p;
const E = {
  o2: (x) => 1 - (1 - x) * (1 - x),
  o3: (x) => 1 - Math.pow(1 - x, 3),
  o4: (x) => 1 - Math.pow(1 - x, 4),
  i3: (x) => x * x * x,
  io3: (x) => (x < 0.5 ? 4 * x * x * x : 1 - Math.pow(-2 * x + 2, 3) / 2),
  oExpo: (x) => (x >= 1 ? 1 : 1 - Math.pow(2, -10 * x)),
  iExpo: (x) => (x <= 0 ? 0 : Math.pow(2, 10 * x - 10)),
  ioExpo: (x) => (x <= 0 ? 0 : x >= 1 ? 1 : x < 0.5 ? Math.pow(2, 20 * x - 10) / 2 : (2 - Math.pow(2, -20 * x + 10)) / 2),
  oBack: (x) => { const c1 = 1.5, c3 = c1 + 1; return 1 + c3 * Math.pow(x - 1, 3) + c1 * Math.pow(x - 1, 2); },
};
const fmtInt = (n) => Math.round(n).toLocaleString("en-US");

// ------------------------------------------------------------------------------------------------- DOM helpers
function mk(tag, cls, parent, html, style) {
  const e = document.createElement(tag);
  if (cls) e.className = cls;
  if (html != null) e.innerHTML = html;
  if (style) Object.assign(e.style, style);
  if (parent) parent.appendChild(e);
  return e;
}
const IMG = {};
function img(name) { return name.includes("/") ? name : `assets/app/${name}.png`; }
function chars(el, text, cls = "ch") {     // one span per character (spaces kept)
  el.innerHTML = "";
  const out = [];
  for (const c of text) {
    const s = mk("span", cls, el, c === " " ? "&nbsp;" : c);
    out.push(s);
  }
  return out;
}
function words(el, text, gradWords = []) {  // one span per word; listed words get the gradient
  el.innerHTML = "";
  const out = [];
  text.split(" ").forEach((w, i, a) => {
    const s = mk("span", "w" + (gradWords.includes(w.replace(/[.,]/g, "")) ? " grad" : ""), el, w);
    out.push(s);
    if (i < a.length - 1) el.appendChild(document.createTextNode(" "));
  });
  return out;
}
function setT(el, tr, op, blur) {
  el.style.transform = tr;
  if (op != null) el.style.opacity = op.toFixed(4);
  el.style.filter = blur && blur > 0.05 ? `blur(${blur.toFixed(2)}px)` : "none";
}
/* in/out reveal of a list of spans: rise in from below, leave upwards */
function animSpans(spans, lt, o) {
  const n = spans.length;
  spans.forEach((s, i) => {
    const a = o.in + i * (o.st ?? 0.02);
    const pin = (o.ein ?? E.oExpo)(lin(lt, a, a + (o.din ?? 0.7)));
    const b = o.out != null ? o.out + (o.rev ? n - 1 - i : i) * (o.sto ?? 0.012) : 1e9;
    const pout = (o.eout ?? E.iExpo)(lin(lt, b, b + (o.dout ?? 0.45)));
    const y = (1 - pin) * (o.fy ?? 60) - pout * (o.ty ?? 50);
    const sc = mix(o.fs ?? 1, 1, pin) * mix(1, o.ts ?? 1, pout);
    const rot = (1 - pin) * (o.fr ?? 0);
    setT(s, `translate3d(0,${y.toFixed(2)}px,0) rotate(${rot.toFixed(2)}deg) scale(${sc.toFixed(4)})`,
      pin * (1 - pout), (1 - pin) * (o.fb ?? 8) + pout * (o.tb ?? 8));
  });
}
function show(el, on) { el.style.display = on ? "block" : "none"; }

// ------------------------------------------------------------------------------------------------- stage
const stage = document.getElementById("stage");
const ui = mk("div", "layer persp", stage);
const fx = mk("div", "layer", stage, null, { pointerEvents: "none" });
const SCENES = [];
const aurora = mk("div", "layer", ui, null, { mixBlendMode: "screen", pointerEvents: "none" });
const blobs = [["rgba(51,64,209,.55)", 1300], ["rgba(138,92,255,.42)", 1100], ["rgba(224,87,58,.30)", 1000]].map(([c, d]) =>
  mk("div", "", aurora, null, { position: "absolute", width: d + "px", height: d + "px", borderRadius: "50%",
    background: `radial-gradient(circle, ${c} 0%, rgba(0,0,0,0) 62%)` }));
function scene(t0, t1, build, update) {
  const root = mk("div", "scene", ui);
  const s = { t0, t1, root, update };
  build(root, s);
  SCENES.push(s);
  return s;
}

// ================================================================================================= 1. opening
scene(TIMES.open, TIMES.fold + 0.3, (r, s) => {
  const a = mk("div", "center", r, null, { top: "735px" });
  s.a1 = chars(mk("div", "h-s", a, null, { color: "#c9cdf0" }), "Every human cell packs");
  const big = mk("div", "h-xl", a, null, { marginTop: "14px" });
  s.a2 = words(big, "2 metres of DNA", ["2", "metres"]);
  const b = mk("div", "center", r, null, { top: "735px" });
  s.b1 = chars(mk("div", "h-s", b, null, { color: "#c9cdf0" }), "into a nucleus about");
  const big2 = mk("div", "h-xl", b, null, { marginTop: "14px" });
  s.b2 = words(big2, "10 micrometres wide.", ["10", "micrometres"]);
  s.kick = mk("div", "center kicker", r, "ChronoCell-5D · 2026", { top: "96px" });
}, (s, t) => {
  animSpans(s.a1, t, { in: 0.9, st: 0.025, out: 3.6, fy: 40 });
  animSpans(s.a2, t, { in: 1.25, st: 0.09, din: 0.9, fy: 90, fs: 1.25, fb: 14, out: 3.75 });
  animSpans(s.b1, t, { in: 4.25, st: 0.022, out: 7.0, fy: 40 });
  animSpans(s.b2, t, { in: 4.6, st: 0.09, din: 0.9, fy: 90, fs: 1.25, fb: 14, out: 7.1 });
  const k = lin(t, 0.3, 1.2) * (1 - lin(t, 6.8, 7.4));
  setT(s.kick, `translateY(${(1 - k) * -12}px)`, k * 0.8, 0);
});

// ================================================================================================= 2. the fold
scene(TIMES.fold, TIMES.title + 0.2, (r, s) => {
  const box = mk("div", "", r, null, { position: "absolute", left: "120px", top: "150px" });
  s.k = mk("div", "kicker", box, "How a chromosome folds");
  const h1 = mk("div", "h-l", box, null, { marginTop: "26px" });
  s.l1 = words(mk("span", "line", h1), "How it folds decides");
  s.l2 = words(mk("span", "line", h1), "which genes switch on.", ["genes", "switch", "on"]);
  const F = D.fold;
  const chipBox = mk("div", "", r, null, { position: "absolute", left: "120px", top: "880px", display: "flex", gap: "16px" });
  s.chips = [
    `<span class="dot"></span><b>${F.chrom}</b> · human, GRCh38`,
    `<b>${fmtInt(F.beads)}</b> beads × <b>${F.resolution_kb} kb</b>`,
    `<b>${F.size_mb}</b> Mb of DNA`,
    `<b>${fmtInt(F.contacts)}</b> contacts`,
  ].map((h) => mk("div", "chip", chipBox, h));
  s.note = mk("div", "mono", r, "The app's synthetic reference model of chr22, drawn from its own coordinates", {
    position: "absolute", right: "120px", top: "960px", font: "400 18px/1 var(--mono)", color: "#7d82a3", letterSpacing: ".04em" });
}, (s, t) => {
  const lt = t - TIMES.fold;
  const k = lin(lt, 2.4, 3.0) * (1 - lin(lt, 7.6, 8.1));
  setT(s.k, `translateX(${(1 - E.o3(k)) * -30}px)`, k, 0);
  animSpans(s.l1, lt, { in: 2.6, st: 0.07, din: 0.8, fy: 70, out: 7.7, sto: 0.03 });
  animSpans(s.l2, lt, { in: 3.0, st: 0.07, din: 0.8, fy: 70, out: 7.75, sto: 0.03 });
  s.chips.forEach((c, i) => {
    const p = E.oBack(lin(lt, 4.6 + i * 0.16, 5.3 + i * 0.16)), q = E.iExpo(lin(lt, 7.6 + i * 0.05, 8.1 + i * 0.05));
    setT(c, `translate3d(0,${(1 - p) * 40 + q * 30}px,0)`, clamp(p) * (1 - q), 0);
  });
  const n = lin(lt, 5.4, 6.0) * (1 - lin(lt, 7.6, 8.0));
  setT(s.note, "none", n * 0.9, 0);
});

// ================================================================================================= 3. title
const LOGO = `<svg width="150" height="150" viewBox="0 0 26 26"><defs><linearGradient id="lg" x1="0" x2="1">
  <stop offset="0" stop-color="#7f8cff"/><stop offset=".5" stop-color="#b38cff"/><stop offset="1" stop-color="#ff8a5c"/></linearGradient></defs>
  <path class="p1" d="M5 20c3-9 5-13 8-13s5 4 8 13" stroke="#f4f3ee" stroke-width="1.6"/>
  <path class="p2" d="M5 6c3 9 5 13 8 13s5-4 8-13" stroke="url(#lg)" stroke-width="1.6"/></svg>`;
scene(TIMES.title - 0.1, TIMES.tour + 0.1, (r, s) => {
  s.logo = mk("div", "center logo", r, LOGO, { top: "250px" });
  s.paths = [...s.logo.querySelectorAll("path")];
  s.paths.forEach((p) => { const L = p.getTotalLength(); p.style.strokeDasharray = L; p.dataset.len = L; });
  const name = mk("div", "center", r, null, { top: "440px", font: "800 168px/1 var(--sans)", letterSpacing: "-0.055em" });
  s.name = chars(name, "ChronoCell-5D");
  s.name.slice(-2).forEach((c) => c.classList.add("grad"));
  s.sub = chars(mk("div", "center h-s", r, null, { top: "650px", color: "#c9cdf0", fontWeight: 500 }),
    "The 3D / 4D genome workstation");
  s.list = mk("div", "center kicker", r, "Structure · Dynamics · Drugs · Genes · Quantum · Evidence", { top: "760px", color: "#8f96d8" });
}, (s, t) => {
  const lt = t - TIMES.title;
  s.paths.forEach((p, i) => {
    const q = E.io3(lin(lt, 0.05 + i * 0.15, 1.15 + i * 0.15));
    p.style.strokeDashoffset = ((1 - q) * p.dataset.len).toFixed(3);
  });
  const lo = E.iExpo(lin(lt, 4.4, 5.0));
  setT(s.logo, `scale(${mix(0.8, 1, E.oExpo(lin(lt, 0, 1.0))) * (1 + lo * 2.2)})`, lin(lt, 0, 0.2) * (1 - lo), lo * 20);
  s.name.forEach((c, i) => {
    const a = 0.35 + i * 0.045, p = E.oExpo(lin(lt, a, a + 0.9));
    const q = E.iExpo(lin(lt, 4.35 + Math.abs(i - 6) * 0.015, 4.95));
    setT(c, `translate3d(0,${(1 - p) * 120}px,${(1 - p) * -400 + q * 900}px) rotateX(${(1 - p) * -80}deg)`, p * (1 - q), (1 - p) * 10 + q * 24);
  });
  animSpans(s.sub, lt, { in: 1.3, st: 0.018, fy: 30, out: 4.3, sto: 0.006, ty: -30 });
  const l = lin(lt, 2.0, 2.6) * (1 - lin(lt, 4.2, 4.6));
  setT(s.list, `translateY(${(1 - l) * 16}px)`, l, 0);
});

// ================================================================================================= 4. the tour
const WS = [
  { n: "01", label: "3D structure", img: "ws01_wide", focus: [0.62, 0.42], cards: [["fold3d", "Reconstructed fold"], ["ws01_plot1", "Polymer scaling law"]] },
  { n: "02", label: "4D dynamics", img: "ws02_top", focus: [0.4, 0.5], gl: true, cards: [["ws02_plot1", "Relaxation"]], slots: [[1490, 250, 330, 200]] },
  { n: "03", label: "Compare", img: "ws03_top", focus: [0.5, 0.7], cards: [["ws03_fold", "Two states"], ["ws03_plot0", "Where they differ"]] },
  { n: "04", label: "Drug lab", img: "ws04_top", focus: [0.5, 0.55], cards: [["ws04_plot1", "Dose–response"], ["ws04_plot0", "Treated fold"]] },
  { n: "05", label: "Genes", img: "ws05_top", focus: [0.4, 0.5], cards: [["genes3d", "Open vs buried genes"], ["ws05_fold", "Gene table"]] },
  { n: "06", label: "Guide", img: "ws06_top", focus: [0.3, 0.6], scroll: "ws06_pages", cards: [["ws06_quantum", "Plain words"]] },
  { n: "07", label: "Quantum lab", img: "ws07_top", focus: [0.4, 0.6], cards: [["q_lattice_plot0", "QAOA"], ["q_noise_plot0", "Error mitigation"], ["q_walk_plot0", "Quantum walk"]] },
  { n: "08", label: "Scoreboard", img: "ws08_top", focus: [0.6, 0.55], cards: [["sb_0", "Before → after"], ["sb_2", "DNA loops"]] },
];
// each workspace's sentence is the app's own one-line purpose (app.py WORKSPACES)
const PURPOSE = {
  "01": "See how one chromosome is folded inside the nucleus, and measure it.",
  "02": "Watch the fold change over time, between conditions, or after a DNA rearrangement.",
  "03": "Put two biological states side by side; rotating one rotates the other.",
  "04": "Apply a virtual epigenetic drug and see how far it pushes the fold back toward healthy.",
  "05": "Find which genes sit in open, active chromatin and which are buried and likely silenced.",
  "06": "What everything means, in plain words, with a 2-minute tour.",
  "07": "Try ChronoCell's problems on a simulated quantum computer, next to the classical answer.",
  "08": "Every accuracy test on held-out real data: what passed, what failed, and by how much.",
};
const CARD_SLOTS = [   // where the floating close-ups sit: [left, top, width, z]
  [720, 600, 520, 160], [1400, 96, 430, 220], [1330, 700, 480, 260]];

const tour = scene(TIMES.tour - 0.45, TIMES.mosaic + 0.1, (r, s) => {
  s.items = WS.map((w, i) => {
    const g = mk("div", "", r, null, { position: "absolute", inset: "0", transformStyle: "preserve-3d" });
    const num = mk("div", "ws-num", g, w.n, { left: "92px", top: "170px" });
    const tag = mk("div", "ws-tag", g, `WORKSPACE ${w.n} / 08`, { left: "120px", top: "520px" });
    const lab = mk("div", "ws-label", g, null, { left: "116px", top: "560px" });
    const labC = chars(lab, w.label);
    const txt = mk("div", "ws-text", g, null, { left: "120px", top: "660px" });
    const txtW = words(txt, PURPOSE[w.n]);
    const big = !w.gl;
    const win = mk("div", "win", g, null, big ? { left: "760px", top: "150px", width: "1060px", height: "640px" }
      : { left: "1260px", top: "596px", width: "560px", height: "359px" });
    const bar = mk("div", "bar", win, "<i></i><i></i><i></i>");
    mk("div", "url", bar, `localhost:8501 · ChronoCell-5D · ${w.n} ${w.label}`);
    const vp = mk("div", "vp", win);
    const im = mk("img", "", vp);
    im.src = img(w.img);
    let im2 = null;
    if (w.scroll) { im2 = mk("img", "", vp); im2.src = img(w.scroll); im2.style.top = "100%"; }
    const sheen = mk("div", "sheen", win);
    const cards = (w.cards || []).map(([name, label], k) => {
      const c = mk("div", "card", g);
      const ci = mk("img", "", c); ci.src = img(name);
      mk("div", "tag", c, label);
      return { el: c, img: ci, slot: (w.slots || CARD_SLOTS)[k], name };
    });
    let cap = null;
    if (w.gl) {
      cap = mk("div", "chip", g, `<span class="dot" style="background:#ff7a4a;box-shadow:0 0 12px #ff7a4a"></span>` +
        `<b>${D.sv.title}</b> · ${D.sv.frames} frames · simulated in the app`, { position: "absolute", left: "860px", top: "150px" });
    }
    return { w, g, num, tag, lab, labC, txt, txtW, win, im, im2, vp, sheen, cards, cap, i };
  });
  s.streak = mk("div", "streak", r);
  const rl = mk("div", "rail-lab", r); s.railLab = WS.map((w) => mk("div", "", rl, `${w.n} ${w.label}`));
  const rail = mk("div", "rail", r); s.rail = WS.map(() => mk("span", "", mk("div", "", rail)));
}, (s, t) => {
  const lt = t - TIMES.tour, W = TIMES.ws;
  const cur = clamp(Math.floor(lt / W), 0, WS.length - 1);
  s.items.forEach((it) => {
    const a = lt - it.i * W;                       // local time of this workspace
    const on = a > -0.42 && a < W + 0.05;          // each entrance overlaps the previous exit
    show(it.g, on);
    if (!on) return;
    const style = it.i % 4;
    const pin = E.oExpo(lin(a, -0.38, 0.45)), pout = E.iExpo(lin(a, W - 0.55, W));
    // ---- text column
    const nIn = E.o3(lin(a, -0.25, 0.3));
    setT(it.num, `translate3d(${(1 - nIn) * -80 - pout * 60}px,0,0)`, nIn * (1 - pout), 0);
    setT(it.tag, `translateY(${(1 - E.o3(lin(a, -0.15, 0.3))) * 20}px)`, lin(a, -0.15, 0.2) * (1 - pout), 0);
    animSpans(it.labC, a, { in: -0.1, st: 0.022, din: 0.6, fy: 80, out: W - 0.5, sto: 0.008, dout: 0.35 });
    animSpans(it.txtW, a, { in: 0.15, st: 0.025, din: 0.6, fy: 30, fb: 6, out: W - 0.48, sto: 0.004, dout: 0.3, ty: 20 });
    // ---- window: entrance style varies, then a slow drift, then the exit
    let tr = "", clip = "none", op = 1;
    const drift = lin(a, 0, W);
    const baseRY = it.w.gl ? -8 : mix(-13, -7, drift), baseZ = mix(0, 60, drift);
    if (style === 0) { tr = `translate3d(0,0,${mix(-1100, baseZ, pin)}px) rotateY(${mix(-38, baseRY, pin)}deg)`; op = pin; }
    else if (style === 1) { tr = `translate3d(${(1 - pin) * 900}px,0,${baseZ}px) rotateY(${mix(-70, baseRY, pin)}deg) skewX(${(1 - pin) * -8}deg)`; op = clamp(pin * 1.5); }
    else if (style === 2) { tr = `translate3d(0,${(1 - pin) * 420}px,${baseZ}px) rotateX(${(1 - pin) * 70}deg) rotateY(${baseRY}deg)`; op = clamp(pin * 1.4); }
    else { tr = `translate3d(0,0,${baseZ}px) rotateY(${baseRY}deg)`; clip = `circle(${(E.io3(lin(a, -0.35, 0.5)) * 140).toFixed(2)}% at 70% 40%)`; }
    if (it.i % 2 === 0) { tr = `translate3d(0,0,${pout * 700}px) ` + tr; op *= 1 - pout; }
    else { tr = `translate3d(${pout * -1400}px,0,0) rotateY(${pout * 35}deg) ` + tr; op *= 1 - pout * 0.9; }
    it.win.style.transform = tr; it.win.style.opacity = op.toFixed(4); it.win.style.clipPath = clip;
    // Ken Burns inside the window; the Guide scrolls down its page instead
    const z = mix(1.0, 1.45, E.io3(lin(a, 0.3, W)));
    const [fxp, fyp] = it.w.focus;
    const vw = it.vp.clientWidth || 1060, vh = it.vp.clientHeight || 596;
    if (it.im2) {
      const sc = E.io3(lin(a, 0.9, W - 0.3));
      const ih = it.im.clientHeight || vh;
      it.im.style.transform = `translate3d(0,${-sc * ih}px,0)`;
      it.im2.style.transform = `translate3d(0,${-sc * ih}px,0)`;
    } else {
      it.im.style.transform = `translate3d(${(-(z - 1) * vw * fxp).toFixed(2)}px,${(-(z - 1) * vh * fyp).toFixed(2)}px,0) scale(${z.toFixed(4)})`;
    }
    it.sheen.style.transform = `translateX(${mix(-120, 120, E.io3(lin(a, 0.2, 1.4)))}%)`;
    // ---- floating close-ups
    it.cards.forEach((c, k) => {
      const [L, T, Wd, Z] = c.slot;
      const asp = (c.img.naturalWidth || 16) / (c.img.naturalHeight || 9);
      c.el.style.left = L + "px"; c.el.style.top = T + "px"; c.el.style.width = Wd + "px";
      c.el.style.height = Math.min(Wd / asp, 340).toFixed(1) + "px";
      const b = 0.3 + k * 0.28, p = E.oBack(lin(a, b, b + 0.75)), q = E.iExpo(lin(a, W - 0.6 + k * 0.03, W - 0.15));
      const fl = Math.sin((a + k) * 1.7) * 8;
      c.el.style.transform = `translate3d(${q * (k % 2 ? 300 : -300)}px,${(1 - clamp(p)) * 120 + fl}px,${mix(-500, Z, clamp(p)) + q * 400}px) rotateY(${(1 - clamp(p)) * 30 - 6}deg) rotateX(${(1 - clamp(p)) * -12}deg) scale(${mix(0.7, 1, clamp(p))})`;
      c.el.style.opacity = (clamp(p * 1.6) * (1 - q)).toFixed(4);
    });
    if (it.cap) { const c = E.o3(lin(a, 0.7, 1.2)) * (1 - pout); setT(it.cap, `translateY(${(1 - c) * 20}px)`, c, 0); }
  });
  // light streak across each change of workspace
  const ph = lt / W, frac = ph - Math.floor(ph), near = frac < 0.5 ? frac : frac - 1;   // time relative to a boundary
  const sx = mix(-700, 2300, clamp((near * W + 0.35) / 0.7));
  const sop = Math.exp(-Math.pow((near * W) / 0.25, 2)) * (lt > 0.1 && lt < W * 8 - 0.1 ? 1 : 0);
  setT(s.streak, `translateX(${sx}px) rotate(18deg)`, sop * 0.9, 0);
  // progress rail
  const rin = lin(lt, 0, 0.6) * (1 - lin(lt, W * 8 - 0.4, W * 8));
  s.rail.forEach((b, i) => { b.style.transform = `scaleX(${clamp((lt - i * W) / W).toFixed(4)})`; });
  s.rail[0].parentNode.parentNode.style.opacity = rin;
  s.railLab.forEach((l, i) => { l.style.color = i === cur ? "rgba(255,255,255,.9)" : "rgba(255,255,255,.3)"; });
  s.railLab[0].parentNode.style.opacity = rin;
});

// ================================================================================================= 5. mosaic
const TILES = ["ws01_wide", "q_lattice_view", "ws04_top", "ws06_quantum", "ws02_fold", "ws08_list", "ws03_fold", "ws05_top",
  "q_noise_view", "ws06_top", "ws01_fold", "q_walk_view", "ws07_top", "ws04_mid", "ws08_top", "ws05_fold", "ws03_top",
  "q_vqe_view", "ws06_pages", "ws02_top", "ws06_new", "ws01_top", "q_combo_view", "q_tad_view", "ws04_mid", "fold3d",
  "ws08_list", "ws02_top", "genes3d", "ws01_wide"];
scene(TIMES.mosaic - 0.1, TIMES.claim + 0.1, (r, s) => {
  const wrap = mk("div", "layer persp", r, null, { perspective: "1800px" });
  s.m = mk("div", "mosaic", wrap);
  s.cols = [];
  for (let c = 0; c < 6; c++) {
    const col = mk("div", "mcol", s.m, null, { left: `${(c - 3) * 590}px`, top: "-1200px" });
    for (let k = 0; k < 5; k++) { const tl = mk("div", "tile", col); const im = mk("img", "", tl); im.src = img(TILES[(c * 5 + k) % TILES.length]); }
    s.cols.push(col);
  }
  const box = mk("div", "center", r, null, { top: "400px" });
  s.l1 = words(mk("div", "h-xl", box), "8 workspaces.", ["8"]);
  s.l2 = words(mk("div", "h-m", box, null, { marginTop: "22px", color: "#dfe2ff" }), "One genome, in 3D and 4D.", ["3D", "4D"]);
  s.dim = mk("div", "layer", r, null, { background: "radial-gradient(ellipse at 50% 50%, rgba(4,5,11,.78), rgba(4,5,11,.35) 70%)" });
  r.appendChild(box);
}, (s, t) => {
  const lt = t - TIMES.mosaic, L = TIMES.claim - TIMES.mosaic;
  const pin = E.oExpo(lin(lt, 0, 1.0)), pout = E.iExpo(lin(lt, L - 0.6, L));
  const sc = mix(2.3, 1.0, E.io3(lin(lt, 0, L))) * (1 + pout * 2.5);
  s.m.style.transform = `translateZ(${mix(-600, 0, pin)}px) rotateX(${mix(58, 48, lin(lt, 0, L))}deg) rotateZ(-26deg) scale(${sc.toFixed(4)})`;
  s.m.style.opacity = (pin * (1 - pout)).toFixed(4);
  s.cols.forEach((c, i) => { c.style.transform = `translate3d(0,${(i % 2 ? 1 : -1) * 240 * lt + (i % 2 ? -300 : 0)}px,0)`; });
  s.dim.style.opacity = (lin(lt, 0.6, 1.3) * (1 - pout)).toFixed(4);
  animSpans(s.l1, lt, { in: 0.8, st: 0.1, din: 0.8, fy: 80, fs: 1.3, fb: 14, out: L - 0.7, ts: 1.6, ty: 0, tb: 20 });
  animSpans(s.l2, lt, { in: 1.3, st: 0.06, din: 0.7, fy: 40, out: L - 0.65, ts: 1.4, ty: 0, tb: 20 });
});

// ================================================================================================= 6. the claim
scene(TIMES.claim - 0.05, TIMES.count + 0.1, (r, s) => {
  const box = mk("div", "center", r, null, { top: "330px" });
  s.l1 = words(mk("div", "h-xl", box), "Every claim is a test.", ["test"]);
  s.l2 = words(mk("div", "h-s", box, null, { marginTop: "40px", color: "#dfe2ff", fontWeight: 500 }), "Written down before its data were read.");
  s.l3 = words(mk("div", "h-m", box, null, { marginTop: "26px" }), "Run once. Failures kept.", ["once", "Failures", "kept"]);
}, (s, t) => {
  const lt = t - TIMES.claim, L = TIMES.count - TIMES.claim;
  const o = { st: 0.11, din: 0.55, fy: 0, fs: 1.7, fb: 16, eout: E.i3, dout: 0.4, ty: -40, tb: 10 };
  animSpans(s.l1, lt, { ...o, in: 0.15, out: L - 0.55 });
  animSpans(s.l2, lt, { ...o, in: 1.35, st: 0.06, fs: 1.3, out: L - 0.5 });
  animSpans(s.l3, lt, { ...o, in: 2.4, st: 0.16, out: L - 0.45 });
});

// ================================================================================================= 7. the count
const STATUS_COL = { pass: "#3ec27a", fail: "#e0573a", other: "#8b8fa8" };
scene(TIMES.count - 0.05, TIMES.loops + 0.1, (r, s) => {
  const T = D.tests;
  const left = mk("div", "", r, null, { position: "absolute", left: "140px", top: "250px" });
  s.k = mk("div", "kicker", left, "The Scoreboard");
  s.n = mk("div", "big-n grad", left, "0", { marginTop: "30px" });
  s.lab = mk("div", "count-lab", left, "tests on held-out real data,<br>each written down first and run once", { marginTop: "26px", width: "560px" });
  s.dots = T.rows.map((row) => mk("div", "tdot", r, null, { background: "#fff", left: "0", top: "0" }));
  s.rows = T.rows;
  const order = { pass: 0, fail: 1, other: 2 };
  const seen = { pass: 0, fail: 0, other: 0 };
  s.target = T.rows.map((row) => { const g = order[row.status], k = seen[row.status]++; return [g, k]; });
  s.groups = [["pass", T.pass, "passed"], ["fail", T.fail, "failed · kept on the record"], ["other", T.other, "mixed, blocked or baseline"]]
    .map(([k, n, l], g) => {
      const e = mk("div", "grp", r, null, { left: "1490px", top: `${238 + g * 210}px` });
      mk("div", "n", e, String(n), { color: STATUS_COL[k] });
      mk("div", "l", e, l, { color: STATUS_COL[k] });
      return e;
    });
}, (s, t) => {
  const lt = t - TIMES.count, L = TIMES.loops - TIMES.count;
  const out = E.iExpo(lin(lt, L - 0.5, L));
  const kk = lin(lt, 0.1, 0.5) * (1 - out);
  setT(s.k, "none", kk, 0);
  const c = E.o4(lin(lt, 0.3, 1.9));
  s.n.textContent = String(Math.round(c * D.tests.n));
  setT(s.n, `scale(${mix(0.8, 1, E.oBack(lin(lt, 0.2, 0.9)))})`, lin(lt, 0.2, 0.5) * (1 - out), (1 - lin(lt, 0.2, 0.6)) * 10 + out * 12);
  setT(s.lab, `translateY(${(1 - E.o3(lin(lt, 0.9, 1.5))) * 20}px)`, lin(lt, 0.9, 1.4) * (1 - out), 0);
  // dots: a grid in test order, then colour by verdict, then sorted into three groups
  const N = s.dots.length;
  s.dots.forEach((d, i) => {
    const gx = 820 + (i % 8) * 70, gy = 260 + Math.floor(i / 8) * 70;
    const [g, k] = s.target[i];
    const tx = 820 + (k % 8) * 70, ty = 250 + g * 210 + Math.floor(k / 8) * 70;
    const app = E.oBack(lin(lt, 0.4 + i * 0.035, 0.85 + i * 0.035));
    const col = lin(lt, 2.1 + i * 0.02, 2.4 + i * 0.02);
    const mv = E.io3(lin(lt, 2.9 + (i % 9) * 0.03, 3.8 + (i % 9) * 0.03));
    const x = mix(gx, tx, mv), y = mix(gy, ty, mv) + Math.sin(mv * Math.PI) * -40;
    d.style.background = col > 0.5 ? STATUS_COL[s.rows[i].status] : "#e9eaf5";
    d.style.boxShadow = col > 0.5 ? `0 0 18px ${STATUS_COL[s.rows[i].status]}` : "0 0 12px rgba(255,255,255,.5)";
    setT(d, `translate3d(${x}px,${y}px,0) scale(${clamp(app, 0, 1.3) * (1 - out)})`, clamp(app) * (1 - out), 0);
  });
  s.groups.forEach((g, i) => {
    const p = E.o3(lin(lt, 3.7 + i * 0.25, 4.3 + i * 0.25));
    setT(g, `translateX(${(1 - p) * 40}px)`, p * (1 - out), 0);
  });
});

// ================================================================================================= 8-11. highlights
function header(r, gate, ok, title, sub) {
  const box = mk("div", "", r, null, { position: "absolute", left: "140px", top: "110px", width: "1640px" });
  const g = mk("div", "gate", box, `${gate} · <span style="color:${ok ? "#3ec27a" : "#e0573a"}">${ok ? "PASS" : "FAIL"}</span>`);
  const h = mk("div", "h-m", box, null, { marginTop: "22px" });
  const tw = words(h, title, ok ? ["better", "accurate", "repaired"] : ["lost", "says", "so"]);
  const b = mk("div", "body", box, sub, { marginTop: "18px", width: "1300px" });
  return { g, tw, b };
}
function animHeader(hd, lt, L) {
  const out = E.iExpo(lin(lt, L - 0.45, L));
  setT(hd.g, `translateY(${(1 - E.o3(lin(lt, 0.05, 0.5))) * 16}px)`, lin(lt, 0.05, 0.35) * (1 - out), 0);
  animSpans(hd.tw, lt, { in: 0.12, st: 0.045, din: 0.6, fy: 60, out: L - 0.45, sto: 0.005, dout: 0.35 });
  setT(hd.b, `translateY(${(1 - E.o3(lin(lt, 0.55, 1.1))) * 16}px)`, lin(lt, 0.55, 0.95) * (1 - out), 0);
  return out;
}

// ---- loops (Gate 6b)
scene(TIMES.loops - 0.05, TIMES.mol + 0.1, (r, s) => {
  s.hd = header(r, "Gate 6b · DNA loop calls", true, "Finds DNA loops better than two standard tools",
    "On two cell types it had never seen (HMEC, HAP-1). Score: F1 against ENCODE's reference loops; higher is better.");
  const ch = mk("div", "chart", r, null, { left: "260px", top: "460px", width: "1400px", height: "470px" });
  mk("div", "axis", ch, null, { left: "0", right: "0", bottom: "0", height: "2px" });
  s.bars = [];
  const cells = Object.entries(D.hi.loops);
  const names = [["chronocell", "ChronoCell"], ["chromosight", "chromosight"], ["mustache", "Mustache"]];
  cells.forEach(([cell, v], gi) => {
    names.forEach(([k, lab], bi) => {
      const x = gi * 720 + bi * 200 + 40, h = v[k] * 430;
      const b = mk("div", "bar", ch, null, { left: `${x}px`, width: "160px", height: `${h}px`,
        background: k === "chronocell" ? "linear-gradient(180deg,#8f9bff,#3340d1)" : (bi === 1 ? "#4a4e66" : "#62667f"),
        boxShadow: k === "chronocell" ? "0 0 40px rgba(111,124,255,.55)" : "none" });
      const val = mk("div", "bar-v", b, "0.00", { top: "-46px", color: k === "chronocell" ? "#fff" : "#b9bcd3" });
      const l = mk("div", "bar-l", ch, lab, { left: `${x - 20}px`, width: "200px", bottom: "-40px" });
      s.bars.push({ b, val, v: v[k], gi, bi, l });
    });
    mk("div", "bar-l", ch, cell === "hmec" ? "HMEC (breast epithelium)" : "HAP-1 (near-haploid line)",
      { left: `${gi * 720 + 40}px`, width: "560px", bottom: "-84px", color: "#fff", fontSize: "21px" });
  });
}, (s, t) => {
  const lt = t - TIMES.loops, L = TIMES.mol - TIMES.loops;
  const out = animHeader(s.hd, lt, L);
  s.bars.forEach((x) => {
    const a = 0.9 + x.gi * 0.35 + x.bi * 0.12, p = E.o4(lin(lt, a, a + 1.1));
    x.b.style.transform = `scaleY(${p.toFixed(4)})`;
    x.b.style.opacity = (1 - out).toFixed(3);
    x.val.textContent = (x.v * p).toFixed(2);
    x.val.style.transform = `scaleY(${p > 0.02 ? (1 / p).toFixed(4) : 0})`;
    x.val.style.transformOrigin = "50% 100%";
    x.l.style.opacity = (lin(lt, a, a + 0.4) * (1 - out)).toFixed(3);
  });
  s.bars[0].l.parentNode.querySelectorAll(".bar-l").forEach((e, i) => { if (!s.bars.some((b) => b.l === e)) e.style.opacity = (lin(lt, 1.0, 1.5) * (1 - out)).toFixed(3); });
});

// ---- molecules (Gate Q6d)
scene(TIMES.mol - 0.05, TIMES.mit + 0.1, (r, s) => {
  const q = D.hi.q6d;
  s.hd = header(r, "Gate Q6d · quantum chemistry", true, `Chemically accurate, ${q.within} of ${q.cases} new molecules`,
    `Bonds stretched toward breaking, on a simulated quantum computer. Each bar is the energy error (log scale); all sit below chemical accuracy, 1.6 mHa. Worst: ${q.worst} mHa.`);
  const ch = mk("div", "chart", r, null, { left: "200px", top: "470px", width: "1520px", height: "440px" });
  mk("div", "axis", ch, null, { left: "0", right: "0", bottom: "0", height: "2px" });
  const y = (e) => ((Math.log10(Math.max(e, 1e-3)) + 3) / 4) * 440;           // 0.001 ... 10 mHa
  s.limit = mk("div", "limit", ch, `<span>chemical accuracy · ${D.hi.chem_accuracy} mHa</span>`, { left: "0", right: "0", bottom: `${y(D.hi.chem_accuracy)}px` });
  [0.001, 0.01, 0.1, 1, 10].forEach((v) => mk("div", "bar-l", ch, String(v), { left: "-90px", width: "70px", bottom: `${y(v) - 9}px`, textAlign: "right", fontSize: "16px" }));
  const sub = (m) => m.replace(/(\d)/g, "<sub>$1</sub>");
  s.bars = q.rows.map((row, i) => {
    const x = 30 + i * 106;
    const b = mk("div", "bar", ch, null, { left: `${x}px`, width: "62px", height: `${Math.max(y(row.err), 6)}px`,
      background: "linear-gradient(180deg,#7ff0b0,#1f8f55)", boxShadow: "0 0 26px rgba(62,194,122,.45)" });
    const l = mk("div", "bar-l", ch, `${sub(row.m)}<br><span style="color:#7d82a3">${row.s}×</span>`, { left: `${x - 24}px`, width: "110px", bottom: "-58px", fontSize: "17px" });
    return { b, l, i };
  });
  s.badge = mk("div", "", r, `<span class="check">✓</span>`, { position: "absolute", right: "150px", top: "120px" });
}, (s, t) => {
  const lt = t - TIMES.mol, L = TIMES.mit - TIMES.mol;
  const out = animHeader(s.hd, lt, L);
  const lp = E.io3(lin(lt, 0.8, 1.6));
  s.limit.style.clipPath = `inset(0 ${((1 - lp) * 100).toFixed(2)}% 0 0)`;
  s.limit.style.opacity = (1 - out).toFixed(3);
  s.bars.forEach((x) => {
    const a = 1.2 + x.i * 0.08, p = E.oBack(lin(lt, a, a + 0.8));
    x.b.style.transform = `scaleY(${clamp(p, 0, 1.15).toFixed(4)})`;
    x.b.style.opacity = (1 - out).toFixed(3);
    x.l.style.opacity = (lin(lt, a, a + 0.3) * (1 - out)).toFixed(3);
  });
  const bp = E.oBack(lin(lt, 2.9, 3.5));
  setT(s.badge, `scale(${clamp(bp, 0, 1.2)}) rotate(${(1 - clamp(bp)) * -90}deg)`, clamp(bp) * (1 - out), 0);
});

// ---- mitigation (Gate Q9)
scene(TIMES.mit - 0.05, TIMES.fail + 0.1, (r, s) => {
  const q = D.hi.q9;
  s.hd = header(r, "Gate Q9 · error mitigation", true, "A noisy quantum chip, repaired",
    `Simulated at today's best noise level, ${q.within} of ${q.cases} new molecules came back within chemical accuracy after symmetry checks and zero-noise extrapolation.`);
  s.svg = mk("div", "", r, `<svg width="1640" height="300" viewBox="0 0 1640 300"><path id="wv" fill="none" stroke-width="5" stroke-linecap="round"/></svg>`,
    { position: "absolute", left: "140px", top: "430px" });
  s.path = s.svg.querySelector("path");
  s.a = mk("div", "", r, null, { position: "absolute", left: "140px", top: "760px" });
  s.aN = mk("div", "mono", s.a, "", { font: "700 96px/1 var(--mono)", color: "#e0573a" });
  mk("div", "l mono", s.a, "median error, noisy chip (symmetry-checked)", { font: "500 22px/1.3 var(--mono)", color: "#e0573a", marginTop: "10px", letterSpacing: ".08em" });
  s.arrow = mk("div", "", r, "→", { position: "absolute", left: "740px", top: "760px", font: "300 110px/1 var(--sans)", color: "#7d82a3" });
  s.b = mk("div", "", r, null, { position: "absolute", left: "940px", top: "760px" });
  s.bN = mk("div", "mono", s.b, "", { font: "700 96px/1 var(--mono)", color: "#3ec27a" });
  mk("div", "l mono", s.b, "after zero-noise extrapolation", { font: "500 22px/1.3 var(--mono)", color: "#3ec27a", marginTop: "10px", letterSpacing: ".08em" });
  s.x = mk("div", "chip", r, `<b>${q.reduction}×</b> smaller`, { position: "absolute", left: "1500px", top: "790px", fontSize: "30px" });
}, (s, t) => {
  const lt = t - TIMES.mit, L = TIMES.fail - TIMES.mit, q = D.hi.q9;
  const out = animHeader(s.hd, lt, L);
  const fix = E.io3(lin(lt, 1.5, 2.6));
  let d = "";
  for (let i = 0; i <= 240; i++) {
    const x = (i / 240) * 1640, ph = i / 240 * Math.PI * 6;
    const n = Math.sin(i * 12.9898 + Math.floor(lt * 30) * 78.233) * 43758.5453;
    const noise = (n - Math.floor(n) - 0.5) * 2;
    const yv = 150 + Math.sin(ph + lt * 2) * 70 + noise * 95 * (1 - fix);
    d += (i ? "L" : "M") + x.toFixed(1) + " " + yv.toFixed(1);
  }
  s.path.setAttribute("d", d);
  s.path.setAttribute("stroke", fix > 0.5 ? "#3ec27a" : "#e0573a");
  s.path.style.filter = `drop-shadow(0 0 14px ${fix > 0.5 ? "rgba(62,194,122,.7)" : "rgba(224,87,58,.7)"})`;
  const wp = E.io3(lin(lt, 0.4, 1.2));
  s.svg.style.clipPath = `inset(0 ${((1 - wp) * 100).toFixed(2)}% 0 0)`;
  s.svg.style.opacity = (1 - out).toFixed(3);
  const an = E.o3(lin(lt, 0.6, 1.2));
  s.aN.textContent = (q.noisy * an).toFixed(1) + " mHa";
  setT(s.a, `translateY(${(1 - an) * 20}px)`, an * (1 - out), 0);
  const ar = lin(lt, 1.6, 2.0);
  setT(s.arrow, `translateX(${(1 - ar) * -30}px)`, ar * (1 - out), 0);
  const bn = E.o3(lin(lt, 1.9, 2.7));
  s.bN.textContent = mix(q.noisy, q.median, bn).toFixed(2) + " mHa";
  setT(s.b, `translateY(${(1 - E.o3(lin(lt, 1.9, 2.3))) * 20}px)`, lin(lt, 1.9, 2.2) * (1 - out), 0);
  const xp = E.oBack(lin(lt, 2.6, 3.0));
  setT(s.x, `scale(${clamp(xp, 0, 1.2)})`, clamp(xp) * (1 - out), 0);
});

// ---- the honest fail (Gate Q7b)
scene(TIMES.fail - 0.05, TIMES.burst + 0.1, (r, s) => {
  const q = D.hi.q7b;
  s.hd = header(r, "Gate Q7b · quantum docking", false, "And where it lost, it says so.",
    `On ${q.complexes} protein–drug complexes it had never seen, the quantum route placed ${q.qaoa} % of the drugs correctly; plain random search placed ${q.random} %. The result stays on the record.`);
  const ch = mk("div", "chart", r, null, { left: "140px", top: "520px", width: "1640px", height: "300px" });
  s.rows = [["Quantum route (QAOA)", q.qaoa, "linear-gradient(90deg,#ff8a5c,#e0573a)"], ["Random search, same score", q.random, "linear-gradient(90deg,#9ea2bd,#d8dae8)"]]
    .map(([lab, v, bg], i) => {
      mk("div", "bar-l", ch, lab, { left: "0", top: `${i * 150}px`, width: "560px", textAlign: "left", fontSize: "24px", color: "#dfe2ff" });
      const b = mk("div", "", ch, null, { position: "absolute", left: "0", top: `${i * 150 + 40}px`, height: "64px", width: `${v * 16}px`,
        background: bg, borderRadius: "8px", transformOrigin: "0 50%" });
      const val = mk("div", "mono", ch, "", { position: "absolute", top: `${i * 150 + 44}px`, font: "700 52px/1 var(--mono)", color: "#fff" });
      return { b, val, v, i };
    });
  s.badge = mk("div", "", r, `<span class="xmark">✕</span>`, { position: "absolute", right: "150px", top: "120px" });
}, (s, t) => {
  const lt = t - TIMES.fail, L = TIMES.burst - TIMES.fail;
  const out = animHeader(s.hd, lt, L);
  s.rows.forEach((x) => {
    const p = E.o4(lin(lt, 0.9 + x.i * 0.25, 2.0 + x.i * 0.25));
    x.b.style.transform = `scaleX(${p.toFixed(4)})`;
    x.b.style.opacity = (1 - out).toFixed(3);
    x.val.textContent = Math.round(x.v * p) + " %";
    x.val.style.left = `${x.v * 16 * p + 24}px`;
    x.val.style.opacity = (lin(lt, 0.9, 1.2) * (1 - out)).toFixed(3);
  });
  s.rows[0].b.parentNode.querySelectorAll(".bar-l").forEach((e) => (e.style.opacity = (lin(lt, 0.7, 1.1) * (1 - out)).toFixed(3)));
  const bp = E.oBack(lin(lt, 2.1, 2.6));
  setT(s.badge, `scale(${clamp(bp, 0, 1.2)})`, clamp(bp) * (1 - out), 0);
});

// ================================================================================================= 12. burst
const BURST = [["3D folds", "fold3d", 0.5, 0.5], ["4D dynamics", "assets/film/sv.jpg", 0.6, 0.5], ["Side by side", "ws03_fold", 0.4, 0.4],
  ["Virtual drugs", "ws04_plot0", 0.5, 0.5], ["Dose–response", "ws04_plot1", 0.5, 0.5], ["Gene activity", "genes3d", 0.5, 0.5],
  ["Loop calling", "sb_2", 0.5, 0.5], ["QAOA", "q_lattice_plot0", 0.5, 0.5], ["VQE", "q_vqe_view", 0.5, 0.35], ["Noisy chips, repaired", "q_noise_plot0", 0.5, 0.5],
  ["Quantum walks", "q_walk_view", 0.4, 0.4], ["Plain words", "ws06_top", 0.3, 0.4], [`${D.tests.n} tests`, "ws08_list", 0.5, 0.3]];
scene(TIMES.burst - 0.05, TIMES.outro + 0.15, (r, s) => {
  s.cuts = BURST.map(([word, name, fxp, fyp]) => {
    const g = mk("div", "layer", r);
    const box = mk("div", "bimg", g);
    const im = mk("img", "", box); im.src = img(name);
    mk("div", "bshade", g);
    const w = mk("div", "bword", g); const ws = chars(w, word);
    return { g, im, ws, fxp, fyp };
  });
  s.flash = mk("div", "flash", r);
}, (s, t) => {
  const lt = t - TIMES.burst, C = TIMES.cut;
  const k = Math.floor(lt / C);
  s.flash.style.opacity = 0;
  s.cuts.forEach((c, i) => {
    const last = i === s.cuts.length - 1;
    const a = lt - i * C, len = last ? TIMES.outro - TIMES.burst - i * C : C;
    const on = a >= 0 && a < len;
    show(c.g, on);
    if (!on) return;
    const iw = c.im.naturalWidth || 1600, ih = c.im.naturalHeight || 900;
    const cover = Math.max(1920 / iw, 1080 / ih);
    const z = cover * mix(1.12, 1.3, a / len) * (i % 2 ? 1 : 1.04);
    const x = clamp(960 - iw * z * c.fxp + mix(-30, 30, a / len) * (i % 2 ? 1 : -1), 1920 - iw * z, 0);
    const y = clamp(540 - ih * z * c.fyp, 1080 - ih * z, 0);                   // always covers the frame
    c.im.style.width = iw + "px"; c.im.style.height = ih + "px"; c.im.style.transformOrigin = "0 0";
    c.im.style.transform = `translate3d(${x.toFixed(1)}px,${y.toFixed(1)}px,0) scale(${z.toFixed(4)})`;
    animSpans(c.ws, a, { in: 0, st: 0.012, din: 0.22, fy: 0, fs: 1.6, fb: 10, out: last ? len - 0.35 : 1e9 });
    const g = (1 - lin(a, 0, 0.12)) * (last ? 1 : 0.85);
    c.g.style.filter = a < 0.08 ? `saturate(1.6) contrast(1.15) hue-rotate(${(1 - a / 0.08) * 40}deg)` : "none";
    if (a < 0.2) s.flash.style.opacity = g.toFixed(3);
    if (last) c.g.style.opacity = (1 - lin(a, len - 0.4, len)).toFixed(3);
    else c.g.style.opacity = 1;
  });
  if (k < 0) s.flash.style.opacity = 0;
});

// ================================================================================================= 13. outro
scene(TIMES.outro - 0.05, TIMES.end + 0.1, (r, s) => {
  s.logo = mk("div", "logo", r, LOGO.replace('width="150" height="150"', 'width="84" height="84"').replace('id="lg"', 'id="lg2"').replace("url(#lg)", "url(#lg2)"), { position: "absolute", left: "0", right: "0", top: "668px", textAlign: "center" });
  s.paths = [...s.logo.querySelectorAll("path")];
  s.paths.forEach((p) => { const L = p.getTotalLength(); p.style.strokeDasharray = L; p.dataset.len = L; });
  s.name = chars(mk("div", "center", r, null, { top: "770px", font: "800 104px/1 var(--sans)", letterSpacing: "-0.05em" }), "ChronoCell-5D");
  s.name.slice(-2).forEach((c) => c.classList.add("grad"));
  s.tag = words(mk("div", "center h-s", r, null, { top: "900px", color: "#dfe2ff", fontWeight: 500 }), "See the genome fold. Test every claim.", ["fold", "claim"]);
  s.fine = mk("div", "center mono", r, `${D.tests.n} tests on held-out real data · each written down before its data were read · run once · failures kept`,
    { top: "980px", font: "400 18px/1 var(--mono)", color: "#7d82a3", letterSpacing: ".04em" });
}, (s, t) => {
  const lt = t - TIMES.outro;
  s.paths.forEach((p, i) => { const q = E.io3(lin(lt, 0.6 + i * 0.12, 1.5 + i * 0.12)); p.style.strokeDashoffset = ((1 - q) * p.dataset.len).toFixed(3); });
  animSpans(s.name, lt, { in: 1.0, st: 0.04, din: 0.9, fy: 60, fb: 10 });
  animSpans(s.tag, lt, { in: 2.2, st: 0.08, din: 0.7, fy: 30 });
  const f = lin(lt, 3.2, 3.8);
  setT(s.fine, "none", f * 0.9, 0);
});

// ================================================================================================= overlays
const streakTop = mk("div", "flash", fx);          // title flash
const vignette = mk("div", "vignette", fx);
const grain = mk("div", "grain", fx);
const lbTop = mk("div", "letterbox", fx, null, { top: "0" });
const lbBot = mk("div", "letterbox", fx, null, { bottom: "0" });
const black = mk("div", "black", fx);
const WIPES = [TIMES.loops, TIMES.mol, TIMES.mit, TIMES.fail, TIMES.burst];
const wipeBars = ["#3340d1", "#8a5cff", "#e0573a"].map((c) => mk("div", "", fx, null, {
  position: "absolute", top: "-10%", height: "120%", width: "2300px", left: "-190px", background: c, opacity: "0" }));
{
  const c = document.createElement("canvas"); c.width = c.height = 256;
  const g = c.getContext("2d"), id = g.createImageData(256, 256);
  let s = 7;
  for (let i = 0; i < id.data.length; i += 4) { s = (s * 16807) % 2147483647; const v = s & 255; id.data[i] = id.data[i + 1] = id.data[i + 2] = v; id.data[i + 3] = 255; }
  g.putImageData(id, 0, 0);
  grain.style.backgroundImage = `url(${c.toDataURL()})`;
}
function overlays(t) {
  // diagonal colour wipes between the proof scenes (the scene changes behind the middle bar)
  const b = WIPES.find((w) => Math.abs(t - w) < 0.42);
  wipeBars.forEach((bar, k) => {
    if (b === undefined) { bar.style.opacity = "0"; return; }
    const p = clamp((t - (b - 0.42)) / 0.84);
    const q = E.io3(clamp(p * 1.3 - (k - 1) * 0.15 - 0.15));
    bar.style.opacity = "1";
    bar.style.transform = `translateX(${mix(-2600, 2600, q).toFixed(1)}px) skewX(-18deg)`;
  });
  // punch on hits and booms: the whole picture kicks and settles
  let pk = 0;
  for (const c of window.CUES || []) {
    const d = t - c.t;
    if (d < 0 || d > 0.7 || (c.k !== "hit" && c.k !== "boom")) continue;
    pk += (c.k === "boom" ? 0.03 : 0.014) * Math.exp(-d * 9);
  }
  stage.style.transform = pk > 1e-4 ? `scale(${(1 + pk).toFixed(5)})` : "none";
  // slow colour glow behind the title and the proof scenes
  const au = 0.9 * Math.max(lin(t, 16.0, 17.0) * (1 - lin(t, 20.4, 21.2)), lin(t, TIMES.claim - 0.4, TIMES.claim + 0.6) * (1 - lin(t, TIMES.burst - 0.3, TIMES.burst)));
  aurora.style.opacity = au.toFixed(3);
  if (au > 0.001) blobs.forEach((bl, i) => {
    const x = 960 + 620 * Math.sin(t * 0.23 + i * 2.1) - bl.offsetWidth / 2, y = 540 + 300 * Math.cos(t * 0.19 + i * 1.7) - bl.offsetHeight / 2;
    bl.style.transform = `translate3d(${x.toFixed(1)}px,${y.toFixed(1)}px,0)`;
  });
  const fl = Math.exp(-Math.pow((t - TIMES.title - 0.1) / 0.18, 2)) * 0.9
    + Math.exp(-Math.pow((t - TIMES.claim) / 0.15, 2)) * 0.35;
  streakTop.style.opacity = fl.toFixed(3);
  const f = Math.floor(t * FPS);
  grain.style.transform = `translate(${(f * 37) % 200 - 100}px, ${(f * 53) % 200 - 100}px)`;
  const lb = 70 * (1 - E.io3(lin(t, 0.2, 1.6))) + 70 * E.io3(lin(t, TIMES.end - 1.1, TIMES.end - 0.2));
  lbTop.style.height = lbBot.style.height = lb.toFixed(1) + "px";
  black.style.opacity = Math.max(1 - lin(t, 0, 0.6), lin(t, TIMES.end - 0.7, TIMES.end - 0.05)).toFixed(3);
}

// ================================================================================================= seek + boot
function seek(t) {
  for (const s of SCENES) {
    const on = t >= s.t0 && t <= s.t1;
    show(s.root, on);
    if (on) s.update(s, t);
  }
  overlays(t);
  GL.render(t);
  return t;
}
window.seek = seek;
window.TIMES = TIMES;
window.DURATION = TIMES.end;
// sound cues (motion/sound.py): impacts on reveals, whooshes on cuts
window.CUES = (() => {
  const c = [{ t: 0.6, k: "rise", d: 3.2 }, { t: 1.3, k: "hit" }, { t: 4.6, k: "hit" }, { t: 7.6, k: "whoosh" }, { t: 9.7, k: "rise", d: 5.6 },
    { t: 12.6, k: "tick" }, { t: 14.2, k: "tick" }, { t: 16.05, k: "boom" }, { t: 20.5, k: "whoosh" }];
  for (let i = 0; i < 8; i++) c.push({ t: TIMES.tour + i * TIMES.ws, k: "whoosh" }, { t: TIMES.tour + i * TIMES.ws + 0.6, k: "tick" });
  c.push({ t: TIMES.mosaic, k: "whoosh" }, { t: TIMES.mosaic + 0.8, k: "hit" }, { t: TIMES.claim - 0.4, k: "rise", d: 0.4 }, { t: TIMES.claim + 0.15, k: "boom" },
    { t: TIMES.claim + 2.4, k: "hit" }, { t: TIMES.count + 0.3, k: "whoosh" }, { t: TIMES.count + 2.9, k: "rise", d: 0.9 }, { t: TIMES.count + 3.8, k: "hit" },
    { t: TIMES.loops, k: "whoosh" }, { t: TIMES.loops + 0.9, k: "tick" }, { t: TIMES.mol, k: "whoosh" }, { t: TIMES.mol + 2.9, k: "hit" },
    { t: TIMES.mit, k: "whoosh" }, { t: TIMES.mit + 1.9, k: "hit" }, { t: TIMES.fail, k: "whoosh" }, { t: TIMES.fail + 2.1, k: "low" },
    { t: TIMES.burst - 0.8, k: "rise", d: 0.8 });
  for (let i = 0; i < BURST.length; i++) c.push({ t: TIMES.burst + i * TIMES.cut, k: i === BURST.length - 1 ? "boom" : "hit" });
  c.push({ t: TIMES.outro, k: "boom" }, { t: TIMES.outro + 1.0, k: "rise", d: 2.5 }, { t: TIMES.outro + 3.5, k: "hit" });
  return c.sort((a, b) => a.t - b.t);
})();
window.READY = (async () => {
  await document.fonts.ready;
  const imgs = [...document.querySelectorAll("img")];
  await Promise.all(imgs.map((i) => (i.complete ? Promise.resolve() : new Promise((ok) => { i.onload = ok; i.onerror = ok; })).then(() => i.decode?.().catch(() => {}))));
  seek(0);
  return true;
})();

// preview: index.html?t=12.5 shows one frame; ?play plays in real time
(async () => {
  const q = new URLSearchParams(location.search);
  await window.READY;
  if (q.has("bare")) { ui.style.display = "none"; fx.style.display = "none"; }   // the 3D layer alone (stills for the app)
  if (q.has("t")) seek(parseFloat(q.get("t")));
  if (q.has("play")) {
    const t0 = performance.now() - 1000 * parseFloat(q.get("play") || "0");
    const loop = () => { const t = (performance.now() - t0) / 1000; seek(t % TIMES.end); requestAnimationFrame(loop); };
    loop();
  }
})();
