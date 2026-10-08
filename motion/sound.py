"""
The film's soundtrack, synthesised from scratch (no samples, no licences): a 120 bpm A-minor bed whose beat grid
lines up with the feature burst's 0.5 s cuts, plus sound effects on the film's own cue list (window.CUES:
whooshes on transitions, hits and booms on reveals, risers before drops).

    python motion/sound.py motion/out/cues.json motion/out/soundtrack.wav
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from scipy import signal
from scipy.io import wavfile

SR = 48_000
BEAT = 0.5                      # 120 bpm
GRID0 = 0.3                     # beats fall on 0.3 + 0.5 k: the burst cuts (77.8 s + 0.5 k) land on beats
RNG = np.random.default_rng(5)

# sections: (start, end, intensity of the groove 0-1, drums on, arpeggio on)
SECTIONS = [(0.0, 7.5, 0.15, False, False), (7.5, 16.0, 0.45, False, True), (16.0, 21.0, 0.55, False, False),
            (21.0, 48.2, 0.9, True, True), (48.2, 52.2, 1.0, True, True), (52.2, 56.6, 0.3, False, False),
            (56.6, 74.8, 0.7, True, True), (74.8, 77.8, 0.45, False, False), (77.8, 84.4, 1.0, True, True),
            (84.4, 92.0, 0.35, False, True)]
CHORDS = [("A", [57, 60, 64]), ("F", [53, 57, 60]), ("C", [48, 52, 55]), ("G", [55, 59, 62])]   # MIDI, 2 s each


def midi(m: float) -> float:
    return 440.0 * 2 ** ((m - 69) / 12)


def env_ad(n: int, a: float, d: float) -> np.ndarray:
    t = np.arange(n) / SR
    return np.minimum(1, t / max(a, 1e-4)) * np.exp(-np.maximum(t - a, 0) / d)


def add(buf: np.ndarray, x: np.ndarray, t0: float, gain: float = 1.0, pan: float = 0.0) -> None:
    i0 = int(t0 * SR)
    if i0 >= buf.shape[1] or i0 + len(x) <= 0:
        return
    if i0 < 0:
        x, i0 = x[-i0:], 0
    x = x[: buf.shape[1] - i0]
    gl, gr = np.cos((pan + 1) * np.pi / 4), np.sin((pan + 1) * np.pi / 4)
    buf[0, i0:i0 + len(x)] += gain * gl * x * np.sqrt(2)
    buf[1, i0:i0 + len(x)] += gain * gr * x * np.sqrt(2)


def lp(x, f, order=2):
    return signal.sosfilt(signal.butter(order, f, "low", fs=SR, output="sos"), x)


def hp(x, f, order=2):
    return signal.sosfilt(signal.butter(order, f, "high", fs=SR, output="sos"), x)


def bp(x, lo, hi, order=2):
    return signal.sosfilt(signal.butter(order, [lo, hi], "band", fs=SR, output="sos"), x)


# ------------------------------------------------------------------------------------------------ instruments
def kick(level=1.0) -> np.ndarray:
    n = int(0.45 * SR)
    t = np.arange(n) / SR
    f = 42 + 120 * np.exp(-t * 30)
    ph = 2 * np.pi * np.cumsum(f) / SR
    x = np.sin(ph) * np.exp(-t * 7.5) + 0.25 * np.exp(-t * 200) * RNG.standard_normal(n)
    return level * np.tanh(1.6 * x)


def hat(level=1.0) -> np.ndarray:
    n = int(0.09 * SR)
    return level * hp(RNG.standard_normal(n), 7000) * np.exp(-np.arange(n) / SR / 0.025)


def snare(level=1.0) -> np.ndarray:
    n = int(0.35 * SR)
    t = np.arange(n) / SR
    body = np.sin(2 * np.pi * 190 * t) * np.exp(-t * 25)
    nz = bp(RNG.standard_normal(n), 1200, 7000) * np.exp(-t * 14)
    return level * (0.5 * body + 0.8 * nz)


def pluck(freq, dur=0.32, level=1.0) -> np.ndarray:
    n = int(dur * SR)
    t = np.arange(n) / SR
    x = np.sin(2 * np.pi * freq * t) + 0.35 * np.sin(4 * np.pi * freq * t) + 0.12 * signal.sawtooth(2 * np.pi * freq * t)
    return level * x * env_ad(n, 0.003, 0.11)


def whoosh(dur=0.7, up=True) -> np.ndarray:
    n = int(dur * SR)
    nz = RNG.standard_normal(n)
    bands = [(150, 400), (300, 900), (700, 2000), (1500, 4500), (3000, 9000)]
    layers = [bp(nz, lo, hi) for lo, hi in bands]
    pos = np.linspace(0, len(bands) - 1, n) if up else np.linspace(len(bands) - 1, 0, n)
    out = np.zeros(n)
    for k, l in enumerate(layers):
        out += l * np.clip(1 - np.abs(pos - k), 0, 1)
    shape = np.sin(np.pi * np.linspace(0, 1, n)) ** 1.5
    return out * shape


def rise(dur: float) -> np.ndarray:
    n = int(dur * SR)
    t = np.arange(n) / SR
    f = 180 * (6 ** (t / dur))
    tone = np.sin(2 * np.pi * np.cumsum(f) / SR) * 0.35
    nz = whoosh(dur, True) * 0.8
    return (tone + nz) * (t / dur) ** 2


def boom() -> np.ndarray:
    n = int(2.6 * SR)
    t = np.arange(n) / SR
    f = 26 + 70 * np.exp(-t * 3.5)
    sub = np.sin(2 * np.pi * np.cumsum(f) / SR) * np.exp(-t * 1.6)
    crack = lp(RNG.standard_normal(n), 1800) * np.exp(-t * 6)
    return np.tanh(1.4 * (sub + 0.5 * crack))


def hit() -> np.ndarray:
    k = kick(0.9)
    s = snare(0.7)
    out = np.zeros(max(len(k), len(s)))
    out[: len(k)] += k
    out[: len(s)] += s
    return out


def tick() -> np.ndarray:
    n = int(0.05 * SR)
    t = np.arange(n) / SR
    return np.sin(2 * np.pi * 2600 * t) * np.exp(-t * 90)


def sting() -> np.ndarray:      # a dark A-minor stab for the honest failure
    n = int(2.0 * SR)
    t = np.arange(n) / SR
    x = sum(signal.sawtooth(2 * np.pi * midi(m) * t * d) for m in (45, 48, 52) for d in (0.997, 1.003))
    return lp(x, 900) * env_ad(n, 0.01, 0.6) * 0.25


def reverb(x: np.ndarray, secs=2.2, seed=1) -> np.ndarray:
    r = np.random.default_rng(seed)
    n = int(secs * SR)
    ir = r.standard_normal(n) * np.exp(-np.arange(n) / SR / (secs / 6.0))
    ir = lp(ir, 6000)
    ir /= np.sqrt((ir ** 2).sum())
    return signal.fftconvolve(x, ir)[: len(x)]


# ------------------------------------------------------------------------------------------------ score
def intensity_at(t: np.ndarray) -> np.ndarray:
    out = np.zeros_like(t)
    for a, b, lev, _, _ in SECTIONS:
        m = (t >= a) & (t < b)
        out[m] = lev
    from scipy.ndimage import uniform_filter1d
    return uniform_filter1d(out, int(0.4 * SR), mode="nearest")      # smooth the steps a little


def section(t: float):
    for s in SECTIONS:
        if s[0] <= t < s[1]:
            return s
    return SECTIONS[-1]


def synth(cues_path: Path, out: Path) -> Path:
    cues = json.loads(Path(cues_path).read_text(encoding="utf-8"))
    dur = float(cues["duration"]) + 0.5
    N = int(dur * SR)
    tt = np.arange(N) / SR
    music = np.zeros((2, N))
    fxbus = np.zeros((2, N))

    # pad: the chord loop, detuned saws, dark and bright versions mixed by intensity (a filter sweep)
    pad = np.zeros(N)
    for i0 in range(0, N, int(2 * SR)):
        name, notes = CHORDS[(i0 // int(2 * SR)) % 4]
        n = min(int(2.4 * SR), N - i0)
        t = np.arange(n) / SR
        seg = sum(signal.sawtooth(2 * np.pi * midi(m - 12) * t * d) for m in notes for d in (0.996, 1.0, 1.004))
        seg += 0.6 * np.sin(2 * np.pi * midi(notes[0] - 24) * t)
        e = np.minimum(1, t / 0.35) * np.minimum(1, np.maximum(0, (2.4 - t) / 0.4))
        pad[i0:i0 + n] += seg * e
    inten = intensity_at(tt)
    pad = (lp(pad, 500) * (1 - inten) + lp(pad, 2600) * inten) * (0.08 + 0.05 * inten)
    music[0] += pad
    music[1] += np.roll(pad, int(0.012 * SR))

    # beats: kick on the beat, hats on the off-beat, snare on 2 and 4, arpeggio in 16ths
    k_ = kick(0.9)
    h_ = hat(0.25)
    s_ = snare(0.35)
    beat = 0
    tb = GRID0
    while tb < dur - 0.5:
        a, b, lev, drums, arp = section(tb)
        if drums:
            add(music, k_, tb, 0.55 * lev)
            add(music, h_, tb + BEAT / 2, 0.8 * lev, pan=0.35)
            if beat % 2 == 1:
                add(music, s_, tb, 0.6 * lev)
        if arp:
            _, notes = CHORDS[int(tb // 2) % 4]
            seq = [notes[0], notes[1], notes[2], notes[1] + 12, notes[2] + 12, notes[0] + 12, notes[2], notes[1]]
            for j in range(4):
                m = seq[(beat * 4 + j) % len(seq)] + 12
                add(music, pluck(midi(m), level=0.12 * (0.5 + lev)), tb + j * BEAT / 4, pan=0.6 * np.sin(beat + j))
        beat += 1
        tb += BEAT
    # ping-pong echo on the music bus (dotted eighth)
    d = int(0.375 * SR)
    echo = np.zeros_like(music)
    echo[0, d:] += 0.22 * music[1, :-d]
    echo[1, d:] += 0.22 * music[0, :-d]
    music += echo

    # sound effects from the film's cue list
    for c in cues["cues"]:
        k, t0 = c["k"], float(c["t"])
        if k == "whoosh":
            add(fxbus, whoosh(0.7), t0 - 0.35, 0.33, pan=-0.2)
            add(fxbus, whoosh(0.7, up=False), t0 - 0.3, 0.2, pan=0.25)
        elif k == "hit":
            add(fxbus, hit(), t0, 0.55)
        elif k == "boom":
            add(fxbus, boom(), t0, 0.8)
        elif k == "rise":
            add(fxbus, rise(float(c.get("d", 1.0))), t0, 0.32)
        elif k == "tick":
            add(fxbus, tick(), t0, 0.12, pan=0.4)
        elif k == "low":
            add(fxbus, sting(), t0, 1.0)

    wet = np.stack([reverb(fxbus[0] + 0.35 * music[0], seed=1), reverb(fxbus[1] + 0.35 * music[1], seed=2)])
    mixd = music + fxbus + 0.32 * wet
    mixd = hp(mixd, 25)
    # fades and a soft limiter
    fade = np.minimum(1, tt / 0.8) * np.clip((dur - 0.5 - tt) / 1.6, 0, 1)
    mixd *= fade
    mixd = np.tanh(1.3 * mixd / (np.abs(mixd).max() + 1e-9)) / np.tanh(1.3)
    mixd *= 10 ** (-1.0 / 20)
    out = Path(out)
    out.parent.mkdir(parents=True, exist_ok=True)
    wavfile.write(out, SR, (mixd.T * 32767).astype(np.int16))
    print(f"soundtrack: {out} ({dur:.1f} s, {len(cues['cues'])} cues)")
    return out


if __name__ == "__main__":
    here = Path(__file__).resolve().parent
    synth(Path(sys.argv[1]) if len(sys.argv) > 1 else here / "out" / "cues.json",
          Path(sys.argv[2]) if len(sys.argv) > 2 else here / "out" / "soundtrack.wav")
