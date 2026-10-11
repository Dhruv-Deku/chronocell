"""
Gate 3b (Phase A4): a learned correction on top of the population model, trained on practice data only.

chronocell/learned_correction.py computes per-pair features from the input and the fitted model and
applies a frozen correction r = log(true / model median). Here: candidates are compared by
leave-one-dataset-out over the practice units of validation/phase_a_units.py (imaging and Hi-C input,
every depth), the chosen one is refitted on all practice units and frozen in
chronocell/data/learned_correction.json. The test is the Gate 3 benchmark plan with the new method
`learned_correction` (validation/benchmark/methods.py), run in parts and scored against every method
already in validation/benchmark/results.json (python validation/learned_correction.py --test).

Candidates: none (the model as is); linear (ridge on the standardised features); mlp (two hidden layers of
64 SiLU units, trained on the GPU with AdamW). Every unit weighs the same in training.
Selection: highest mean held-out trend-removed Spearman (the pattern beyond the separation trend, the
numerator of "% of the reproducible pattern"); a candidate must beat "none" on every held-out group to
be chosen, else "none" is frozen (no correction).
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent
sys.path.insert(0, str(ROOT.parent))
sys.path.insert(0, str(ROOT))

import phase_a as A             # noqa: E402
import protocol as PR           # noqa: E402
from chronocell import learned_correction as LC     # noqa: E402


def unit_data(spec) -> dict:
    u = A.load(spec)
    n = len(u["median"])
    meta = u["meta"]
    f = LC.features(u["median"].astype(np.float64), u["freq"], u["counts"], u["sep"], meta["r_c_nm"],
                    meta["protocol"] != "imaging", meta["n_eff"], meta["step_bp"])
    iu = np.triu_indices(n, 1)
    m, t = u["median"][iu].astype(np.float64), u["truth"][iu].astype(np.float64)
    ok = np.isfinite(m) & np.isfinite(t) & (m > 0) & (t > 0) & np.isfinite(f[iu][:, 1])
    return {"X": f[iu][ok], "y": np.log(t[ok] / m[ok]), "model": m[ok], "truth": t[ok], "sep": u["sep"][iu][ok].astype(float)}


def _standardiser(X: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Feature means and SDs ignoring missing values (a pair's input can be unobserved)."""
    mu, sd = np.nanmean(X, axis=0), np.nanstd(X, axis=0)
    return np.nan_to_num(mu), np.where(np.isfinite(sd) & (sd > 0), sd, 1.0)


def _std(X: np.ndarray, mu: np.ndarray, sd: np.ndarray) -> np.ndarray:
    """Standardised features; a missing value becomes 0 (the mean), as in chronocell.learned_correction."""
    z = (X - mu) / sd
    return np.where(np.isfinite(z), z, 0.0)


def fit_linear(units: list[dict], lam: float = 1.0) -> dict:
    X = np.concatenate([u["X"] for u in units])
    y = np.concatenate([u["y"] for u in units])
    w = np.concatenate([np.full(len(u["y"]), 1.0 / len(u["y"])) for u in units])
    mu, sd = _standardiser(X)
    Z = np.column_stack([np.ones(len(X)), _std(X, mu, sd)])
    A_ = Z.T @ (w[:, None] * Z) + lam * np.diag([0.0] + [1.0] * X.shape[1]) * w.sum() / len(units)
    beta = np.linalg.solve(A_, Z.T @ (w * y))
    return {"kind": "linear", "mu": mu.tolist(), "sd": sd.tolist(), "intercept": float(beta[0]), "beta": beta[1:].tolist()}


def fit_mlp(units: list[dict], epochs: int = 15, seed: int = 0) -> dict:
    import torch
    torch.manual_seed(seed)
    dev = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    X = np.concatenate([u["X"] for u in units])
    y = np.concatenate([u["y"] for u in units])
    w = np.concatenate([np.full(len(u["y"]), 1.0 / len(u["y"])) for u in units])
    mu, sd = _standardiser(X)
    Xt = torch.as_tensor(_std(X, mu, sd), dtype=torch.float32, device=dev)
    yt = torch.as_tensor(y, dtype=torch.float32, device=dev)
    wt = torch.as_tensor(w / w.mean(), dtype=torch.float32, device=dev)
    net = torch.nn.Sequential(torch.nn.Linear(X.shape[1], 64), torch.nn.SiLU(), torch.nn.Linear(64, 64), torch.nn.SiLU(),
                              torch.nn.Linear(64, 1)).to(dev)
    opt = torch.optim.AdamW(net.parameters(), lr=1e-3, weight_decay=1e-4)
    g = torch.Generator(device="cpu").manual_seed(seed)
    for _ in range(epochs):
        perm = torch.randperm(len(yt), generator=g).to(dev)
        for k in range(0, len(yt), 8192):
            idx = perm[k:k + 8192]
            loss = (wt[idx] * (net(Xt[idx])[:, 0] - yt[idx]) ** 2).mean()
            opt.zero_grad()
            loss.backward()
            opt.step()
    layers = [m for m in net if isinstance(m, torch.nn.Linear)]
    return {"kind": "mlp", "mu": mu.tolist(), "sd": sd.tolist(),
            "layers": [{"W": l.weight.detach().cpu().numpy().tolist(), "b": l.bias.detach().cpu().numpy().tolist()} for l in layers],
            "device_trained": str(dev)}


def predict_r(params: dict | None, X: np.ndarray) -> np.ndarray:
    if params is None:
        return np.zeros(len(X))
    if params["kind"] == "linear":
        z = (X - np.asarray(params["mu"])) / np.asarray(params["sd"])
        return np.where(np.isfinite(z), z, 0.0) @ np.asarray(params["beta"]) + params["intercept"]
    return LC._mlp(params, X)


def evaluate(units: list[dict], params) -> dict:
    rho, ccc = [], []
    for u in units:
        pred = u["model"] * np.exp(predict_r(params, u["X"]))
        s = PR._flat_scores(pred, u["truth"], u["sep"])
        rho.append(s["spearman_distance_corrected"])
        ccc.append(s["lin_ccc_nm"])
    return {"pattern_rho": float(np.nanmean(rho)), "ccc": float(np.nanmean(ccc))}


def practice() -> None:
    import phase_a_units as U
    from hic_size_calibration import group_of
    specs = U.PRACTICE
    data = {s: unit_data(s) for s in specs}
    groups = sorted({group_of(s) for s in specs})
    table = {}
    for cand in ("none", "linear", "mlp"):
        per = {}
        for g in groups:
            train = [d for s, d in data.items() if group_of(s) != g]
            held = [d for s, d in data.items() if group_of(s) == g]
            params = None if cand == "none" else (fit_linear(train) if cand == "linear" else fit_mlp(train))
            per[g] = evaluate(held, params)
        table[cand] = {"per_group": per, "mean_pattern_rho": float(np.mean([v["pattern_rho"] for v in per.values()])),
                       "mean_ccc": float(np.mean([v["ccc"] for v in per.values()]))}
        print(f"{cand:7s} held-out pattern rho {table[cand]['mean_pattern_rho']:+.4f}  CCC {table[cand]['mean_ccc']:.3f}  " +
              "  ".join(f"{g}: {v['pattern_rho']:+.3f}/{v['ccc']:.2f}" for g, v in per.items()), flush=True)
    best = "none"
    for cand in ("linear", "mlp"):
        beats = all(table[cand]["per_group"][g]["pattern_rho"] > table["none"]["per_group"][g]["pattern_rho"] for g in groups)
        if beats and table[cand]["mean_pattern_rho"] > table[best]["mean_pattern_rho"]:
            best = cand
    allu = list(data.values())
    params = None if best == "none" else (fit_linear(allu) if best == "linear" else fit_mlp(allu))
    out = {"chosen": best, "features": list(LC.FEATURES), "params": params, "trained_on": [s.name for s in specs],
           "selection": "highest mean held-out pattern rho; must beat 'none' on every held-out group",
           "fitted_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds"),
           "note": "Practice data only. Held-out test: Gate 3b (validation/RESULTS.md)."}
    if params is not None:
        out.update(params)
    LC.PATH.write_text(json.dumps(out) + "\n", encoding="utf-8")
    (ROOT / "results_learned_correction_practice.json").write_text(json.dumps(
        {"mode": "practice (leave-one-dataset-out)", "table": table, "chosen": best}, indent=1, default=float))
    print(f"chosen: {best}")


PARTS = ("part_lc_bintu", "part_lc_su_chr21", "part_lc_su_chr21_rep", "part_lc_su_genome")


def test() -> None:
    """Gate 3b: the method's benchmark parts (python -m validation.benchmark.run --datasets ... --methods
    learned_correction --out part_lc_*) scored with every method of validation/benchmark/results.json; that file
    is not changed. Rule in frozen.LEARNED_CORRECTION."""
    from frozen import LEARNED_CORRECTION as R
    from validation.benchmark import run as BR
    if R.get("status") == "not run":
        raise SystemExit(f"Gate 3b is not run: {R['reason']} (validation/frozen.py).")
    bench = BR.HERE
    base = json.loads((bench / "results.json").read_text())
    rows = [r for r in base["rows"] if r["method"] != "learned_correction"]
    for p in PARTS:
        rows += json.loads((bench / f"{p}.json").read_text())["rows"]
    summ = BR.summarise(rows)
    units = []
    for ds, inputs in summ.items():
        for inp, msets in inputs.items():
            ap = msets.get("all_pairs", {})
            if "learned_correction" not in ap:
                continue
            learned = ap["learned_correction"]["percent_of_ceiling"]
            others = {m: v["percent_of_ceiling"] for m, v in ap.items() if m != "learned_correction"}
            best_name = max(others, key=others.get)
            current = ap["v4_whole"]["percent_of_ceiling"] if "v4_whole" in ap else ap["v3_3_windowed"]["percent_of_ceiling"]
            ok = learned >= others[best_name] if inp == "hic" else learned >= current - R["imaging_max_loss_points"]
            units.append({"unit": ds, "input": inp, "learned": learned, "best_other": others[best_name],
                          "best_other_name": best_name, "current": current, "ok": bool(ok)})
    hic = [u for u in units if u["input"] == "hic"]
    img = [u for u in units if u["input"] == "imaging"]
    frac = sum(u["ok"] for u in hic) / max(len(hic), 1)
    verdict = "pass" if frac >= R["hic_best_fraction"] and all(u["ok"] for u in img) else "fail"
    line = (f"Test (run once): learned correction best on {sum(u['ok'] for u in hic)} of {len(hic)} Hi-C units "
            f"(needed: ≥ {100 * R['hic_best_fraction']:.0f} %); within {R['imaging_max_loss_points']:g} point of the current "
            f"model on {sum(u['ok'] for u in img)} of {len(img)} imaging units (needed: all).")
    (ROOT / "results_gate3b.json").write_text(json.dumps(
        {"mode": "test (run once; correction frozen)", "rule": R, "units": units, "summary_line": line, "verdict": verdict,
         "parts": list(PARTS), "summary": summ, "run_utc": dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")},
        indent=1, default=float))
    print(line, "->", verdict)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--practice", action="store_true")
    g.add_argument("--test", action="store_true")
    a = ap.parse_args()
    practice() if a.practice else test()


if __name__ == "__main__":
    main()
