"""
Side check: how much does the West Africa R0 estimate depend on the fitting window?
Run from the project folder:  python3 scripts/west_africa_window_check.py
Writes results/west_africa_windows.json
"""
import json, sys
from pathlib import Path
import numpy as np, pandas as pd
from scipy.optimize import minimize
ROOT = Path(__file__).resolve().parents[1]; sys.path.insert(0, str(ROOT))
from episim_ebola import EbolaParams, simulate_deterministic

wa = pd.read_csv(ROOT / "data" / "west_africa_2014_cumulative.csv", parse_dates=["date"])
BASE = EbolaParams(onset_to_admission=5.0, fatality_home=0.708, fatality_hospital=0.708)
BURN = 60
POP = {"Guinea": 11_800_000, "Liberia": 4_400_000, "Sierra Leone": 7_100_000}

def fit(country, start, end):
    s = wa[(wa.country == country) & (wa.date >= start) & (wa.date <= end)].copy()
    t = (s.date - pd.Timestamp(start)).dt.days.values
    y = s.cumulative_cases.values; base = BASE.with_changes(population=POP[country])
    def loss(th):
        r0, i0 = np.exp(th)
        if not (0.5 < r0 < 6 and 1e-3 < i0 < 5000): return 1e9
        d = simulate_deterministic(base.with_changes(r0=r0), BURN + int(t.max()), initial_infectious=i0, steps_per_day=4)
        cum = d["onsets"].cumsum(); m = y[0] + cum[BURN:][t] - cum[BURN]
        return float(np.sum((np.log(m) - np.log(y)) ** 2))
    res = min((minimize(loss, np.log([r, i]), method="Nelder-Mead", options={"xatol":1e-5,"fatol":1e-9,"maxiter":3000})
               for r in (1.5, 2.2) for i in (5.0, 50.0)), key=lambda m: m.fun)
    return float(np.exp(res.x[0])), int(len(s))

out = []
for country, windows in {
    "Guinea": [("2014-07-20", "2014-09-14"), ("2014-07-20", "2014-08-20"), ("2014-08-01", "2014-09-14")],
    "Liberia": [("2014-06-20", "2014-09-14"), ("2014-06-20", "2014-08-01"), ("2014-07-15", "2014-09-14")],
    "Sierra Leone": [("2014-06-01", "2014-09-14"), ("2014-06-01", "2014-07-15"), ("2014-07-15", "2014-09-14")],
}.items():
    for start, end in windows:
        r0, n = fit(country, start, end)
        out.append({"country": country, "start": start, "end": end, "points": n, "r0": r0})
        print(f"{country:13s} {start} to {end}  ({n:2d} points)  R0 = {r0:.2f}")
json.dump(out, open(ROOT / "results" / "west_africa_windows.json", "w"), indent=1)
