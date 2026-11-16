"""
Demo: how EpiSim would use the model to compare response scenarios.

Run from the project folder:   python3 scripts/demo_scenarios.py

It loads the parameter set calibrated to Kikwit 1995, applies it to a district of
200,000 people with 5 first cases, and compares three response timings.
Writes results/demo_scenarios.json and results/fig_demo_scenarios.png.

These are scenario estimates, not forecasts.
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episim_ebola import EbolaParams, Response, simulate_stochastic, summarise  # noqa: E402

preset = json.load(open(ROOT / "calibrated" / "kikwit_1995.json"))
params = EbolaParams(**preset["parameters"]).with_changes(reporting_fraction=1.0)
obs = preset["response_as_observed"]
DAYS, RUNS, SEED, FIRST_CASES = 365, 500, 2026, 5


def response(start_day):
    """The response seen in Kikwit, started on a day of our choosing."""
    return (Response(start_day=start_day, ramp_days=obs["ramp_days"],
                     community=obs["community"], hospital=obs["hospital"],
                     funeral=obs["funeral"], label=f"Response from day {start_day}"),)


scenarios = {
    "Response from day 60": response(60),
    "Response from day 90": response(90),
    "Response from day 120": response(120),
}

output = {"preset": preset["name"], "status": preset["status"], "days": DAYS, "runs": RUNS,
          "seed": SEED, "first_cases": FIRST_CASES, "population": params.population,
          "scenarios": {}}
print(f"{'Scenario':<24}{'Cases (median, 5-95%)':>30}{'Deaths (median)':>18}{'Peak beds (median, 5-95%)':>30}")
for name, resp in scenarios.items():
    result = simulate_stochastic(params, DAYS, resp, runs=RUNS, seed=SEED,
                                 initial_infectious=FIRST_CASES)
    s = summarise(result)
    output["scenarios"][name] = s
    c, d, b = (s["totals_all_runs"][k] for k in ("cases", "deaths", "peak_in_hospital"))
    print(f"{name:<24}{c['median']:>10,.0f} ({c['p05']:,.0f} to {c['p95']:,.0f})".ljust(54)
          + f"{d['median']:>12,.0f}" + f"{b['median']:>14,.0f} ({b['p05']:,.0f} to {b['p95']:,.0f})")

json.dump(output, open(ROOT / "results" / "demo_scenarios.json", "w"))

import matplotlib   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
for (name, s), colour in zip(output["scenarios"].items(), ("#1e8449", "#2471a3", "#c0392b")):
    for a, key, title in ((ax[0], "cumulative_cases", "Cumulative cases"),
                          (ax[1], "in_hospital", "People in hospital (beds needed)")):
        band = s["daily"][key]
        a.plot(band["median"], color=colour, lw=2, label=name)
        a.fill_between(range(DAYS + 1), band["p05"], band["p95"], color=colour, alpha=0.13)
        a.set_title(title + ": median and 5-95% range of 500 runs", fontsize=10)
        a.set_xlabel("Day")
ax[0].legend(fontsize=9)
fig.suptitle("Scenario estimates, not forecasts. District of 200,000, 5 first cases, "
             "parameters calibrated to Kikwit 1995", fontsize=9)
plt.tight_layout()
plt.savefig(ROOT / "results" / "fig_demo_scenarios.png", dpi=150)
print("\nSaved results/demo_scenarios.json and results/fig_demo_scenarios.png")
