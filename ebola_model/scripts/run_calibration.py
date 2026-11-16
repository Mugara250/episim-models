"""
Calibrate the EpiSim Ebola model and test it against real outbreak data.

Run from the project folder:   python3 scripts/run_calibration.py

Writes  results/results.json            every number quoted in MODEL.md
        results/fig_*.png               the figures
        calibrated/kikwit_1995.json     the calibrated parameter set for the demo
"""
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episim_ebola import (EbolaParams, Response, simulate_deterministic,   # noqa: E402
                          simulate_stochastic, summarise)
from episim_ebola.calibrate import (weekly_sums, fit_outbreak, profile_interval,  # noqa: E402
                                    predict_weekly, negbin_interval)

RESULTS = {}
RAMP_GRID = (0, 7, 14, 21)

# =============================================================================
# PART 1. KIKWIT, DEMOCRATIC REPUBLIC OF THE CONGO, 1995
# =============================================================================
# Data: Khan et al. (1999), daily onsets and deaths, day 0 = 6 January 1995.
kik = pd.read_csv(ROOT / "data" / "kikwit_1995_daily.csv")
FIRST_DAY = 52        # 27 February 1995; weekly bins end on the last day of data
WEEKS = 20
RESPONSE_DAY = 118    # 4 May 1995: first interventions (Legrand et al. 2007)
obs_onsets = weekly_sums(kik.onsets.values, FIRST_DAY, WEEKS)
obs_deaths = weekly_sums(kik.deaths.values, FIRST_DAY, WEEKS)

# Fixed values for this outbreak, each with its source (see MODEL.md, section 4)
KIKWIT = EbolaParams(
    population=200_000,                 # Kikwit area, Legrand et al. 2007
    share_community=0.5 / 2.7,          # Legrand's central estimate of R0 by route
    share_hospital=0.4 / 2.7,
    share_funeral=1.8 / 2.7,
    incubation=8.5, onset_to_death=9.3, onset_to_recovery=13.0,   # Nash et al. 2024
    burial_delay=2.0,                   # Legrand et al. 2007
    onset_to_admission=5.0,             # Legrand et al. 2007 (DRC 1995)
    share_hospitalised=0.80,            # Legrand et al. 2007
    fatality_home=0.81, fatality_hospital=0.81,   # 81% overall, Khan et al. 1999
    reporting_fraction=292 / 315,       # cases with a known onset date
)


def fit_best_ramp(base, obs, fit_weeks=None, response_day=RESPONSE_DAY, fixed=None):
    """Fit for each ramp length in the grid and keep the most likely one."""
    fits = [fit_outbreak(base, obs, FIRST_DAY, response_day, fit_weeks=fit_weeks,
                         fixed=fixed, ramp_days=r) for r in RAMP_GRID]
    return max(fits, key=lambda f: f["loglik"]), fits


def predict(base, f, what="onsets", response_day=RESPONSE_DAY):
    return predict_weekly(base, f["r0"], f["remaining"], f["i0"], response_day,
                          FIRST_DAY, WEEKS, what=what, ramp_days=f["ramp_days"])


print("PART 1: Kikwit 1995")
best, by_ramp = fit_best_ramp(KIKWIT, obs_onsets)
r0_ci = profile_interval(KIKWIT, obs_onsets, FIRST_DAY, RESPONSE_DAY, best, "r0", 1.2, 6.0)
rem_ci = profile_interval(KIKWIT, obs_onsets, FIRST_DAY, RESPONSE_DAY, best,
                          "remaining", 0.005, 0.6)
mu_onsets = predict(KIKWIT, best)
lo, hi = negbin_interval(mu_onsets, best["size"])
inside = (obs_onsets >= lo) & (obs_onsets <= hi)
print(f"  R0 = {best['r0']:.2f} ({r0_ci[0]:.2f} to {r0_ci[1]:.2f}), "
      f"ramp {best['ramp_days']:.0f} days, R after response = {best['r_after']:.2f}")

RESULTS["kikwit_fit"] = {
    "r0": best["r0"], "r0_ci95": r0_ci,
    "transmission_remaining": best["remaining"], "transmission_remaining_ci95": rem_ci,
    "r_after_response": best["r_after"],
    "r_after_response_ci95": [best["r0"] * rem_ci[0], best["r0"] * rem_ci[1]],
    "ramp_days": best["ramp_days"], "initial_infectious": best["i0"],
    "scatter_size": best["size"], "loglik": best["loglik"],
    "observed_weekly_onsets": obs_onsets.tolist(),
    "model_weekly_onsets": mu_onsets.tolist(),
    "observed_total": float(obs_onsets.sum()), "model_total": float(mu_onsets.sum()),
    "weeks_inside_90pct_range": int(inside.sum()), "weeks": WEEKS,
    "mean_abs_error_per_week": float(np.abs(obs_onsets - mu_onsets).mean()),
    "generation_time_days": KIKWIT.generation_time,
}

# ---- 1b. Sensitivity: how much does R0 move under other reasonable choices? ----
print("  sensitivity analysis...")
sens = []


def add(label, base=KIKWIT, response_day=RESPONSE_DAY, fixed=None, ramps=None):
    if ramps is None:
        f, _ = fit_best_ramp(base, obs_onsets, response_day=response_day, fixed=fixed)
    else:
        f = fit_outbreak(base, obs_onsets, FIRST_DAY, response_day, fixed=fixed,
                         ramp_days=ramps)
    sens.append({"case": label, "r0": f["r0"], "r_after_response": f["r_after"],
                 "ramp_days": f["ramp_days"], "loglik": f["loglik"]})


add("Main calibration")
for r, f in zip(RAMP_GRID, by_ramp):
    sens.append({"case": f"Response builds up over {r} days", "r0": f["r0"],
                 "r_after_response": f["r_after"], "ramp_days": f["ramp_days"],
                 "loglik": f["loglik"]})
add("Start fixed at 1 case on 6 January", fixed={"i0": 1.0})
add("Response dated 9 May instead of 4 May", response_day=123)
add("Route shares 55/15/30 (generic preset)",
    base=KIKWIT.with_changes(share_community=0.55, share_hospital=0.15, share_funeral=0.30))
add("All transmission in the community",
    base=KIKWIT.with_changes(share_community=1.0, share_hospital=0.0, share_funeral=0.0))
add("Survivors infectious for 10 days (Legrand)",
    base=KIKWIT.with_changes(onset_to_recovery=10.0))
add("Incubation 6 days", base=KIKWIT.with_changes(incubation=6.0))
add("Incubation 11 days", base=KIKWIT.with_changes(incubation=11.0))
RESULTS["kikwit_sensitivity"] = sens
r0s = [s["r0"] for s in sens]
RESULTS["kikwit_sensitivity_r0_range"] = [min(r0s), max(r0s)]
print(f"  R0 across all cases: {min(r0s):.2f} to {max(r0s):.2f}")

# ---- 1c. Test on deaths (never used in the fit) --------------------------------
mu_deaths = predict(KIKWIT, best, what="deaths")
peak_obs, peak_mod = int(obs_deaths.argmax()), int(mu_deaths.argmax())
RESULTS["kikwit_deaths_test"] = {
    "observed_weekly_deaths": obs_deaths.tolist(),
    "model_weekly_deaths": mu_deaths.tolist(),
    "observed_total": float(obs_deaths.sum()), "model_total": float(mu_deaths.sum()),
    "observed_peak_week": peak_obs, "model_peak_week": peak_mod,
    "mean_abs_error_per_week": float(np.abs(obs_deaths - mu_deaths).mean()),
    "correlation": float(np.corrcoef(obs_deaths, mu_deaths)[0, 1]),
}
print(f"  deaths test: observed {obs_deaths.sum():.0f}, model {mu_deaths.sum():.0f}, "
      f"peak week observed {peak_obs}, model {peak_mod}")

# ---- 1d. Hold-out test: fit on the early weeks, predict the rest ---------------
print("  hold-out tests...")
holdout = []
for fit_weeks in (10, 11, 12, 13, 14):
    f, _ = fit_best_ramp(KIKWIT, obs_onsets, fit_weeks=fit_weeks)
    mu = predict(KIKWIT, f)
    last_day_used = FIRST_DAY + 7 * fit_weeks - 1
    holdout.append({
        "weeks_used": fit_weeks,
        "last_date_used": str(kik.date[last_day_used]),
        "days_after_response_start": int(last_day_used - RESPONSE_DAY),
        "r0": f["r0"], "r_after_response": f["r_after"], "ramp_days": f["ramp_days"],
        "observed_remaining_cases": float(obs_onsets[fit_weeks:].sum()),
        "predicted_remaining_cases": float(mu[fit_weeks:].sum()),
        "observed_final_total": float(obs_onsets.sum()),
        "predicted_final_total": float(obs_onsets[:fit_weeks].sum() + mu[fit_weeks:].sum()),
        "predicted_weekly": mu.tolist(),
    })
    h = holdout[-1]
    print(f"    data to {h['last_date_used']}: predicted {h['predicted_remaining_cases']:.0f} "
          f"more cases, observed {h['observed_remaining_cases']:.0f}")
RESULTS["kikwit_holdout"] = holdout

# Pre-response data only: what R0 would we have estimated before any response?
pre = fit_outbreak(KIKWIT, obs_onsets, FIRST_DAY, RESPONSE_DAY, fit_weeks=9,
                   fixed={"remaining": 1.0})
RESULTS["kikwit_pre_response_r0"] = pre["r0"]

# ---- 1e. Stochastic check: do the random runs cover what happened? -------------
START_WEEK = 5                                   # start once cases are established
start_day = FIRST_DAY + 7 * START_WEEK
resp = (Response(start_day=RESPONSE_DAY, community=best["remaining"],
                 hospital=best["remaining"], funeral=best["remaining"],
                 ramp_days=best["ramp_days"]),)
fitted = KIKWIT.with_changes(r0=best["r0"])
det = simulate_deterministic(fitted, FIRST_DAY + 7 * WEEKS, resp,
                             initial_infectious=best["i0"])
state = {k: det[k][start_day] for k in "SEIHFR"}
state["S"] = KIKWIT.population - sum(round(state[k]) for k in "EIHFR")
shift = (Response(start_day=RESPONSE_DAY - start_day, community=best["remaining"],
                  hospital=best["remaining"], funeral=best["remaining"],
                  ramp_days=best["ramp_days"]),)
RUNS = 1000
sto = simulate_stochastic(fitted, 7 * (WEEKS - START_WEEK), shift, runs=RUNS, seed=2026,
                          initial_state=state)
rng = np.random.default_rng(7)
reported = rng.binomial(sto["onsets"], KIKWIT.reporting_fraction)
weekly_runs = reported[1:].reshape(WEEKS - START_WEEK, 7, RUNS).sum(axis=1)
s_lo, s_mid, s_hi = np.percentile(weekly_runs, [5, 50, 95], axis=1)
obs_tail = obs_onsets[START_WEEK:]
cover = (obs_tail >= s_lo) & (obs_tail <= s_hi)
totals = weekly_runs.sum(axis=0)
RESULTS["kikwit_stochastic_check"] = {
    "runs": RUNS, "start_date": str(kik.date[start_day]),
    "start_state": {k: int(round(v)) for k, v in state.items()},
    "weeks_checked": int(len(obs_tail)), "weeks_inside_5_95_range": int(cover.sum()),
    "observed_total_from_start": float(obs_tail.sum()),
    "simulated_total_p05_median_p95": [float(x) for x in np.percentile(totals, [5, 50, 95])],
    "weekly_p05": s_lo.tolist(), "weekly_median": s_mid.tolist(), "weekly_p95": s_hi.tolist(),
}
print(f"  stochastic check: {cover.sum()} of {len(obs_tail)} weeks inside the 5-95% range")

# ---- 1f. A scenario for the demo: what if the response had come later? ---------
scen = {}
for label, delay in (("Response as it happened (4 May)", 0),
                     ("Response 2 weeks later", 14), ("Response 4 weeks later", 28)):
    r = (Response(start_day=RESPONSE_DAY + delay, community=best["remaining"],
                  hospital=best["remaining"], funeral=best["remaining"],
                  ramp_days=best["ramp_days"]),)
    d = simulate_deterministic(fitted, 400, r, initial_infectious=best["i0"])
    scen[label] = {"total_cases": float(d["onsets"].sum()),
                   "total_deaths": float(d["deaths"].sum()),
                   "peak_in_hospital": float(d["H"].max())}
RESULTS["kikwit_delay_scenarios"] = scen

# =============================================================================
# PART 2. WEST AFRICA 2014 (a second, independent outbreak)
# =============================================================================
# Data: WHO situation reports, cumulative reported cases by country.
print("PART 2: West Africa 2014")
wa = pd.read_csv(ROOT / "data" / "west_africa_2014_cumulative.csv", parse_dates=["date"])
FIT_END = pd.Timestamp("2014-09-14")       # same cut-off as WHO Ebola Response Team
COUNTRIES = {   # fit window start, approximate population, WHO ERT estimate of R0
    "Guinea":       (pd.Timestamp("2014-07-20"), 11_800_000, 1.71),
    "Liberia":      (pd.Timestamp("2014-06-20"), 4_400_000, 1.83),
    "Sierra Leone": (pd.Timestamp("2014-06-01"), 7_100_000, 2.02),
}
BURN_IN = 60     # days the model runs before the window, to settle into steady growth
WEST_AFRICA = EbolaParams(onset_to_admission=5.0, fatality_home=0.708,
                          fatality_hospital=0.708)   # WHO Ebola Response Team 2014

from scipy.optimize import minimize   # noqa: E402


def wa_curve(base, r0, i0, days):
    d = simulate_deterministic(base.with_changes(r0=r0), BURN_IN + days,
                               initial_infectious=i0, steps_per_day=4)
    cum = d["onsets"].cumsum()
    return cum[BURN_IN:] - cum[BURN_IN]


west = {}
for country, (start, pop, who_r0) in COUNTRIES.items():
    s = wa[(wa.country == country) & (wa.date >= start)].copy()
    s["t"] = (s.date - start).dt.days
    train = s[s.date <= FIT_END]
    c0 = train.cumulative_cases.iloc[0]
    base = WEST_AFRICA.with_changes(population=pop)
    horizon = int(s.t.max())

    def loss(theta, train=train, base=base, c0=c0):
        r0, i0 = np.exp(theta)
        if not (0.5 < r0 < 6 and 1e-3 < i0 < 5000):
            return 1e9
        model = c0 + wa_curve(base, r0, i0, int(train.t.max()))[train.t.values]
        return float(np.sum((np.log(model) - np.log(train.cumulative_cases.values)) ** 2))

    res = min((minimize(loss, np.log([r, i]), method="Nelder-Mead",
                        options={"xatol": 1e-5, "fatol": 1e-9, "maxiter": 3000})
               for r in (1.5, 2.2) for i in (5.0, 50.0)), key=lambda m: m.fun)
    r0, i0 = np.exp(res.x)
    curve = c0 + wa_curve(base, r0, i0, horizon)
    checks = []
    for weeks_ahead in (0, 4, 8):
        target = FIT_END + pd.Timedelta(weeks=weeks_ahead)
        row = s.iloc[(s.date - target).abs().argsort().iloc[0]]
        checks.append({"weeks_after_cutoff": weeks_ahead, "date": str(row.date.date()),
                       "reported": int(row.cumulative_cases),
                       "model": float(curve[int(row.t)])})
    west[country] = {
        "window_start": str(start.date()), "fit_end": str(FIT_END.date()),
        "points_fitted": int(len(train)), "r0": float(r0), "who_ert_r0": who_r0,
        "rms_log_error": float(np.sqrt(res.fun / len(train))), "checks": checks,
        "dates": [str(d.date()) for d in s.date], "reported": s.cumulative_cases.tolist(),
        "model": [float(curve[int(t)]) for t in s.t],
    }
    c = checks
    print(f"  {country}: R0 {r0:.2f} (WHO team {who_r0}); +4 weeks model {c[1]['model']:.0f} "
          f"vs reported {c[1]['reported']}; +8 weeks model {c[2]['model']:.0f} "
          f"vs reported {c[2]['reported']}")
RESULTS["west_africa"] = west
RESULTS["west_africa_generation_time_days"] = WEST_AFRICA.generation_time

# =============================================================================
# PART 3. SAVE RESULTS, THE CALIBRATED PRESET AND THE FIGURES
# =============================================================================
(ROOT / "results").mkdir(exist_ok=True)
(ROOT / "calibrated").mkdir(exist_ok=True)
json.dump(RESULTS, open(ROOT / "results" / "results.json", "w"), indent=1)

preset = {
    "name": "Ebola, calibrated to Kikwit 1995",
    "status": "Calibrated to one outbreak. Not reviewed by an epidemiologist. "
              "Scenario estimates, not forecasts.",
    "parameters": fitted.to_dict(),
    "response_as_observed": {"start_day": RESPONSE_DAY, "ramp_days": best["ramp_days"],
                             "community": best["remaining"], "hospital": best["remaining"],
                             "funeral": best["remaining"]},
    "initial_infectious": best["i0"],
    "day_zero": "1995-01-06",
    "calibration": {"r0_ci95": r0_ci,
                    "r0_range_across_assumptions": RESULTS["kikwit_sensitivity_r0_range"],
                    "data": "Khan et al. 1999, J Infect Dis 179:S76-S86",
                    "method": "Maximum likelihood on weekly onsets, negative binomial"},
}
json.dump(preset, open(ROOT / "calibrated" / "kikwit_1995.json", "w"), indent=1)

import matplotlib   # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt   # noqa: E402

BLUE, RED, GREY, GREEN = "#2471a3", "#c0392b", "#7f8c8d", "#1e8449"
week_dates = pd.to_datetime(kik.date[[FIRST_DAY + 7 * w for w in range(WEEKS)]].values)
resp_date = pd.Timestamp(kik.date[RESPONSE_DAY])

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
ax[0].bar(week_dates, obs_onsets, width=5.5, color=GREY, alpha=0.55, label="Reported cases")
ax[0].plot(week_dates, mu_onsets, color=BLUE, lw=2, marker="o", ms=4, label="Model (fitted)")
ax[0].fill_between(week_dates, lo, hi, color=BLUE, alpha=0.15, label="90% range")
ax[0].axvline(resp_date, color=RED, ls="--", lw=1); ax[0].text(resp_date, ax[0].get_ylim()[1] * 0.93, " response starts", color=RED, fontsize=8)
ax[0].set_title("Kikwit 1995: new cases per week (used for fitting)")
ax[0].set_ylabel("Cases per week"); ax[0].legend(fontsize=8, loc="upper left")
ax[1].bar(week_dates, obs_deaths, width=5.5, color=GREY, alpha=0.55, label="Reported deaths")
ax[1].plot(week_dates, mu_deaths, color=GREEN, lw=2, marker="o", ms=4, label="Model (not fitted to deaths)")
ax[1].axvline(resp_date, color=RED, ls="--", lw=1)
ax[1].set_title("Kikwit 1995: deaths per week (test)")
ax[1].set_ylabel("Deaths per week"); ax[1].legend(fontsize=8, loc="upper left")
for a in ax:
    a.tick_params(axis="x", rotation=30)
plt.tight_layout(); plt.savefig(ROOT / "results" / "fig_kikwit_fit.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 2, figsize=(12.5, 4.4))
ax[0].bar(week_dates, obs_onsets, width=5.5, color=GREY, alpha=0.55, label="Reported cases")
for h, colour in zip(holdout, ("#d35400", "#e67e22", "#8e44ad", "#2471a3", "#1e8449")):
    n = h["weeks_used"]
    ax[0].plot(week_dates[n - 1:], np.array(h["predicted_weekly"])[n - 1:], color=colour, lw=1.8,
               label=f"Data to {pd.Timestamp(h['last_date_used']):%d %b}: {h['predicted_final_total']:.0f} total")
ax[0].axvline(resp_date, color=RED, ls="--", lw=1)
ax[0].set_title(f"Hold-out test: predictions made with early data only\n(actual total {obs_onsets.sum():.0f})", fontsize=10)
ax[0].set_ylabel("Cases per week"); ax[0].legend(fontsize=7.5); ax[0].set_ylim(0, 120)
tail_dates = week_dates[START_WEEK:]
ax[1].fill_between(tail_dates, s_lo, s_hi, color=BLUE, alpha=0.2, label="5-95% of 1,000 random runs")
ax[1].plot(tail_dates, s_mid, color=BLUE, lw=2, label="Median run")
ax[1].plot(tail_dates, obs_tail, "o", color="black", ms=4, label="Reported cases")
ax[1].axvline(resp_date, color=RED, ls="--", lw=1)
ax[1].set_title("Stochastic check: do the random runs cover what happened?", fontsize=10)
ax[1].set_ylabel("Cases per week"); ax[1].legend(fontsize=8)
for a in ax:
    a.tick_params(axis="x", rotation=30)
plt.tight_layout(); plt.savefig(ROOT / "results" / "fig_kikwit_tests.png", dpi=150); plt.close()

fig, ax = plt.subplots(1, 3, figsize=(13.5, 4.2))
for a, (country, w) in zip(ax, west.items()):
    dates = pd.to_datetime(w["dates"])
    a.plot(dates, w["reported"], "o", color="black", ms=3, label="Reported (WHO)")
    a.plot(dates, w["model"], color=BLUE, lw=2, label=f"Model, R0 = {w['r0']:.2f}")
    a.axvline(FIT_END, color=RED, ls="--", lw=1)
    a.set_yscale("log"); a.set_title(f"{country} 2014 (WHO team R0: {w['who_ert_r0']})", fontsize=10)
    a.tick_params(axis="x", rotation=30); a.legend(fontsize=8)
ax[0].set_ylabel("Cumulative cases (log scale)")
fig.suptitle("West Africa 2014: fitted to the left of the dashed line, projected to the right", fontsize=10)
plt.tight_layout(); plt.savefig(ROOT / "results" / "fig_west_africa.png", dpi=150); plt.close()
print("Saved results/results.json, calibrated/kikwit_1995.json and 3 figures")
