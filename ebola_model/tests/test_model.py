"""
Automated checks of the model and the calibration code.

Run from the project folder:   python3 -m pytest tests -q
"""
import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from episim_ebola import (EbolaParams, Response, simulate_deterministic,   # noqa: E402
                          simulate_stochastic, summarise, reproduction_numbers)
from episim_ebola.calibrate import fit_outbreak, predict_weekly, weekly_sums  # noqa: E402

P = EbolaParams()


# ---- the model is built correctly -------------------------------------------
def test_route_contributions_add_up_to_r0():
    b_i, b_h, b_f = P.betas
    total = (b_i * P.time_in_i + b_h * P.share_hospitalised * P.time_in_h
             + b_f * P.fatality_overall * P.burial_delay)
    assert total == pytest.approx(P.r0, rel=1e-12)
    assert reproduction_numbers(P)["total"] == pytest.approx(P.r0)


def test_deterministic_conserves_people():
    d = simulate_deterministic(P, 600)
    total = sum(d[k] for k in "SEIHFR")
    assert np.allclose(total, P.population)
    assert min(d[k].min() for k in "SEIHFR") >= -1e-9


def test_stochastic_conserves_people_and_stays_whole():
    s = simulate_stochastic(P, 300, runs=50, seed=1)
    total = sum(s[k] for k in "SEIHFR")
    assert (total == P.population).all()
    assert min(s[k].min() for k in "SEIHFR") >= 0


def test_final_size_matches_theory():
    """For any model of this kind, the share ever infected z solves z = 1 - exp(-R0 z)."""
    z = 0.5
    for _ in range(500):
        z = 1 - math.exp(-P.r0 * z)
    d = simulate_deterministic(P, 1500)
    assert (d["onsets"].sum() + 1) / P.population == pytest.approx(z, rel=0.005)


def test_deaths_match_the_fatality_set():
    d = simulate_deterministic(P, 1500)
    assert d["deaths"].sum() / (d["onsets"].sum() + 1) == pytest.approx(P.fatality_overall, rel=0.002)


def test_outbreak_fades_when_r0_is_below_one():
    d = simulate_deterministic(P.with_changes(r0=0.8), 600, initial_infectious=10)
    assert d["onsets"].sum() < 60
    assert d["I"][-1] < 0.01


# ---- the two solvers agree ---------------------------------------------------
def test_stochastic_average_matches_deterministic():
    d = simulate_deterministic(P, 150, initial_infectious=200)
    s = simulate_stochastic(P, 150, runs=600, seed=3, initial_infectious=200)
    assert s["onsets"].sum(axis=0).mean() == pytest.approx(d["onsets"].sum(), rel=0.04)
    assert s["deaths"].sum(axis=0).mean() == pytest.approx(d["deaths"].sum(), rel=0.04)


def test_same_seed_gives_same_result():
    a = simulate_stochastic(P, 120, runs=20, seed=42)
    b = simulate_stochastic(P, 120, runs=20, seed=42)
    c = simulate_stochastic(P, 120, runs=20, seed=43)
    assert (a["onsets"] == b["onsets"]).all()
    assert not (a["onsets"] == c["onsets"]).all()


def test_some_introductions_die_out_by_chance():
    s = summarise(simulate_stochastic(P, 400, runs=400, seed=5))
    assert 0.2 < s["share_died_out"] < 0.8


# ---- responses behave sensibly -----------------------------------------------
def test_response_reduces_cases_and_earlier_is_better():
    def total(day):
        r = (Response(start_day=day, community=0.4, funeral=0.2),)
        return simulate_deterministic(P, 700, r, initial_infectious=5)["onsets"].sum()
    none = simulate_deterministic(P, 700, initial_infectious=5)["onsets"].sum()
    assert total(30) < total(90) < none


def test_response_with_no_effect_changes_nothing():
    a = simulate_deterministic(P, 200)
    b = simulate_deterministic(P, 200, (Response(start_day=10),))
    assert np.allclose(a["onsets"], b["onsets"])


def test_response_window_and_ramp():
    r = Response(start_day=10, end_day=20, community=0.5, ramp_days=4)
    assert r.strength(9.9) == 0 and r.strength(12) == 0.5 and r.strength(15) == 1 and r.strength(20) == 0
    assert reproduction_numbers(P, (r,), day=15)["community"] == pytest.approx(P.r0 * P.share_community * 0.5)


def test_bad_inputs_are_rejected():
    with pytest.raises(ValueError):
        EbolaParams(share_community=0.5, share_hospital=0.5, share_funeral=0.5)
    with pytest.raises(ValueError):
        EbolaParams(fatality_home=1.2)
    with pytest.raises(ValueError):
        Response(start_day=0, community=1.5)


# ---- the calibration code can recover a known answer -------------------------
def test_calibration_recovers_known_parameters():
    """Make fake data from known parameters, add noise, and check the fit finds them."""
    base = EbolaParams(population=200_000)
    true = dict(r0=2.2, remaining=0.25, i0=1.5)
    mu = predict_weekly(base, true["r0"], true["remaining"], true["i0"], 100, 40, 20)
    rng = np.random.default_rng(11)
    fake = rng.poisson(mu)
    f = fit_outbreak(base, fake, 40, 100, starts=3)
    assert f["r0"] == pytest.approx(true["r0"], rel=0.15)
    assert f["r_after"] == pytest.approx(true["r0"] * true["remaining"], abs=0.25)


def test_weekly_sums():
    assert weekly_sums(np.arange(21), 7, 2).tolist() == [sum(range(7, 14)), sum(range(14, 21))]


# ---- the saved calibration still matches the data ----------------------------
def test_calibrated_preset_reproduces_kikwit():
    preset = json.load(open(ROOT / "calibrated" / "kikwit_1995.json"))
    p = EbolaParams(**preset["parameters"])
    r = preset["response_as_observed"]
    resp = (Response(start_day=r["start_day"], community=r["community"],
                     hospital=r["hospital"], funeral=r["funeral"], ramp_days=r["ramp_days"]),)
    d = simulate_deterministic(p, 192, resp, initial_infectious=preset["initial_infectious"])
    reported = p.reporting_fraction * weekly_sums(d["onsets"], 52, 20).sum()
    assert reported == pytest.approx(291, rel=0.05)         # cases with a known onset date
    reported_deaths = p.reporting_fraction * weekly_sums(d["deaths"], 52, 20).sum()
    assert reported_deaths == pytest.approx(236, rel=0.10)  # deaths were never fitted
    lo, hi = preset["calibration"]["r0_ci95"]
    assert lo < 2.7 < hi                                    # Legrand et al. 2007 estimate
