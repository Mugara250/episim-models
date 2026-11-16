"""
Calibration tools: fit the model to weekly case counts.

Method
    The deterministic solver gives the expected number of new cases per week.
    Observed weekly counts are assumed to scatter around that expectation with a
    negative binomial distribution (like Poisson, but allowing extra scatter).
    The best parameters are the ones that make the observed counts most likely
    (maximum likelihood). Ranges come from the profile likelihood.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import minimize, brentq
from scipy.special import gammaln
from scipy.stats import chi2, nbinom

from .model import EbolaParams, Response, simulate_deterministic

__all__ = ["weekly_sums", "negbin_loglik", "fit_outbreak", "profile_interval",
           "predict_weekly", "negbin_interval"]


def weekly_sums(daily: np.ndarray, first_day: int, weeks: int) -> np.ndarray:
    """Add daily values into `weeks` consecutive 7-day bins starting at first_day."""
    d = np.asarray(daily, dtype=float)[first_day:first_day + 7 * weeks]
    return d.reshape(weeks, 7).sum(axis=1)


def negbin_loglik(observed: np.ndarray, expected: np.ndarray, size: float) -> float:
    """Log-likelihood of counts with mean `expected` and variance mean + mean^2/size."""
    y = np.asarray(observed, dtype=float)
    mu = np.maximum(np.asarray(expected, dtype=float), 1e-9)
    return float(np.sum(gammaln(y + size) - gammaln(size) - gammaln(y + 1)
                        + size * np.log(size / (size + mu))
                        + y * np.log(mu / (size + mu))))


def negbin_interval(expected: np.ndarray, size: float, level: float = 0.90):
    """Range that should contain `level` of observations around each expected value."""
    mu = np.maximum(np.asarray(expected, dtype=float), 1e-9)
    prob = size / (size + mu)
    lo = nbinom.ppf((1 - level) / 2, size, prob)
    hi = nbinom.ppf(1 - (1 - level) / 2, size, prob)
    return lo, hi


def predict_weekly(base: EbolaParams, r0: float, remaining: float, i0: float,
                   response_day: float, first_day: int, weeks: int,
                   what: str = "onsets", ramp_days: float = 0.0) -> np.ndarray:
    """Expected reported weekly counts for a given parameter set."""
    p = base.with_changes(r0=r0)
    resp = (Response(start_day=response_day, community=remaining,
                     hospital=remaining, funeral=remaining, ramp_days=ramp_days),)
    sim = simulate_deterministic(p, first_day + 7 * weeks, resp, initial_infectious=i0)
    return base.reporting_fraction * weekly_sums(sim[what], first_day, weeks)


def _unpack(theta, fixed):
    """theta holds the free parameters on a transformed scale."""
    names = [n for n in ("r0", "remaining", "i0", "size") if n not in fixed]
    vals = dict(fixed)
    for n, t in zip(names, theta):
        vals[n] = 1 / (1 + np.exp(-t)) if n == "remaining" else float(np.exp(t))
    return vals


def fit_outbreak(base: EbolaParams, observed_weekly: np.ndarray, first_day: int,
                 response_day: float, fit_weeks: int | None = None,
                 fixed: dict | None = None, starts: int = 4, seed: int = 0,
                 ramp_days: float = 0.0) -> dict:
    """Fit r0, the share of transmission remaining after the response, the
    initial number infectious and the scatter (`size`).

    `fit_weeks` limits the fit to the first weeks (used for hold-out tests).
    `fixed` pins any of the four parameters, e.g. {"remaining": 1.0}.
    """
    fixed = dict(fixed or {})
    y_all = np.asarray(observed_weekly, dtype=float)
    weeks = len(y_all) if fit_weeks is None else fit_weeks
    y = y_all[:weeks]
    names = [n for n in ("r0", "remaining", "i0", "size") if n not in fixed]

    def nll(theta):
        v = _unpack(theta, fixed)
        if not (0.2 < v["r0"] < 15 and 1e-3 < v["i0"] < 500 and 0.05 < v["size"] < 1e4):
            return 1e12
        mu = predict_weekly(base, v["r0"], v["remaining"], v["i0"],
                            response_day, first_day, weeks, ramp_days=ramp_days)
        return -negbin_loglik(y, mu, v["size"])

    rng = np.random.default_rng(seed)
    guess = {"r0": np.log(2.0), "remaining": 0.0, "i0": np.log(1.0), "size": np.log(5.0)}
    best = None
    for s in range(starts):
        x0 = np.array([guess[n] for n in names])
        if s > 0:
            x0 = x0 + rng.normal(0, 0.6, size=len(names))
        res = minimize(nll, x0, method="Nelder-Mead",
                       options={"xatol": 1e-4, "fatol": 1e-6, "maxiter": 4000})
        if best is None or res.fun < best.fun:
            best = res
    v = _unpack(best.x, fixed)
    v = {k: float(x) for k, x in v.items()}
    v.update(loglik=-float(best.fun), weeks_fitted=int(weeks), free=names,
             ramp_days=float(ramp_days),
             r_after=v["r0"] * v["remaining"])
    return v


def profile_interval(base: EbolaParams, observed_weekly: np.ndarray, first_day: int,
                     response_day: float, fit: dict, name: str,
                     lo: float, hi: float, level: float = 0.95) -> tuple[float, float]:
    """Range of one parameter that the data cannot rule out (profile likelihood).

    For each trial value of the parameter, the other parameters are re-fitted.
    The range is where the likelihood stays within the chi-square cut-off.
    """
    cut = fit["loglik"] - chi2.ppf(level, 1) / 2

    def gap(value):
        f = fit_outbreak(base, observed_weekly, first_day, response_day,
                         fit_weeks=fit["weeks_fitted"], fixed={name: value}, starts=2,
                         ramp_days=fit.get("ramp_days", 0.0))
        return f["loglik"] - cut

    best = fit[name]
    low = brentq(gap, lo, best, xtol=1e-3) if gap(lo) < 0 else lo
    high = brentq(gap, best, hi, xtol=1e-3) if gap(hi) < 0 else hi
    return float(low), float(high)
