"""
EpiSim Ebola model: SEIHFR compartmental model with three routes of transmission.

Structure follows Legrand et al. (2007), "Understanding the dynamics of Ebola
epidemics", Epidemiology and Infection 135(4):610-621.

Compartments
    S  susceptible
    E  exposed (infected, not yet infectious)
    I  infectious in the community
    H  hospitalised
    F  dead but not yet buried (funeral)
    R  removed (recovered, or buried)

    S -> E -> I -> H -> F -> R
                \\-> F     \\-> R
                 \\-> R

Everyone who dies (at home or in hospital) passes through F, as in Legrand et al.
Safe burial is modelled as a reduction of funeral transmission, not as a bypass of F.

Two solvers share the same parameters:
    simulate_deterministic  one smooth average curve (used for calibration)
    simulate_stochastic     many random runs (used for scenario ranges)

This is a scenario model. It is not a forecasting tool. See MODEL.md for the
calibration, the tests, the limits and the sources of every value.
"""
from __future__ import annotations

from dataclasses import dataclass, field, asdict, replace
from typing import Iterable, Optional

import numpy as np

__all__ = ["EbolaParams", "Response", "simulate_deterministic",
           "simulate_stochastic", "summarise", "reproduction_numbers"]


@dataclass(frozen=True)
class EbolaParams:
    """All inputs of the model. Durations are in days, shares are 0 to 1."""
    population: int = 100_000
    r0: float = 1.8
    # Shares of R0 by route. They must add up to 1.
    share_community: float = 0.55
    share_hospital: float = 0.15
    share_funeral: float = 0.30
    incubation: float = 8.5          # infection -> symptoms
    onset_to_admission: float = 4.5  # symptoms -> hospital, for those admitted
    onset_to_death: float = 9.3
    onset_to_recovery: float = 13.0
    burial_delay: float = 2.0        # death -> burial
    share_hospitalised: float = 0.60
    fatality_home: float = 0.646     # cases never admitted
    fatality_hospital: float = 0.50  # admitted cases
    reporting_fraction: float = 1.0  # share of cases that appear in the data

    def __post_init__(self):
        total = self.share_community + self.share_hospital + self.share_funeral
        if abs(total - 1.0) > 1e-9:
            raise ValueError(f"route shares must add up to 1, got {total}")
        for name in ("share_hospitalised", "fatality_home", "fatality_hospital",
                     "reporting_fraction"):
            if not 0.0 <= getattr(self, name) <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")
        if self.onset_to_admission >= min(self.onset_to_death, self.onset_to_recovery):
            raise ValueError("onset_to_admission must be shorter than time to outcome")
        if self.share_hospital > 0 and self.share_hospitalised == 0:
            raise ValueError("hospital transmission needs some hospitalised cases")

    # ---- derived quantities -------------------------------------------------
    @property
    def time_in_i(self) -> float:
        """Average days a case is infectious in the community."""
        home = (self.fatality_home * self.onset_to_death
                + (1 - self.fatality_home) * self.onset_to_recovery)
        return (self.share_hospitalised * self.onset_to_admission
                + (1 - self.share_hospitalised) * home)

    @property
    def time_in_h(self) -> float:
        """Average days an admitted case spends in hospital."""
        return (self.fatality_hospital * self.onset_to_death
                + (1 - self.fatality_hospital) * self.onset_to_recovery
                - self.onset_to_admission)

    @property
    def exits_from_i(self) -> tuple[float, float, float]:
        """Shares leaving I to (hospital, funeral, recovery)."""
        h = self.share_hospitalised
        return h, (1 - h) * self.fatality_home, (1 - h) * (1 - self.fatality_home)

    @property
    def fatality_overall(self) -> float:
        h = self.share_hospitalised
        return h * self.fatality_hospital + (1 - h) * self.fatality_home

    @property
    def betas(self) -> tuple[float, float, float]:
        """Transmission rate per day for (community, hospital, funeral).
        rate = route's part of R0 / average days one case spends in that setting."""
        b_i = self.r0 * self.share_community / self.time_in_i
        b_h = (self.r0 * self.share_hospital / (self.share_hospitalised * self.time_in_h)
               if self.share_hospital > 0 else 0.0)
        b_f = (self.r0 * self.share_funeral / (self.fatality_overall * self.burial_delay)
               if self.share_funeral > 0 else 0.0)
        return b_i, b_h, b_f

    @property
    def generation_time(self) -> float:
        """Average days between a case's infection and the infections it causes."""
        t_i = self.time_in_i
        return (self.incubation
                + self.share_community * t_i
                + self.share_hospital * (t_i + self.time_in_h)
                + self.share_funeral * (t_i + self.burial_delay))

    def with_changes(self, **changes) -> "EbolaParams":
        return replace(self, **changes)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass(frozen=True)
class Response:
    """A response that changes transmission between two days.

    Each value is the share of transmission that REMAINS on that route
    (1.0 = no change, 0.2 = cut by 80%). Overlapping responses multiply.
    With `ramp_days` above zero, the effect builds up evenly over that many days
    instead of switching on at once."""
    start_day: float
    end_day: Optional[float] = None
    community: float = 1.0
    hospital: float = 1.0
    funeral: float = 1.0
    ramp_days: float = 0.0
    label: str = ""

    def __post_init__(self):
        for name in ("community", "hospital", "funeral"):
            if not 0.0 <= getattr(self, name) <= 1.0:
                raise ValueError(f"{name} must be between 0 and 1")

    def active(self, day: float) -> bool:
        return day >= self.start_day and (self.end_day is None or day < self.end_day)

    def strength(self, day: float) -> float:
        """0 before the start, rising to 1 once the ramp is complete."""
        if not self.active(day):
            return 0.0
        if self.ramp_days <= 0:
            return 1.0
        return min(1.0, (day - self.start_day) / self.ramp_days)


def _multipliers(responses: Iterable[Response], day: float) -> tuple[float, float, float]:
    m_i = m_h = m_f = 1.0
    for r in responses:
        w = r.strength(day)
        if w > 0:
            m_i *= 1 - w * (1 - r.community)
            m_h *= 1 - w * (1 - r.hospital)
            m_f *= 1 - w * (1 - r.funeral)
    return m_i, m_h, m_f


def reproduction_numbers(params: EbolaParams,
                         responses: Iterable[Response] = (), day: float = 0.0) -> dict:
    """R by route at a given day, in a fully susceptible population."""
    m_i, m_h, m_f = _multipliers(responses, day)
    c = params.r0 * params.share_community * m_i
    h = params.r0 * params.share_hospital * m_h
    f = params.r0 * params.share_funeral * m_f
    return {"community": c, "hospital": h, "funeral": f, "total": c + h + f}


def simulate_deterministic(params: EbolaParams, days: int,
                           responses: Iterable[Response] = (),
                           initial_infectious: float = 1.0,
                           initial_exposed: float = 0.0,
                           steps_per_day: int = 10) -> dict:
    """Average behaviour of the model. Returns daily arrays of length days + 1.

    'onsets' and 'deaths' are the new events during each day (index 0 is zero).
    """
    responses = tuple(responses)
    p = params
    n = float(p.population)
    b_i, b_h, b_f = p.betas
    to_h, to_f, to_r = p.exits_from_i
    dt = 1.0 / steps_per_day
    # Rates per step, worked out once
    r_ei, r_i = dt / p.incubation, dt / p.time_in_i
    r_h, r_f = dt / p.time_in_h, dt / p.burial_delay
    f_h = p.fatality_hospital
    # Response multipliers only change at a few moments, so look them up per step
    # only when a response exists.
    if responses:
        mult = [_multipliers(responses, k * dt) for k in range(days * steps_per_day)]
    else:
        mult = None

    S = n - initial_infectious - initial_exposed
    E, I, H, F, R = float(initial_exposed), float(initial_infectious), 0.0, 0.0, 0.0
    names = ("S", "E", "I", "H", "F", "R", "onsets", "deaths", "admissions")
    rows = [[S, E, I, H, F, R, 0.0, 0.0, 0.0]]
    k = 0
    for day in range(1, days + 1):
        onsets = deaths = admissions = 0.0
        for _ in range(steps_per_day):
            if mult is None:
                pressure = b_i * I + b_h * H + b_f * F
            else:
                m_i, m_h, m_f = mult[k]
                pressure = b_i * m_i * I + b_h * m_h * H + b_f * m_f * F
            k += 1
            new_e = pressure * S / n * dt
            new_i = E * r_ei
            leave_i = I * r_i
            leave_h = H * r_h
            buried = F * r_f
            die_h = leave_h * f_h
            die_home = leave_i * to_f
            adm = leave_i * to_h

            S -= new_e
            E += new_e - new_i
            I += new_i - leave_i
            H += adm - leave_h
            F += die_home + die_h - buried
            R += leave_i * to_r + (leave_h - die_h) + buried
            onsets += new_i
            deaths += die_home + die_h
            admissions += adm
        rows.append([S, E, I, H, F, R, onsets, deaths, admissions])
    arr = np.array(rows)
    return {name: arr[:, j] for j, name in enumerate(names)}


def simulate_stochastic(params: EbolaParams, days: int,
                        responses: Iterable[Response] = (),
                        runs: int = 200, seed: int = 2026,
                        initial_infectious: int = 1, initial_exposed: int = 0,
                        initial_state: Optional[dict] = None,
                        steps_per_day: int = 4) -> dict:
    """`runs` independent random outbreaks. Returns arrays of shape (days + 1, runs).

    Time moves in short steps (4 per day by default). In each step, each person
    has a chance of moving, drawn at random. The chance of leaving a compartment
    in one step is step / average stay, which makes the average stay exactly
    right, so the model's R0 equals the R0 that was set.

    `initial_state` (keys S, E, I, H, F, R) overrides the initial counts, which
    lets a simulation start from the middle of an outbreak.
    """
    responses = tuple(responses)
    p = params
    rng = np.random.default_rng(seed)
    n = p.population
    b_i, b_h, b_f = p.betas
    exits = p.exits_from_i
    dt = 1.0 / steps_per_day
    if dt > min(p.incubation, p.time_in_i, p.time_in_h, p.burial_delay):
        raise ValueError("steps_per_day is too small for the shortest duration")

    if initial_state is None:
        initial_state = {"S": n - initial_infectious - initial_exposed,
                         "E": initial_exposed, "I": initial_infectious,
                         "H": 0, "F": 0, "R": 0}
    start = {k: int(round(initial_state[k])) for k in ("S", "E", "I", "H", "F", "R")}
    total = sum(start.values())
    S, E, I, H, F, R = (np.full(runs, start[k], dtype=np.int64)
                        for k in ("S", "E", "I", "H", "F", "R"))

    p_ei = dt / p.incubation
    p_i = dt / p.time_in_i
    p_h = dt / p.time_in_h
    p_f = dt / p.burial_delay
    home_death = exits[1] / (exits[1] + exits[2]) if exits[1] + exits[2] > 0 else 0.0

    out = {k: np.zeros((days + 1, runs), dtype=np.int64) for k in
           ("S", "E", "I", "H", "F", "R", "onsets", "deaths", "admissions")}
    for k, v in (("S", S), ("E", E), ("I", I), ("H", H), ("F", F), ("R", R)):
        out[k][0] = v

    for day in range(1, days + 1):
        onsets = np.zeros(runs, dtype=np.int64)
        deaths = np.zeros(runs, dtype=np.int64)
        admissions = np.zeros(runs, dtype=np.int64)
        for step in range(steps_per_day):
            m_i, m_h, m_f = _multipliers(responses, day - 1 + step * dt)
            force = (b_i * m_i * I + b_h * m_h * H + b_f * m_f * F) / total
            new_e = rng.binomial(S, 1.0 - np.exp(-force * dt))
            new_i = rng.binomial(E, p_ei)
            leave_i = rng.binomial(I, p_i)
            leave_h = rng.binomial(H, p_h)
            buried = rng.binomial(F, p_f)
            to_h = rng.binomial(leave_i, exits[0])
            to_f = rng.binomial(leave_i - to_h, home_death)
            to_r = leave_i - to_h - to_f
            die_h = rng.binomial(leave_h, p.fatality_hospital)

            S = S - new_e
            E = E + new_e - new_i
            I = I + new_i - leave_i
            H = H + to_h - leave_h
            F = F + to_f + die_h - buried
            R = R + to_r + (leave_h - die_h) + buried
            onsets += new_i
            deaths += to_f + die_h
            admissions += to_h
        if not np.all(S + E + I + H + F + R == total):
            raise RuntimeError("people were created or lost")
        for k, v in (("S", S), ("E", E), ("I", I), ("H", H), ("F", F), ("R", R),
                     ("onsets", onsets), ("deaths", deaths), ("admissions", admissions)):
            out[k][day] = v
    return out


def summarise(result: dict, takeoff_threshold: int = 20) -> dict:
    """Turn stochastic output into the numbers a dashboard shows.

    Returns daily 5th, 50th and 95th percentiles, plus totals. Runs that end
    with fewer than `takeoff_threshold` cases are counted as 'died out'.
    """
    onsets = result["onsets"]
    cum_cases = onsets.cumsum(axis=0) + (result["I"][0] + result["H"][0])
    cum_deaths = result["deaths"].cumsum(axis=0)
    final = cum_cases[-1]
    took_off = final >= takeoff_threshold

    def bands(a):
        lo, mid, hi = np.percentile(a, [5, 50, 95], axis=1)
        return {"p05": lo.tolist(), "median": mid.tolist(), "p95": hi.tolist()}

    def stats(v):
        if len(v) == 0:
            return None
        lo, mid, hi = np.percentile(v, [5, 50, 95])
        return {"p05": float(lo), "median": float(mid), "p95": float(hi)}

    return {
        "runs": int(onsets.shape[1]),
        "days": int(onsets.shape[0] - 1),
        "share_died_out": float(1 - took_off.mean()),
        "daily": {
            "new_cases": bands(onsets),
            "cumulative_cases": bands(cum_cases),
            "cumulative_deaths": bands(cum_deaths),
            "in_hospital": bands(result["H"]),
        },
        "totals_all_runs": {
            "cases": stats(final), "deaths": stats(cum_deaths[-1]),
            "peak_in_hospital": stats(result["H"].max(axis=0)),
        },
        "totals_runs_that_took_off": {
            "cases": stats(final[took_off]), "deaths": stats(cum_deaths[-1][took_off]),
            "peak_in_hospital": stats(result["H"].max(axis=0)[took_off]),
        },
    }
