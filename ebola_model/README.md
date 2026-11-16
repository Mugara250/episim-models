# EpiSim Ebola model

A compartmental Ebola model (SEIHFR, three routes of transmission), calibrated to the Kikwit 1995 outbreak and tested against real data. Built for the EpiSim demo.

**Read `MODEL.md` first.** It describes the model, the calibration, every test result (including the ones that failed), the limits and the sources.

Status: calibrated to one outbreak; not a forecasting tool; not reviewed by an epidemiologist.

## Contents

| Path | What it is |
|---|---|
| `episim_ebola/model.py` | The model: parameters, responses, deterministic and stochastic solvers, summary for dashboards |
| `episim_ebola/calibrate.py` | Fitting tools |
| `calibrated/kikwit_1995.json` | The calibrated parameter set |
| `data/` | The two outbreak datasets and `SOURCES.md` |
| `data_raw/` | The files as downloaded |
| `scripts/run_calibration.py` | Reproduces the calibration and all tests (about 5 minutes) |
| `scripts/west_africa_window_check.py` | Sensitivity of the West Africa estimates to the fitting window |
| `scripts/demo_scenarios.py` | Example scenario comparison for the demo |
| `tests/test_model.py` | 16 automated tests |
| `results/` | Numbers (`results.json`) and figures |

## Set-up

```
pip install numpy scipy pandas matplotlib pytest
python3 -m pytest tests -q
python3 scripts/demo_scenarios.py
```

The model itself needs only NumPy. SciPy and pandas are used for calibration, matplotlib for figures.

## Using the model

```python
import json
from episim_ebola import EbolaParams, Response, simulate_stochastic, summarise

preset = json.load(open("calibrated/kikwit_1995.json"))
params = EbolaParams(**preset["parameters"]).with_changes(
    population=350_000,          # the district being studied
    reporting_fraction=1.0)

responses = [
    Response(start_day=45, funeral=0.2, label="Safe burials"),            # cut by 80%
    Response(start_day=60, end_day=180, community=0.5, ramp_days=14,
             label="Isolation and contact tracing"),
]

result = simulate_stochastic(params, days=365, responses=responses,
                             runs=200, seed=2026, initial_infectious=5)
summary = summarise(result)      # plain dict, ready to return as JSON
```

`summary` contains, for each day, the 5th, 50th and 95th percentile of new cases, cumulative cases, cumulative deaths and people in hospital, plus totals and the share of runs that died out.

A response value is the share of transmission that **remains** on a route: 1.0 means no change, 0.2 means cut by 80%.

## Connecting it to the EpiSim backend (FastAPI)

The backend's simulation worker can call the model directly:

```python
def run_scenario(scenario: dict) -> dict:
    params = EbolaParams(**scenario["disease_parameters"]).with_changes(
        population=scenario["population"])
    responses = [Response(**r) for r in scenario["responses"]]
    result = simulate_stochastic(params, days=scenario["days"], responses=responses,
                                 runs=scenario["runs"], seed=scenario["seed"],
                                 initial_infectious=scenario["first_cases"])
    return summarise(result)
```

- Store the scenario's inputs and the seed with every run, so a result can be reproduced.
- 200 runs of 365 days take under a second.
- Hospital capacity is applied afterwards, by comparing `in_hospital` with the beds available.
- The intervention effect sizes in the EpiSim library are placeholders. This model does not estimate the effect of individual measures.

## Wording for the demo

"Model calibrated to the Kikwit 1995 outbreak. Scenario estimates, not forecasts. Not yet reviewed by an epidemiologist."

## Data credits

Kikwit 1995: Khan et al. (1999), via the `outbreaks` R package. West Africa 2014: WHO Disease Outbreak News and situation reports, compiled by C. Rivers. See `data/SOURCES.md`.
