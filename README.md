# LEO Orbital Decay Simulator

A first-principles Python simulation of atmospheric-drag-induced orbital decay for satellites in Low Earth Orbit (LEO) — built from scratch with **no reliance on commercial simulation packages** (STK, GMAT, etc.).

Originally developed as part of my graduate coursework/research at the University of Florida and published as a reproducible artifact on Zenodo: **[DOI: 10.5281/zenodo.19704660](https://doi.org/10.5281/zenodo.19704660)**.

## What it does

- Models atmospheric density with an **exponential (barometric) atmosphere** approximation, piecewise over 300–600 km.
- Propagates altitude decay for a circular LEO orbit using the classical **King-Hele drag-decay relation**.
- Implements two numerical integrators from scratch: **explicit Euler** and **4th-order Runge-Kutta (RK4)**.
- Runs a **parametric sweep** across:
  - Ballistic coefficient `B = m / (Cd·A)`: **4.5 – 36.4 kg/m²**
  - Drag coefficient `Cd`: **1.5 – 3.0**
  - Initial altitude: **300 – 600 km**
- Includes a **case-study mode** reproducing the qualitative Feb 2022 Starlink loss event, where a geomagnetic storm inflates thermospheric density and sharply accelerates drag decay on a low insertion orbit.

## Usage

```bash
pip install -r requirements.txt   # no third-party deps beyond the stdlib

# Single-satellite propagation
python orbit_decay.py --altitude 400 --mass 400 --area 10 --cd 2.2 --days 30 --integrator rk4

# Parametric sweep -> CSV
python orbit_decay.py --sweep-csv sweep_results.csv
```

Example output:

```
Ballistic coefficient B = 18.18 kg/m^2
Start altitude: 400.0 km
After 30.00 days -> altitude: 374.409 km (Δh = -25.591 km)
```

`sweep_results.csv` in this repo contains a full 4×3×3 sweep (altitude × Cd × ballistic coefficient), with re-entry lifetime estimates ranging from **~8 days** (300 km, low B) to **multi-year** timescales (600 km, high B) — consistent with the "weeks to decades" spread reported in the original study.

## Background

This model contextualizes results against space-debris mitigation policy: the **IADC 25-year post-mission disposal guideline** and the **FCC 5-year disposal rule** for US-licensed LEO satellites, and was used to quantitatively analyze the February 2022 Starlink deployment, where 38 of 49 newly launched satellites were lost to unexpectedly high drag during a geomagnetic storm.

## Files

- `orbit_decay.py` — simulator, integrators, sweep, and Starlink case-study function
- `sweep_results.csv` — example sweep output
- `requirements.txt` — dependencies (stdlib only)

## Author

Arvind Kanagasabapathi Chandirakala — M.S. Aerospace & Mechanical Engineering, University of Florida
[LinkedIn](https://www.linkedin.com/in/arvind-kc-/)
