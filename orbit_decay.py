"""
LEO Orbital Decay Simulator
============================
A first-principles Python model of atmospheric-drag-induced orbital decay
for satellites in Low Earth Orbit (LEO), built from scratch with no
reliance on commercial simulation packages (e.g. STK, GMAT).

Physics
-------
Drag deceleration on a circular orbit:
    a_drag = 0.5 * rho(h) * v^2 * (Cd * A / m)
           = 0.5 * rho(h) * v^2 / B          where B = m / (Cd * A) is the
                                              ballistic coefficient [kg/m^2]

Atmospheric density is modelled with a simple exponential (barometric)
approximation:
    rho(h) = rho0 * exp(-(h - h0) / H)

Two numerical integrators are provided for the altitude/semi-major-axis
decay ODE: explicit Euler and classical 4th-order Runge-Kutta (RK4).

Author: Arvind Kanagasabapathi Chandirakala
"""

from __future__ import annotations

import argparse
import csv
import math
from dataclasses import dataclass

MU_EARTH = 3.986004418e14  # m^3/s^2, Earth's gravitational parameter
R_EARTH = 6378137.0        # m, mean equatorial radius

# Reference exponential-atmosphere table (base altitude km : (rho0 kg/m^3, scale height km))
# Loosely follows the piecewise exponential model used in NRLMSISE-00-style summaries.
ATMOSPHERE_TABLE = [
    (300, 1.916e-11, 53.628),
    (350, 7.014e-12, 53.298),
    (400, 2.803e-12, 58.515),
    (450, 1.184e-12, 60.828),
    (500, 5.215e-13, 63.822),
    (550, 2.384e-13, 71.835),
    (600, 1.137e-13, 88.667),
]


def atmospheric_density(h_km: float) -> float:
    """Exponential-model atmospheric density [kg/m^3] at altitude h_km [km]."""
    h_km = max(h_km, ATMOSPHERE_TABLE[0][0])
    band = ATMOSPHERE_TABLE[0]
    for entry in ATMOSPHERE_TABLE:
        if h_km >= entry[0]:
            band = entry
        else:
            break
    h0, rho0, scale_h = band
    return rho0 * math.exp(-(h_km - h0) / scale_h)


def orbital_velocity(h_km: float) -> float:
    """Circular orbital velocity [m/s] at altitude h_km."""
    r = R_EARTH + h_km * 1000.0
    return math.sqrt(MU_EARTH / r)


@dataclass
class SatelliteConfig:
    mass_kg: float
    area_m2: float
    cd: float
    storm_multiplier: float = 1.0  # density scaling for a geomagnetic-storm case

    @property
    def ballistic_coefficient(self) -> float:
        """B = m / (Cd * A)  [kg/m^2]"""
        return self.mass_kg / (self.cd * self.area_m2)


def altitude_decay_rate(h_km: float, sat: SatelliteConfig) -> float:
    """dh/dt [km/s] due to atmospheric drag at circular altitude h_km."""
    rho = atmospheric_density(h_km) * sat.storm_multiplier
    v = orbital_velocity(h_km)
    b = sat.ballistic_coefficient
    # Energy-balance approximation for circular-orbit decay:
    # dh/dt = -rho * v * (mu / (2*B)) / v  simplifies to the standard King-Hele form
    r = R_EARTH + h_km * 1000.0
    # Simplified King-Hele decay-rate relation for a circular orbit:
    dh_dt_m_s = -(rho * v) / b * r
    return dh_dt_m_s / 1000.0  # convert m/s -> km/s


def propagate_euler(h0_km: float, sat: SatelliteConfig, dt_s: float, t_end_s: float):
    t, h = 0.0, h0_km
    history = [(t, h)]
    while t < t_end_s and h > ATMOSPHERE_TABLE[0][0] * 0.5:
        dh = altitude_decay_rate(h, sat) * dt_s
        h += dh
        t += dt_s
        history.append((t, h))
    return history


def propagate_rk4(h0_km: float, sat: SatelliteConfig, dt_s: float, t_end_s: float):
    t, h = 0.0, h0_km
    history = [(t, h)]
    f = altitude_decay_rate
    while t < t_end_s and h > ATMOSPHERE_TABLE[0][0] * 0.5:
        k1 = f(h, sat)
        k2 = f(h + 0.5 * dt_s * k1, sat)
        k3 = f(h + 0.5 * dt_s * k2, sat)
        k4 = f(h + dt_s * k3, sat)
        h += (dt_s / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt_s
        history.append((t, h))
    return history


def estimate_lifetime_days(h0_km: float, sat: SatelliteConfig, dt_s: float = 3600.0,
                            max_days: float = 365 * 25) -> float:
    """Estimate re-entry lifetime in days using RK4, capped at max_days.

    Uses a coarse (1-hour) step for the long-horizon lifetime sweep -- fine
    enough given the exponential-atmosphere model's own uncertainty, and fast
    enough to sweep many (B, Cd, altitude) cases.
    """
    t, h = 0.0, h0_km
    t_end_s = max_days * 86400.0
    floor_km = ATMOSPHERE_TABLE[0][0] * 0.5
    f = altitude_decay_rate
    while t < t_end_s and h > floor_km:
        k1 = f(h, sat)
        k2 = f(h + 0.5 * dt_s * k1, sat)
        k3 = f(h + 0.5 * dt_s * k2, sat)
        k4 = f(h + dt_s * k3, sat)
        h += (dt_s / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4)
        t += dt_s
    return t / 86400.0


def starlink_feb2022_case_study():
    """
    Reproduces the qualitative Feb 2022 Starlink loss scenario: a geomagnetic
    storm inflates thermospheric density, sharply increasing drag on newly
    deployed satellites in a low (~210 km) insertion orbit.
    """
    quiet = SatelliteConfig(mass_kg=260.0, area_m2=8.0, cd=2.2, storm_multiplier=1.0)
    storm = SatelliteConfig(mass_kg=260.0, area_m2=8.0, cd=2.2, storm_multiplier=4.0)
    h0 = 210.0
    quiet_hist = propagate_rk4(h0, quiet, dt_s=30.0, t_end_s=10 * 86400.0)
    storm_hist = propagate_rk4(h0, storm, dt_s=30.0, t_end_s=10 * 86400.0)
    return quiet_hist, storm_hist


def sweep_and_export(csv_path: str):
    """Parametric sweep across ballistic coefficient, Cd, and altitude, per resume scope:
    B in [4.5, 36.4] kg/m^2, Cd in [1.5, 3.0], altitude in [300, 600] km."""
    rows = []
    for h0 in (300, 400, 500, 600):
        for cd in (1.5, 2.2, 3.0):
            for mass_area_ratio in (4.5, 15.0, 36.4):  # stand-in for m/A giving target B
                area = 1.0
                mass = mass_area_ratio * cd * area  # yields B = mass_area_ratio
                sat = SatelliteConfig(mass_kg=mass, area_m2=area, cd=cd)
                lifetime = estimate_lifetime_days(h0, sat)
                rows.append({
                    "altitude_km": h0,
                    "cd": cd,
                    "ballistic_coefficient_kg_m2": round(sat.ballistic_coefficient, 2),
                    "lifetime_days": round(lifetime, 1),
                })
    with open(csv_path, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    return rows


def main():
    parser = argparse.ArgumentParser(description="LEO orbital decay simulator")
    parser.add_argument("--altitude", type=float, default=400.0, help="initial altitude [km]")
    parser.add_argument("--mass", type=float, default=400.0, help="satellite mass [kg]")
    parser.add_argument("--area", type=float, default=10.0, help="cross-sectional area [m^2]")
    parser.add_argument("--cd", type=float, default=2.2, help="drag coefficient")
    parser.add_argument("--days", type=float, default=30.0, help="simulation duration [days]")
    parser.add_argument("--integrator", choices=["euler", "rk4"], default="rk4")
    parser.add_argument("--sweep-csv", type=str, default=None,
                         help="run the parametric sweep and export results to this CSV path")
    args = parser.parse_args()

    if args.sweep_csv:
        rows = sweep_and_export(args.sweep_csv)
        print(f"Wrote {len(rows)} sweep rows to {args.sweep_csv}")
        return

    sat = SatelliteConfig(mass_kg=args.mass, area_m2=args.area, cd=args.cd)
    propagate = propagate_rk4 if args.integrator == "rk4" else propagate_euler
    hist = propagate(args.altitude, sat, dt_s=60.0, t_end_s=args.days * 86400.0)
    t_final, h_final = hist[-1]
    print(f"Ballistic coefficient B = {sat.ballistic_coefficient:.2f} kg/m^2")
    print(f"Start altitude: {args.altitude:.1f} km")
    print(f"After {t_final/86400:.2f} days -> altitude: {h_final:.3f} km "
          f"(Δh = {h_final - args.altitude:.3f} km)")


if __name__ == "__main__":
    main()
