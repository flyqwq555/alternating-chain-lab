"""Full h in (0, 1] phase diagram and critical point solver for N=4."""

from __future__ import annotations

from dataclasses import dataclass
from typing import List

import numpy as np
from scipy.optimize import brentq, minimize

from .energy import (
    energy_per_particle,
    exact_n4_small_h_energy,
    periodic_gaussian_kernel_fourier,
)


@dataclass
class PhasePoint:
    h: float
    best_energy: float
    symmetric_energy: float
    gaps: List[float]
    regime: str


def compute_n4_critical_points() -> tuple[float, float]:
    """Find exact critical thresholds: h_cross (first-order transition) and h_instab (Hessian instability)."""
    def diff_sym_vs_3gap(h: float) -> float:
        return exact_n4_small_h_energy(h) - energy_per_particle([0.0, h, 2.0 * h, 3.0 * h])

    def d2e_ds2(h: float) -> float:
        k2_at_2 = periodic_gaussian_kernel_fourier(2.0, 4.0, order=2)
        k2_at_2_minus_h = periodic_gaussian_kernel_fourier(2.0 - h, 4.0, order=2)
        return -0.5 * (k2_at_2 - 2.0 * k2_at_2_minus_h)

    h_cross = float(brentq(diff_sym_vs_3gap, 0.60, 0.70, xtol=1e-14))
    h_instab = float(brentq(d2e_ds2, 0.60, 0.75, xtol=1e-14))
    return h_cross, h_instab


def classify_n4_gaps(gaps: List[float], h: float, tol: float = 1e-4) -> str:
    """Classify a 4-gap vector into its physical phase."""
    if abs(h - 1.0) < 1e-6:
        return "equidistant"
    s = sorted(gaps)
    locked = sum(1 for g in s if abs(g - h) < tol)
    if locked >= 3:
        return "three_gap_locked"
    if locked == 2 and abs(s[2] - s[3]) < tol:
        return "symmetric_1_3"
    return "symmetry_broken_1_2_1"


def optimize_n4_at_h(h: float, n_random: int = 24, seed: int = 42) -> tuple[float, List[float]]:
    """Globally minimize N=4 energy subject to minimum periodic spacing >= h."""
    if h >= 1.0 - 1e-12:
        pts = [0.0, 1.0, 2.0, 3.0]
        return energy_per_particle(pts), [1.0, 1.0, 1.0, 1.0]

    rng = np.random.default_rng(seed)
    candidates = [
        [h, h, 2.0 - h],
        [h, h, h],
        [h, 2.0 - h, 2.0 - h],
    ]
    rem = 4.0 - 4.0 * h
    for _ in range(n_random):
        r = rng.dirichlet([1.0, 1.0, 1.0, 1.0]) * rem
        candidates.append([h + r[0], h + r[1], h + r[2]])

    best_gaps = [h, h, 2.0 - h, 2.0 - h]
    best_e = energy_per_particle([0.0, 2.0 - h, 2.0, 2.0 + h])
    cons = ({"type": "ineq", "fun": lambda d: 4.0 - float(np.sum(d)) - h},)
    bnds = [(h, 4.0 - 3.0 * h)] * 3

    def obj(d: np.ndarray) -> float:
        return energy_per_particle([0.0, d[0], d[0] + d[1], d[0] + d[1] + d[2]])

    for d0 in candidates:
        res = minimize(
            obj,
            np.asarray(d0, dtype=float),
            method="SLSQP",
            bounds=bnds,
            constraints=cons,
            options={"ftol": 1e-13, "maxiter": 300},
        )
        d = [float(res.x[0]), float(res.x[1]), float(res.x[2])]
        gaps = d + [float(4.0 - sum(d))]
        # Reject infeasible SLSQP endpoints (all gaps must be >= h)
        if min(gaps) < h - 1e-9:
            continue
        e = energy_per_particle([0.0, d[0], d[0] + d[1], d[0] + d[1] + d[2]])
        if e < best_e:
            best_e = e
            best_gaps = gaps

    return best_e, best_gaps


def scan_n4_phase_diagram(
    h_min: float = 0.01, h_max: float = 1.0, num_points: int = 60
) -> tuple[list[PhasePoint], tuple[float, float]]:
    """Scan h in [h_min, h_max] and return phase points along with (h_c1, h_c2)."""
    h_c1, h_c2 = compute_n4_critical_points()
    hs = np.linspace(h_min, h_max, num_points)
    results: list[PhasePoint] = []
    for h in hs:
        h_val = float(h)
        best_e, gaps = optimize_n4_at_h(h_val)
        sym_e = exact_n4_small_h_energy(h_val)
        regime = classify_n4_gaps(gaps, h_val)
        results.append(
            PhasePoint(
                h=h_val,
                best_energy=best_e,
                symmetric_energy=sym_e,
                gaps=gaps,
                regime=regime,
            )
        )
    return results, (h_c1, h_c2)
