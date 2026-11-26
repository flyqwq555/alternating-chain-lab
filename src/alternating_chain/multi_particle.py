"""Extension to N=4, 6, 8, 10 alternating particle chains."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Sequence

import numpy as np
from scipy.optimize import minimize

from .energy import energy_per_particle


@dataclass
class MultiNResult:
    n: int
    h: float
    best_energy: float
    equidistant_energy: float
    best_gaps: List[float]
    ansatz_energies: Dict[str, float]


def _gaps_to_positions(gaps: Sequence[float]) -> list[float]:
    pts = [0.0]
    for g in gaps[:-1]:
        pts.append(pts[-1] + float(g))
    return pts


def compare_cluster_ansatz(n: int, h: float) -> Dict[str, float]:
    """Evaluate canonical cluster ansatz energies for even N and spacing h."""
    if n % 2 != 0 or n < 4:
        raise ValueError(f"N must be even >= 4, got {n}.")
    out: Dict[str, float] = {}
    out["equidistant"] = energy_per_particle(list(range(n)))

    m = n // 2
    if (n - 2) * h < n:
        g_1 = [n / 2.0 - (n - 2) * h / 2.0] + [h] * (n - 2) + [n / 2.0 - (n - 2) * h / 2.0]
        if min(g_1) >= h - 1e-12:
            out[f"1+{n-1}"] = energy_per_particle(_gaps_to_positions(g_1))

    g_eq = [h] * (m - 1) + [n / 2.0 - (m - 1) * h] + [h] * (m - 1) + [n / 2.0 - (m - 1) * h]
    if min(g_eq) >= h - 1e-12:
        out[f"{m}+{m}"] = energy_per_particle(_gaps_to_positions(g_eq))

    if m % 2 == 0 and m >= 2:
        g_odd = [h] * (m - 2) + [n / 2.0 - (m - 2) * h] + [h] * m + [n / 2.0 - m * h]
        if min(g_odd) >= h - 1e-12:
            out[f"{m-1}+{m+1}"] = energy_per_particle(_gaps_to_positions(g_odd))

    return out


def optimize_chain(n: int, h: float, n_starts: int = 30, seed: int = 42) -> tuple[float, List[float]]:
    """Multi-start SLSQP optimization for N alternating particles with min gap h."""
    if h >= 1.0 - 1e-12:
        return energy_per_particle(list(range(n))), [1.0] * n

    rng = np.random.default_rng(seed)
    inits: list[list[float]] = []
    rem = float(n) - n * h
    for _ in range(n_starts):
        r = rng.dirichlet(np.ones(n)) * rem
        inits.append(list((h + r)[:-1]))

    m = n // 2
    g_eq = [h] * (m - 1) + [n / 2.0 - (m - 1) * h] + [h] * (m - 1) + [n / 2.0 - (m - 1) * h]
    if min(g_eq) >= h - 1e-12:
        inits.append(g_eq[:-1])
    if m % 2 == 0 and m >= 2:
        g_odd = [h] * (m - 2) + [n / 2.0 - (m - 2) * h] + [h] * m + [n / 2.0 - m * h]
        if min(g_odd) >= h - 1e-12:
            inits.append(g_odd[:-1])

    best_e = float("inf")
    best_gaps = [1.0] * n
    cons = ({"type": "ineq", "fun": lambda d: float(n) - float(np.sum(d)) - h},)
    bnds = [(h, float(n) - (n - 1) * h)] * (n - 1)

    def obj(d: np.ndarray) -> float:
        pts = [0.0]
        for g in d:
            pts.append(pts[-1] + float(g))
        return energy_per_particle(pts)

    for d0 in inits:
        res = minimize(
            obj,
            np.asarray(d0, dtype=float),
            method="SLSQP",
            bounds=bnds,
            constraints=cons,
            options={"ftol": 1e-12, "maxiter": 300},
        )
        if res.fun < best_e:
            best_e = float(res.fun)
            best_gaps = [float(x) for x in res.x] + [float(n - np.sum(res.x))]

    return best_e, best_gaps


def scan_multi_n(
    n_values: Sequence[int] = (4, 6, 8, 10),
    h_values: Sequence[float] = (0.05, 0.2, 0.5, 0.75, 1.0),
) -> list[MultiNResult]:
    """Run global optimization and ansatz comparison across multiple N and h."""
    results: list[MultiNResult] = []
    for n in n_values:
        for h in h_values:
            best_e, gaps = optimize_chain(n, h)
            ansatz = compare_cluster_ansatz(n, h)
            results.append(
                MultiNResult(
                    n=n,
                    h=h,
                    best_energy=best_e,
                    equidistant_energy=ansatz["equidistant"],
                    best_gaps=gaps,
                    ansatz_energies=ansatz,
                )
            )
    return results
