"""Automated verification tests for alternating-chain-lab."""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from alternating_chain.energy import (
    energy_per_particle,
    exact_n4_small_h_energy,
    n4_taylor_coefficients,
    periodic_gaussian_kernel_direct,
    periodic_gaussian_kernel_fourier,
)
from alternating_chain.multi_particle import compare_cluster_ansatz, optimize_chain
from alternating_chain.phase_diagram import compute_n4_critical_points, optimize_n4_at_h
from alternating_chain.rigour import run_rigorous_checks


def test_kernel_direct_vs_fourier() -> None:
    for x in [0.0, 0.37, 1.5, 2.0]:
        k_dir = periodic_gaussian_kernel_direct(x, period=4.0, m_max=20)
        k_fou = periodic_gaussian_kernel_fourier(x, period=4.0, n_max=25)
        assert abs(k_dir - k_fou) < 1e-14


def test_n4_literature_counterexample() -> None:
    e_eq = energy_per_particle([0.0, 1.0, 2.0, 3.0])
    e_cl = energy_per_particle([0.0, 1.5, 2.0, 2.5])
    assert abs(e_eq - 0.999633292158) < 1e-11
    assert abs(e_cl - 0.924751872822) < 1e-11
    assert e_eq - e_cl > 0.07478


def test_n4_taylor_coefficients() -> None:
    c0, c2, c4 = n4_taylor_coefficients()
    assert abs(c0 - 0.849687099565508) < 1e-14
    assert abs(c2 - 0.370882223854927) < 1e-14
    assert abs(c4 - (-0.307271003509848)) < 1e-14

    h = 0.01
    e_exact = exact_n4_small_h_energy(h)
    e_series = c0 + c2 * (h ** 2) + c4 * (h ** 4)
    assert abs(e_exact - e_series) < 1e-11


def test_n4_phase_transitions() -> None:
    h_cross, h_instab = compute_n4_critical_points()
    assert 0.6680 < h_cross < 0.6681
    assert 0.6684 < h_instab < 0.6685

    e_low, _ = optimize_n4_at_h(0.65)
    assert abs(e_low - exact_n4_small_h_energy(0.65)) < 1e-12

    e_high, _ = optimize_n4_at_h(0.75)
    assert abs(e_high - energy_per_particle([0.0, 0.75, 1.5, 2.25])) < 1e-12


def test_multi_n_beats_equidistant() -> None:
    for n in [6, 8, 10]:
        e_best, _ = optimize_chain(n, h=0.1, n_starts=15)
        ansatz = compare_cluster_ansatz(n, h=0.1)
        assert e_best < ansatz["equidistant"] - 0.05


def test_rigorous_mpmath_checks() -> None:
    rep = run_rigorous_checks(dps=50)
    assert rep.all_passed
