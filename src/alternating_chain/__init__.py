"""Alternating Gaussian particle chains: phase diagrams, multi-N optimization, and rigorous verification."""

from .energy import (
    energy_per_particle,
    exact_n4_small_h_energy,
    n4_taylor_coefficients,
    periodic_gaussian_kernel_direct,
    periodic_gaussian_kernel_fourier,
)
from .phase_diagram import (
    compute_n4_critical_points,
    scan_n4_phase_diagram,
)
from .multi_particle import (
    compare_cluster_ansatz,
    optimize_chain,
    scan_multi_n,
)
from .rigour import run_rigorous_checks

__all__ = [
    "energy_per_particle",
    "exact_n4_small_h_energy",
    "n4_taylor_coefficients",
    "periodic_gaussian_kernel_direct",
    "periodic_gaussian_kernel_fourier",
    "compute_n4_critical_points",
    "scan_n4_phase_diagram",
    "compare_cluster_ansatz",
    "optimize_chain",
    "scan_multi_n",
    "run_rigorous_checks",
]
