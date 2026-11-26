"""Core periodic Gaussian kernel and energy definitions for alternating particle chains."""

from __future__ import annotations

import math
from typing import Sequence

import numpy as np


def periodic_gaussian_kernel_direct(x: float, period: float, m_max: int = 20) -> float:
    """Compute K_L(x) = sum_{m in Z} exp(-(x + m*L)^2 / 4) by direct lattice summation."""
    total = 0.0
    for m in range(-m_max, m_max + 1):
        d = x + m * period
        total += math.exp(-0.25 * d * d)
    return total


def periodic_gaussian_kernel_fourier(
    x: float, period: float, n_max: int = 25, order: int = 0
) -> float:
    """Compute K_L(x) or its derivative of given order via Poisson summation."""
    c = 2.0 * math.sqrt(math.pi) / period
    w0 = 2.0 * math.pi / period
    alpha = w0 * w0

    if order == 0:
        s = 1.0
        for n in range(1, n_max + 1):
            s += 2.0 * math.exp(-alpha * n * n) * math.cos(w0 * n * x)
        return c * s

    s = 0.0
    for n in range(1, n_max + 1):
        wn = w0 * n
        coeff = 2.0 * math.exp(-alpha * n * n) * (wn ** order)
        mod = order % 4
        if mod == 1:
            s -= coeff * math.sin(wn * x)
        elif mod == 2:
            s -= coeff * math.cos(wn * x)
        elif mod == 3:
            s += coeff * math.sin(wn * x)
        else:
            s += coeff * math.cos(wn * x)
    return c * s


def energy_per_particle(positions: Sequence[float], m_max: int = 20) -> float:
    """Compute normalized energy per particle for N alternating particles on [0, N)."""
    pts = np.asarray(positions, dtype=float)
    n = len(pts)
    if n % 2 != 0 or n < 2:
        raise ValueError(f"Particle count N must be an even integer >= 2, got {n}.")
    period = float(n)
    signs = np.where(np.arange(n) % 2 == 0, 1.0, -1.0)
    m_vals = np.arange(-m_max, m_max + 1, dtype=float)[:, None, None]
    diffs = pts[:, None] - pts[None, :]
    kernel_mat = np.sum(np.exp(-0.25 * (diffs[None, :, :] + m_vals * period) ** 2), axis=0)
    total = float(signs @ kernel_mat @ signs)
    return 1.0 - total / n


def exact_n4_small_h_energy(h: float) -> float:
    """Exact energy of the symmetric (0, 2-h, 2, 2+h) configuration for N=4, L=4."""
    k0 = periodic_gaussian_kernel_fourier(0.0, 4.0)
    k2 = periodic_gaussian_kernel_fourier(2.0, 4.0)
    kh = periodic_gaussian_kernel_fourier(h, 4.0)
    k2h = periodic_gaussian_kernel_fourier(2.0 * h, 4.0)
    k2_minus_h = periodic_gaussian_kernel_fourier(2.0 - h, 4.0)
    return 1.0 - k0 - 0.5 * (k2 + k2h - 2.0 * k2_minus_h - 2.0 * kh)


def n4_taylor_coefficients() -> tuple[float, float, float]:
    """Return (c0, c2, c4) in E_min(h) = c0 + c2*h^2 + c4*h^4 + O(h^6) for N=4."""
    k0 = periodic_gaussian_kernel_fourier(0.0, 4.0, order=0)
    k2 = periodic_gaussian_kernel_fourier(2.0, 4.0, order=0)
    k0_2 = periodic_gaussian_kernel_fourier(0.0, 4.0, order=2)
    k2_2 = periodic_gaussian_kernel_fourier(2.0, 4.0, order=2)
    k0_4 = periodic_gaussian_kernel_fourier(0.0, 4.0, order=4)
    k2_4 = periodic_gaussian_kernel_fourier(2.0, 4.0, order=4)

    c0 = 1.0 - 0.5 * k0 + 0.5 * k2
    c2 = -0.5 * (k0_2 - k2_2)
    c4 = -0.5 * ((14.0 / 24.0) * k0_4 - (2.0 / 24.0) * k2_4)
    return c0, c2, c4
