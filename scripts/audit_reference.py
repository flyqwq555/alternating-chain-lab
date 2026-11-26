"""Independent brute-force reference implementation, used only to audit energy_per_particle."""

from __future__ import annotations

import math


def ref_kernel(x: float, period: float, m_max: int = 200) -> float:
    total = 0.0
    for m in range(-m_max, m_max + 1):
        d = x + m * period
        total += math.exp(-0.25 * d * d)
    return total


def ref_energy(positions: list[float]) -> float:
    n = len(positions)
    period = float(n)
    total = 0.0
    for i in range(n):
        si = 1.0 if i % 2 == 0 else -1.0
        for j in range(n):
            sj = 1.0 if j % 2 == 0 else -1.0
            total += si * sj * ref_kernel(positions[i] - positions[j], period)
    return 1.0 - total / n
