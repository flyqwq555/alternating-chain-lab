"""High-precision mpmath and interval-arithmetic verification of all analytical inequalities."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict

import mpmath as mp
from mpmath import iv


@dataclass
class RigorousReport:
    q_less_than_one_eighth: bool
    tail_bound_holds: bool
    stability_constant_holds: bool
    envelope_second_deriv_positive: bool
    c0_str: str
    c2_str: str
    c4_str: str
    h_c1_str: str
    h_c2_str: str
    details: Dict[str, str]

    @property
    def all_passed(self) -> bool:
        return (
            self.q_less_than_one_eighth
            and self.tail_bound_holds
            and self.stability_constant_holds
            and self.envelope_second_deriv_positive
        )


def _mp_k4(x: mp.mpf, order: int = 0, n_max: int = 40) -> mp.mpf:
    c = mp.sqrt(mp.pi) / 2
    w0 = mp.pi / 2
    q = mp.exp(-(mp.pi ** 2) / 4)
    if order == 0:
        s = mp.mpf(1) + 2 * mp.fsum((q ** (n * n)) * mp.cos(w0 * n * x) for n in range(1, n_max + 1))
        return c * s
    s = mp.mpf(0)
    for n in range(1, n_max + 1):
        wn = w0 * n
        coeff = 2 * (q ** (n * n)) * (wn ** order)
        mod = order % 4
        if mod == 1:
            s -= coeff * mp.sin(wn * x)
        elif mod == 2:
            s -= coeff * mp.cos(wn * x)
        elif mod == 3:
            s += coeff * mp.sin(wn * x)
        else:
            s += coeff * mp.cos(wn * x)
    return c * s


def run_rigorous_checks(dps: int = 55) -> RigorousReport:
    """Verify all proof inequalities with interval arithmetic and 55-digit precision."""
    mp.mp.dps = dps
    iv.dps = max(40, dps - 5)

    q_iv = iv.exp(-(iv.pi ** 2) / 4)
    q_less_than_one_eighth = bool(q_iv < iv.mpf(1) / 8)

    finite_tail_iv = sum((n * n) * (q_iv ** (n * n)) for n in range(2, 13))
    rem_iv = iv.mpf(["0", "1e-120"])
    tail_iv = finite_tail_iv + rem_iv

    mid_bound_iv = 4 * (q_iv ** 4) / (1 - 4 * (q_iv ** 2))
    q_over_120_iv = q_iv / 120
    tail_bound_holds = bool(tail_iv <= mid_bound_iv) and bool(mid_bound_iv <= q_over_120_iv)

    lhs_iv = 2 * q_iv - 36 * tail_iv
    rhs_iv = iv.mpf(17) / 10 * q_iv
    stability_constant_holds = bool(lhs_iv >= rhs_iv)

    c_iv = iv.sqrt(iv.pi) / 2
    err_iv = iv.mpf(["-1e-120", "1e-120"])
    g2_pi_iv = -2 * c_iv * sum(((-1) ** n) * n * n * (q_iv ** (n * n)) for n in range(1, 13)) + err_iv
    g2_0_iv = -2 * c_iv * sum(n * n * (q_iv ** (n * n)) for n in range(1, 13)) + err_iv
    diff_g2_iv = g2_pi_iv - g2_0_iv
    envelope_second_deriv_positive = bool(diff_g2_iv > 0)

    c0 = 1 - mp.mpf("0.5") * _mp_k4(mp.mpf(0), 0) + mp.mpf("0.5") * _mp_k4(mp.mpf(2), 0)
    c2 = -mp.mpf("0.5") * (_mp_k4(mp.mpf(0), 2) - _mp_k4(mp.mpf(2), 2))
    c4 = -mp.mpf("0.5") * (mp.mpf(14) / 24 * _mp_k4(mp.mpf(0), 4) - mp.mpf(2) / 24 * _mp_k4(mp.mpf(2), 4))

    def _e_sym_mp(h: mp.mpf) -> mp.mpf:
        return 1 - _mp_k4(mp.mpf(0), 0) - mp.mpf("0.5") * (
            _mp_k4(mp.mpf(2), 0) + _mp_k4(2 * h, 0) - 2 * _mp_k4(2 - h, 0) - 2 * _mp_k4(h, 0)
        )

    def _e_3gap_mp(h: mp.mpf) -> mp.mpf:
        return 1 - _mp_k4(mp.mpf(0), 0) - mp.mpf("0.5") * (
            2 * _mp_k4(2 * h, 0) - 3 * _mp_k4(h, 0) - _mp_k4(3 * h, 0)
        )

    h_c1 = mp.findroot(lambda h: _e_sym_mp(h) - _e_3gap_mp(h), mp.mpf("0.668"))
    h_c2 = mp.findroot(lambda h: _mp_k4(mp.mpf(2), 2) - 2 * _mp_k4(2 - h, 2), mp.mpf("0.6685"))

    return RigorousReport(
        q_less_than_one_eighth=q_less_than_one_eighth,
        tail_bound_holds=tail_bound_holds,
        stability_constant_holds=stability_constant_holds,
        envelope_second_deriv_positive=envelope_second_deriv_positive,
        c0_str=mp.nstr(c0, 35),
        c2_str=mp.nstr(c2, 35),
        c4_str=mp.nstr(c4, 35),
        h_c1_str=mp.nstr(h_c1, 35),
        h_c2_str=mp.nstr(h_c2, 35),
        details={
            "q": mp.nstr(mp.exp(-(mp.pi ** 2) / 4), 20),
            "g2_pi_minus_g2_0": mp.nstr(diff_g2_iv.mid, 20),
        },
    )
