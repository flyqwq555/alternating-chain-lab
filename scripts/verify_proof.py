"""Machine verification of every numerical constant in docs/proof-n4-phase-transition.md.

Each constant quoted in the note is produced here, and every truncation error is
controlled by an explicit geometric majorant (not estimated numerically), so the
quoted digits are certified rather than merely computed.

Run directly (prints a report, exit code 1 on any failure) or via pytest.
"""

from __future__ import annotations

import mpmath as mp

DPS = 50          # working precision
N_TRUNC = 30      # Fourier truncation; tail majorant is < 1e-30

mp.mp.dps = DPS

q = mp.e ** (-(mp.pi ** 2) / 4)
c = mp.sqrt(mp.pi) / 2
W = mp.pi / 2     # x -> theta scaling factor


def _series_bound(k: int) -> tuple[mp.mpf, mp.mpf]:
    """Return (partial, upper bound) for sum_{n>=1} n^k q^{n^2}.

    Tail: for n=N+1+j, n^2 >= (N+1)n and n <= (N+1)2^j, hence
    n^k q^{n^2} <= (N+1)^k 2^{kj} q^{(N+1)^2} (q^{N+1})^j.  Summing the
    geometric series gives a finite, rigorous majorant.
    """
    part = mp.fsum([q ** (n * n) * n ** k for n in range(1, N_TRUNC + 1)])
    base = q ** ((N_TRUNC + 1) ** 2)
    ratio = q ** (N_TRUNC + 1)
    denominator = 1 - (2 ** k) * ratio
    tail = (N_TRUNC + 1) ** k * base / denominator
    return part, part + tail


def g(theta, k: int = 0):
    """k-th derivative of g(theta) = K(2 theta / pi), period 2*pi."""
    return 2 * c * mp.fsum(
        [q ** (n * n) * n ** k * mp.cos(n * theta + k * mp.pi / 2) for n in range(1, N_TRUNC + 1)]
    )


def Gd(x, k: int = 0):
    """k-th derivative of the periodised Gaussian kernel K(x) = g(pi x / 2)."""
    return W ** k * g(W * x, k)


SUP = {k: 2 * c * _series_bound(k)[1] for k in range(0, 7)}


# ------------------------------------------------------------------- energies

def E_A(h):
    """Symmetric 1+3 configuration (0, 2-h, 2, 2+h)."""
    return 1 - Gd(0) - mp.mpf(1) / 2 * (Gd(2) + Gd(2 * h) - 2 * Gd(2 - h) - 2 * Gd(h))


def E_B(h):
    """Three-gap-locked configuration (0, h, 2h, 3h)."""
    return 1 - Gd(0) - mp.mpf(1) / 2 * (2 * Gd(2 * h) - 3 * Gd(h) - Gd(3 * h))


def F(h):
    """E_A - E_B; the sign says which candidate is lower."""
    return E_A(h) - E_B(h)


def Fp(h):
    return mp.mpf(1) / 2 * (2 * Gd(2 * h, 1) - 2 * Gd(2 - h, 1) - Gd(h, 1) - 3 * Gd(3 * h, 1))


def Fpp(h):
    return mp.mpf(1) / 2 * (4 * Gd(2 * h, 2) + 2 * Gd(2 - h, 2) - Gd(h, 2) - 9 * Gd(3 * h, 2))


# analytic sup-bounds for the Lipschitz margins
M2 = 8 * W ** 2 * SUP[2]      # |F''| <= (1/2)(4+2+1+9) sup|G''|
M3 = 19 * W ** 3 * SUP[3]     # |F'''| <= (1/2)(8+2+1+27) sup|G'''|


# ------------------------------------------------------- 1. the crossing h_c

def check_crossing():
    h_c = mp.findroot(F, mp.mpf("0.668"))
    assert F(h_c - mp.mpf("1e-14")) < 0 < F(h_c + mp.mpf("1e-14"))
    bracket_left = mp.mpf("0.6680191154868614433")
    bracket_right = mp.mpf("0.6680191154868614435")
    assert F(bracket_left) < 0 < F(bracket_right)

    slope = Fp(h_c)
    assert slope > mp.mpf("0.09"), "crossing is not transversal enough"

    h_star = mp.findroot(Fp, mp.mpf("0.79"))
    assert Fpp(h_star) < 0, "h* must be a strict maximum of F"

    # (a) near h = 0: F'(0) = 0 by symmetry, and F'' > 0 on [0, 1e-3]
    steps = 40
    s = mp.mpf("1e-3") / steps
    for i in range(steps + 1):
        assert Fpp(i * s) > M3 * s / 2, "F'' not certified positive near 0"

    # (b) F' > 0 on the bulk [1e-3, h*-rho]
    rho = mp.mpf("2e-3")
    steps = 2500
    a, b = mp.mpf("1e-3"), h_star - rho
    s = (b - a) / steps
    for i in range(steps + 1):
        h = a + i * s
        assert Fp(h) > M2 * s / 2, f"F' not certified positive at h={mp.nstr(h, 12)}"

    # (c) inside [h*-rho, h*+rho]: F'' < 0 => F' strictly decreasing; with opposite
    #     signs at the ends this pins the zero of F' uniquely at h*
    steps = 40
    s = 2 * rho / steps
    for i in range(steps + 1):
        h = h_star - rho + i * s
        assert Fpp(h) < -M3 * s / 2, f"F'' not certified negative at h={mp.nstr(h, 12)}"
    assert Fp(h_star - rho) > 0 > Fp(h_star + rho)

    # (d) F' < 0 on [h*+rho, 1-eta]
    eta = mp.mpf("1e-3")
    steps = 2000
    a, b = h_star + rho, 1 - eta
    s = (b - a) / steps
    for i in range(steps + 1):
        h = a + i * s
        assert Fp(h) < -M2 * s / 2, f"F' not certified negative at h={mp.nstr(h, 12)}"

    # (e) near h = 1: F'(1) = 0 exactly, and F'' > 0 on [1-eta, 1]
    assert abs(Fp(mp.mpf(1))) < mp.mpf("1e-40"), "F'(1) must vanish"
    steps = 40
    s = eta / steps
    for i in range(steps + 1):
        assert Fpp(1 - i * s) > M3 * s / 2, "F'' not certified positive near 1"

    # conclusions: F < 0 on (0,h_c), F > 0 on (h_c,1), F(1) = 0
    assert F(mp.mpf("0.5")) < 0 and F(mp.mpf("0.9")) > 0
    assert abs(F(mp.mpf(1))) < mp.mpf("1e-40")

    return {
        "h_c": h_c, "F_prime_at_h_c": slope, "h_star": h_star,
        "F_second_at_h_star": Fpp(h_star), "F_at_1": F(mp.mpf(1)),
        "bracket_left": bracket_left, "bracket_right": bracket_right,
        "M2": M2, "M3": M3,
    }


# ------------------------------------------ 2. first-order / instability data

def check_first_order():
    h_c = mp.findroot(F, mp.mpf("0.668"))
    h_instab = mp.findroot(lambda h: Gd(2, 2) - 2 * Gd(2 - h, 2), mp.mpf("0.6685"))
    assert h_c < h_instab, "the crossing must precede the Hessian instability"

    hess = -mp.mpf(1) / 2 * (Gd(2, 2) - 2 * Gd(2 - h_c, 2))
    assert hess > 0, "the symmetric branch must stay locally stable at the crossing"
    return {"h_c": h_c, "h_instab": h_instab, "hessian_at_h_c": hess, "gap": h_instab - h_c}


# ------------------------------------------------- 3. t-elimination constants

def check_t_elimination():
    """Explicit constants for the envelope-theorem (t-elimination) step.

    Angular variables, period 2*pi, delta = pi*h/2, positions 0, pi+t-u, pi+t,
    pi+t+v with u, v >= delta. Box used for the certificate: u, v in [delta, 3 delta],
    |t| <= 4 delta.
    """
    Delta = g(mp.pi, 2) - g(mp.mpf(0), 2)          # g''(pi) - g''(0) > 0
    assert Delta > 0
    g2_pi = g(mp.pi, 2)
    assert g2_pi > 0

    kappa, tmax = 3, 4
    s5 = SUP[5]
    g4_0, g4_pi = abs(g(mp.mpf(0), 4)), abs(g(mp.pi, 4))

    cubic_1 = 3 * kappa ** 3 + 3 * kappa ** 3 + kappa ** 3   # |3u^2v+3uv^2+v^3| / delta^3
    cubic_2 = (tmax + kappa) ** 3                            # |t-u|^3 / delta^3
    rem = 3 * (s5 / 24) * (tmax + kappa) ** 4                # Taylor remainder / delta^4

    gross = tmax ** 3 + (tmax + kappa) ** 3 + (tmax + kappa) ** 3
    eta_cubic = g4_pi / 6 * gross / g2_pi
    eta_quartic = rem / g2_pi

    b3 = g2_pi * eta_cubic + g4_0 / 6 * cubic_1 + g4_pi / 6 * cubic_2
    b4 = g2_pi * eta_quartic + rem

    def positive(dl):
        return Delta / 2 - b3 * dl ** 2 / 2 - b4 * dl ** 3 / 2 > 0

    lo, hi = mp.mpf("1e-6"), mp.mpf("0.5")
    assert positive(lo) and not positive(hi)
    for _ in range(200):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if positive(mid) else (lo, mid)
    root_estimate = lo
    delta0 = mp.mpf("0.0868")
    assert positive(delta0), "fixed delta_0=0.0868 must have strict margin"
    assert delta0 > mp.mpf("1e-3"), "certified range must contain the h < 1/(500 pi) regime"

    box_lower = c * q * (
        mp.mpf(119) / 120 - mp.mpf(475) / 8 * delta0 ** 2
    )
    assert box_lower > 0

    # independent grid check of dE/du > 0 (no algebra reused)
    def dE_du(u, v, t):
        return -mp.mpf(1) / 2 * (g(u + v, 1) + g(mp.pi + t - u, 1) - g(u, 1))

    worst = None
    for dl in (mp.mpf("1e-3"), mp.mpf("0.01"), delta0 * mp.mpf("0.9")):
        for i in range(7):
            for j in range(7):
                u = dl * (1 + 2 * mp.mpf(i) / 6)
                v = dl * (1 + 2 * mp.mpf(j) / 6)
                for k in range(7):
                    t = (u - v) + (2 * mp.mpf(k) / 6 - 1) * dl
                    val = dE_du(u, v, t)
                    if worst is None or val < worst[0]:
                        worst = (val, dl, u, v, t)
    assert worst[0] > 0, "grid check of dE/du failed"

    # dE/dt = 0 at t = 0 when u = v (oddness of g' about pi), and strict convexity
    dl = mp.mpf("1e-3")
    dt0 = -mp.mpf(1) / 2 * (g(mp.pi, 1) - g(mp.pi - dl, 1) - g(mp.pi + dl, 1))
    assert abs(dt0) < mp.mpf("1e-30"), "dE/dt must vanish at t=0 when u=v"
    d2t = -mp.mpf(1) / 2 * (g(mp.pi, 2) - g(mp.pi - dl, 2) - g(mp.pi + dl, 2))
    assert d2t > 0, "strict convexity in t"

    return {
        "Delta": Delta, "g2_pi": g2_pi, "sup_g4": SUP[4], "sup_g5": s5,
        "b3": b3, "b4": b4, "root_estimate": root_estimate,
        "delta0": delta0, "h0": 2 * delta0 / mp.pi,
        "box_Ett_lower": box_lower,
        "worst_dE_du": worst[0], "d2E_dt2": d2t,
    }


# ------------------------------------------------ 4. global small-h theorem

def check_global_small_h():
    """Constants and regression checks for Theorem 5 (0 < h <= 1/(2*pi))."""
    delta_max = mp.mpf(1) / 4

    T = _series_bound(2)[1] - q
    U = _series_bound(4)[1] - q
    assert T <= q / 120, "T bound failed"
    assert U <= q / 24, "U bound failed"

    # g''' is positive on (0, pi); check the analytic lower bound on a grid.
    for i in range(1, 200):
        theta = mp.pi * i / 200
        lower = 2 * c * (q - U) * mp.sin(theta)
        assert g(theta, 3) >= lower - mp.mpf("1e-40")

    # Localization radius at delta=1/4 and the refined defect coefficient.
    radius_upper = (1 + 3 * mp.pi / 5) * delta_max
    assert radius_upper < mp.pi / 4
    d_coeff = (mp.mpf(121) / 60) / (
        mp.mpf(37) / 20 - mp.mpf(5) / 2 * delta_max ** 2
    )
    assert d_coeff <= mp.mpf(968) / 813
    assert d_coeff < mp.mpf(6) / 5

    # Trial-energy and final derivative constants.
    trial_coeff = 2 * (q + T) / q
    assert trial_coeff <= mp.mpf(121) / 60
    final_margin = 2 * mp.cos(delta_max) - 1 - 3 * T / q
    assert final_margin >= mp.mpf(39) / 40 - delta_max ** 2
    assert final_margin > mp.mpf(73) / 80

    # Regression: the old uniform E_tt >= 0.0749 claim is false at this box point.
    delta = mp.mpf("0.0868")
    box_corner = -(g(mp.pi + 4 * delta, 2)
                   - g(mp.pi + 3 * delta, 2)
                   - g(mp.pi + 5 * delta, 2)) / 2
    corrected_lower = c * q * (
        mp.mpf(119) / 120 - mp.mpf(475) / 8 * delta ** 2
    )
    assert box_corner < mp.mpf("0.0749")
    assert corrected_lower > 0

    # A legal configuration with D < 5*delta^2 but u=4*delta > 3*delta.
    counter_delta = mp.mpf("0.001")
    counter_s1 = (1 + mp.exp(-1j * counter_delta)
                  - mp.exp(3j * counter_delta)
                  + mp.exp(4j * counter_delta))
    counter_D = 2 - abs(counter_s1)
    assert counter_D < 5 * counter_delta ** 2

    return {
        "delta_max": delta_max, "h_max": 2 * delta_max / mp.pi,
        "T_bound": T, "U_bound": U, "radius_upper": radius_upper,
        "d_coeff": d_coeff, "trial_coeff": trial_coeff,
        "final_margin": final_margin, "box_corner_Ett": box_corner,
        "corrected_Ett_lower": corrected_lower,
        "counterexample_D": counter_D,
    }


# ----------------------------------------------------------------- reporting

def run_all():
    assert q < mp.mpf(1) / 8
    return {
        "kernel": {"q": q, "q_lt_one_eighth": True},
        "crossing": check_crossing(),
        "first_order": check_first_order(),
        "t_elimination": check_t_elimination(),
        "global_small_h": check_global_small_h(),
    }


def main() -> int:
    rep = run_all()
    k, cr, fo, te, gs = (
        rep["kernel"], rep["crossing"], rep["first_order"],
        rep["t_elimination"], rep["global_small_h"],
    )

    print("=" * 68)
    print(" verification of docs/proof-n4-phase-transition.md")
    print("=" * 68)
    print(f" q = {mp.nstr(k['q'], 18)}   (< 1/8)")
    print(" Fourier truncation majorant < 1e-30\n")

    print(f" h_c        = {mp.nstr(cr['h_c'], 20)}")
    print(f" F'(h_c)    = {mp.nstr(cr['F_prime_at_h_c'], 15)}")
    print(f" h*         = {mp.nstr(cr['h_star'], 15)}   F''(h*) = {mp.nstr(cr['F_second_at_h_star'], 10)}")
    print(f" |F''| <= {mp.nstr(cr['M2'], 8)}    |F'''| <= {mp.nstr(cr['M3'], 8)}\n")

    print(f" h_instab       = {mp.nstr(fo['h_instab'], 20)}")
    print(f" h_instab - h_c = {mp.nstr(fo['gap'], 8)}")
    print(f" Hessian at h_c = {mp.nstr(fo['hessian_at_h_c'], 12)}  (> 0)\n")

    print(f" Delta = g''(pi)-g''(0) = {mp.nstr(te['Delta'], 15)}")
    print(f" sup|g''''| <= {mp.nstr(te['sup_g4'], 10)}   sup|g^(5)| <= {mp.nstr(te['sup_g5'], 10)}")
    print(f" b3 = {mp.nstr(te['b3'], 10)}   b4 = {mp.nstr(te['b4'], 10)}")
    print(f" delta_0 = {mp.nstr(te['delta0'], 10)}  =>  h_0 = {mp.nstr(te['h0'], 10)}")
    print(f" delta_0 root estimate = {mp.nstr(te['root_estimate'], 18)}")
    print(f" corrected box E_tt lower = {mp.nstr(te['box_Ett_lower'], 10)}  (> 0)")
    print(f" 1/(500 pi) = {mp.nstr(1 / (500 * mp.pi), 10)}  (far inside the certified range)")
    print(f" min dE/du on grid = {mp.nstr(te['worst_dE_du'], 10)}  (> 0)")
    print(f" d2E/dt2 at (delta,delta,0) = {mp.nstr(te['d2E_dt2'], 10)}  (> 0)")
    print(f" theorem 5 range h <= {mp.nstr(gs['h_max'], 12)}")
    print(f" radius at delta=1/4 = {mp.nstr(gs['radius_upper'], 12)}  (< pi/4)")
    print(f" refined D coefficient = {mp.nstr(gs['d_coeff'], 12)}  (< 6/5)")
    print(f" final Et/sin(t) margin = {mp.nstr(gs['final_margin'], 12)}  (> 73/80)")
    print("\n ALL CHECKS PASSED")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
