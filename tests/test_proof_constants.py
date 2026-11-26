"""Guards every numerical constant quoted in docs/proof-n4-phase-transition.md."""

from __future__ import annotations

import pathlib
import sys

import mpmath as mp

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

import verify_proof as vp  # noqa: E402


def test_proof_constants_match_document() -> None:
    rep = vp.run_all()
    cr, fo, te = rep["crossing"], rep["first_order"], rep["t_elimination"]

    # theorem 2: unique crossing
    assert abs(cr["h_c"] - mp.mpf("0.66801911548686144337")) < mp.mpf("1e-20")
    assert cr["F_prime_at_h_c"] > mp.mpf("0.09")
    assert abs(cr["h_star"] - mp.mpf("0.784532788955013")) < mp.mpf("1e-12")
    assert cr["F_second_at_h_star"] < 0
    assert vp.F(cr["bracket_left"]) < 0 < vp.F(cr["bracket_right"])

    # theorem 4: first order, transition precedes instability
    assert abs(fo["h_instab"] - mp.mpf("0.66846732835982469912")) < mp.mpf("1e-20")
    assert fo["h_c"] < fo["h_instab"]
    assert fo["hessian_at_h_c"] > 0

    # theorem 3: envelope-theorem step and its certified radius
    assert abs(te["Delta"] - mp.mpf("0.300625807303069")) < mp.mpf("1e-14")
    assert te["delta0"] > mp.mpf("0.08")
    assert te["worst_dE_du"] > 0
    assert te["d2E_dt2"] > 0

    # the documented h < 1/(500 pi) regime sits far inside the certified range
    assert te["h0"] > 50 / (500 * mp.pi)


def test_global_small_h_constants() -> None:
    gs = vp.check_global_small_h()

    assert gs["delta_max"] == mp.mpf(1) / 4
    assert gs["T_bound"] <= vp.q / 120
    assert gs["U_bound"] <= vp.q / 24
    assert gs["radius_upper"] < mp.pi / 4
    assert gs["d_coeff"] < mp.mpf(6) / 5
    assert gs["final_margin"] > mp.mpf(73) / 80
    assert gs["counterexample_D"] < 5 * mp.mpf("0.001") ** 2


def test_lemma7_counterexample() -> None:
    gs = vp.check_global_small_h()

    # The old uniform lower bound 0.0749 fails at this legal box point.
    assert gs["box_corner_Ett"] < mp.mpf("0.0749")
    assert gs["corrected_Ett_lower"] > 0
