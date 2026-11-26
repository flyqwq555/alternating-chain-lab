"""Reproduce all numerical phase diagrams, multi-N comparisons, and rigorous checks."""

from __future__ import annotations

import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

import matplotlib.pyplot as plt

from alternating_chain.energy import n4_taylor_coefficients
from alternating_chain.multi_particle import scan_multi_n
from alternating_chain.phase_diagram import scan_n4_phase_diagram
from alternating_chain.rigour import run_rigorous_checks


def main() -> None:
    fig_dir = ROOT / "figures"
    fig_dir.mkdir(exist_ok=True)

    print("=== 1. Rigorous mpmath / Interval Arithmetic Checks ===")
    rep = run_rigorous_checks(dps=55)
    print(f"All rigorous checks passed: {rep.all_passed}")
    print(f"  c0 = {rep.c0_str}")
    print(f"  c2 = {rep.c2_str}")
    print(f"  c4 = {rep.c4_str}")
    print(f"  h_cross  = {rep.h_c1_str}")
    print(f"  h_instab = {rep.h_c2_str}")

    c0, c2, c4 = n4_taylor_coefficients()
    print(f"Float Taylor check: c0={c0:.15f}, c2={c2:.15f}, c4={c4:.15f}")

    print("\n=== 2. N=4 Full Phase Diagram h in (0, 1] ===")
    pts, (hc1, hc2) = scan_n4_phase_diagram(h_min=0.01, h_max=1.0, num_points=80)
    print(f"Critical points: h_cross = {hc1:.12f}, h_instab = {hc2:.12f}")

    hs = [p.h for p in pts]
    e_best = [p.best_energy for p in pts]
    e_sym = [p.symmetric_energy for p in pts]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4.8))
    ax1.plot(hs, e_sym, "--", label="Symmetric (0, 2-h, 2, 2+h)", color="tab:orange")
    ax1.plot(hs, e_best, "-", label="Global Minimum E_min(h)", color="tab:blue", lw=2)
    ax1.axvline(hc1, color="tab:red", ls=":", label=f"h_cross = {hc1:.4f}")
    ax1.set_xlabel("Minimum spacing h")
    ax1.set_ylabel("Energy per particle E")
    ax1.set_title("N=4 Ground State Energy vs h")
    ax1.legend(fontsize=8)
    ax1.grid(alpha=0.3)

    g_sorted = [sorted(p.gaps) for p in pts]
    for idx in range(4):
        ax2.plot(hs, [g[idx] for g in g_sorted], label=f"gap[{idx+1}]")
    ax2.axvline(hc1, color="tab:red", ls=":")
    ax2.axvline(hc2, color="tab:green", ls=":")
    ax2.set_xlabel("Minimum spacing h")
    ax2.set_ylabel("Sorted periodic gaps")
    ax2.set_title("N=4 Optimal Gap Bifurcation")
    ax2.legend(fontsize=8)
    ax2.grid(alpha=0.3)

    fig.tight_layout()
    fig.savefig(fig_dir / "n4_phase_diagram.png", dpi=160)
    plt.close(fig)

    print("\n=== 3. Multi-N (N=4, 6, 8, 10) Comparison ===")
    multi = scan_multi_n(n_values=(4, 6, 8, 10), h_values=(0.05, 0.2, 0.5, 0.75, 1.0))
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    for n in (4, 6, 8, 10):
        sub = [r for r in multi if r.n == n]
        ax.plot([r.h for r in sub], [r.best_energy for r in sub], "o-", label=f"N={n}")
        for r in sub:
            print(f"N={r.n:2d}, h={r.h:.2f}: best_E={r.best_energy:.9f}, eq_E={r.equidistant_energy:.9f}")
    ax.set_xlabel("Minimum spacing h")
    ax.set_ylabel("Minimal energy per particle")
    ax.set_title("Minimal Energy across N=4, 6, 8, 10")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(fig_dir / "multi_n_comparison.png", dpi=160)
    plt.close(fig)
    print(f"\nSaved figures to {fig_dir}")


if __name__ == "__main__":
    main()
