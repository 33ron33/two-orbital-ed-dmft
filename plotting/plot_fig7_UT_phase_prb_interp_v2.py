#!/usr/bin/env python3
"""
PRB-style interpolated U-T phase diagram for Figure 7.

Updates in v2:
- removes blank low-temperature region by setting y-limits to actual Tmin..Tmax
- capitalizes axis labels: Temperature, Interaction strength
- keeps smooth interpolated phase map with no black grid lines
- right panel shows Z1, Z2 at the lowest temperature / largest beta
"""

import argparse
import glob
from pathlib import Path
import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.patches import Patch
from matplotlib.ticker import MultipleLocator, FormatStrFormatter
from scipy.interpolate import griddata

# ============================================================
# PRB-LIKE STYLE
# ============================================================
mpl.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "text.latex.preamble": r"\usepackage{amsmath}",
    "axes.labelsize": 22,
    "font.size": 20,
    "legend.fontsize": 15,
    "xtick.labelsize": 18,
    "ytick.labelsize": 18,
})
mpl.rcParams["axes.xmargin"] = 0.02
mpl.rcParams["axes.ymargin"] = 0.04


def load_rows(data_dir, jr_target, jr_tol=2e-3):
    rows = []
    files = sorted(glob.glob(f"{data_dir}/*.npz"))

    for f in files:
        d = np.load(f)

        U = float(d["U"])
        J = float(d["J"])
        beta = float(d["beta"])
        T = 1.0 / beta

        if U > 1e-12:
            Jr = J / U
            if abs(Jr - jr_target) > jr_tol:
                continue
        else:
            Jr = jr_target

        rows.append(dict(
            file=f,
            U=U,
            J=J,
            Jr=Jr,
            beta=beta,
            T=T,
            Z1=float(d["Z1"]),
            Z2=float(d["Z2"]),
            n1=float(d["n1"]),
            n2=float(d["n2"]),
            conv=bool(d["converged"]),
            delta=float(d["last_delta"]),
        ))

    if not rows:
        raise RuntimeError(f"No matching rows found in {data_dir}")

    return rows


def classify_phase(Z1, Z2, Zc):
    phase = np.full_like(Z1, np.nan, dtype=float)

    both_metal = (Z1 > Zc) & (Z2 > Zc)
    both_mott = (Z1 <= Zc) & (Z2 <= Zc)
    osmt = np.isfinite(Z1) & np.isfinite(Z2) & ~(both_metal | both_mott)

    phase[both_metal] = 0
    phase[osmt] = 1
    phase[both_mott] = 2

    return phase


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/nb1_scanD_fig7_beta5_100_dense")
    ap.add_argument("--out", default="figures_nb1_fig7_beta5_100_interp_v2")
    ap.add_argument("--t2", type=float, default=0.3)
    ap.add_argument("--Jr", type=float, default=0.10)
    ap.add_argument("--z_thresh", type=float, default=0.15)
    ap.add_argument("--nx", type=int, default=650)
    ap.add_argument("--ny", type=int, default=450)
    ap.add_argument("--show_points", action="store_true")
    args = ap.parse_args()

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = load_rows(args.data, args.Jr)

    U_raw = np.array([r["U"] for r in rows])
    T_raw = np.array([r["T"] for r in rows])
    Z1_raw = np.array([r["Z1"] for r in rows])
    Z2_raw = np.array([r["Z2"] for r in rows])

    U_min, U_max = float(np.min(U_raw)), float(np.max(U_raw))
    T_min, T_max = float(np.min(T_raw)), float(np.max(T_raw))

    U_fine = np.linspace(U_min, U_max, args.nx)
    T_fine = np.linspace(T_min, T_max, args.ny)
    UU, TT = np.meshgrid(U_fine, T_fine)

    points = np.column_stack([U_raw, T_raw])

    Z1_lin = griddata(points, Z1_raw, (UU, TT), method="linear")
    Z2_lin = griddata(points, Z2_raw, (UU, TT), method="linear")

    Z1_near = griddata(points, Z1_raw, (UU, TT), method="nearest")
    Z2_near = griddata(points, Z2_raw, (UU, TT), method="nearest")

    Z1_grid = np.where(np.isfinite(Z1_lin), Z1_lin, Z1_near)
    Z2_grid = np.where(np.isfinite(Z2_lin), Z2_lin, Z2_near)

    phase = classify_phase(Z1_grid, Z2_grid, args.z_thresh)

    cmap = ListedColormap([
        "#3B8FD8",   # both metal
        "#F0B000",   # OSMT
        "#D94A4A",   # both Mott
    ])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5], cmap.N)

    fig, (ax0, ax1) = plt.subplots(
        1, 2,
        figsize=(12.8, 5.0),
        dpi=300,
        gridspec_kw={"width_ratios": [1.08, 1.0], "wspace": 0.28}
    )

    # ========================================================
    # Left panel: interpolated U-T phase map
    # ========================================================
    ax0.pcolormesh(
        UU, TT, phase,
        cmap=cmap,
        norm=norm,
        shading="auto",
        linewidth=0.0,
        edgecolors="none",
        antialiased=True,
        zorder=1,
    )

    if args.show_points:
        ax0.scatter(
            U_raw, T_raw,
            s=5,
            c="black",
            alpha=0.15,
            linewidths=0,
            zorder=3,
        )

    ax0.set_xlabel(r"Interaction strength $(U/t_1)$")
    ax0.set_ylabel(r"Temperature $(T/t_1 = 1/\beta)$")
    ax0.set_title(rf"$J/U = {args.Jr:.2f}$, $Z_c = {args.z_thresh:.2f}$", pad=10)

    ax0.set_xlim(U_min, U_max)

    # Important: remove blank strip below lowest computed temperature.
    ax0.set_ylim(T_min, T_max)

    ax0.xaxis.set_major_locator(MultipleLocator(1.0))
    ax0.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax0.yaxis.set_major_locator(MultipleLocator(0.025))
    ax0.yaxis.set_minor_locator(MultipleLocator(0.0125))
    ax0.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax0.yaxis.set_major_formatter(FormatStrFormatter("%.3f"))

    legend_handles = [
        Patch(facecolor="#3B8FD8", edgecolor="black", label=r"Both Metal"),
        Patch(facecolor="#F0B000", edgecolor="black", label=r"OSMT"),
        Patch(facecolor="#D94A4A", edgecolor="black", label=r"Both Insulator"),
    ]

    ax0.legend(
        handles=legend_handles,
        frameon=False,
        loc="upper left",
        handlelength=2.1,
        borderpad=0.3,
        labelspacing=0.45,
    )

    # ========================================================
    # Right panel: lowest-temperature Z(U) cut
    # ========================================================
    beta_max = max(r["beta"] for r in rows)
    cut = sorted(
        [r for r in rows if abs(r["beta"] - beta_max) < 1e-8],
        key=lambda r: r["U"]
    )

    U_cut = np.array([r["U"] for r in cut])
    Z1_cut = np.array([r["Z1"] for r in cut])
    Z2_cut = np.array([r["Z2"] for r in cut])

    ax1.plot(
        U_cut, Z1_cut,
        marker="o",
        lw=2.4,
        ms=6.5,
        label=r"$Z_1$",
    )
    ax1.plot(
        U_cut, Z2_cut,
        marker="s",
        lw=2.4,
        ms=6.5,
        label=r"$Z_2$",
    )
    ax1.axhline(
        args.z_thresh,
        color="0.45",
        lw=1.8,
        ls="--",
        label=r"$Z_c$",
    )

    ax1.set_xlabel(r"Interaction strength $(U/t_1)$")
    ax1.set_ylabel(r"$Z_m$")
    ax1.set_title(rf"Lowest temperature: $\beta = {beta_max:.0f}$", pad=10)

    ax1.set_xlim(U_min, U_max)
    ax1.set_ylim(-0.04, 1.02)

    ax1.xaxis.set_major_locator(MultipleLocator(1.0))
    ax1.xaxis.set_minor_locator(MultipleLocator(0.25))
    ax1.yaxis.set_major_locator(MultipleLocator(0.2))
    ax1.yaxis.set_minor_locator(MultipleLocator(0.1))
    ax1.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
    ax1.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    ax1.legend(
        frameon=False,
        loc="upper right",
        handlelength=2.2,
        borderpad=0.3,
        labelspacing=0.45,
    )

    # ========================================================
    # Common PRB tick/spine styling
    # ========================================================
    for ax in (ax0, ax1):
        for spine in ax.spines.values():
            spine.set_linewidth(1.7)

        ax.tick_params(
            which="major",
            direction="in",
            top=True,
            right=True,
            length=8,
            width=1.4,
            pad=8,
        )
        ax.tick_params(
            which="minor",
            direction="in",
            top=True,
            right=True,
            length=4.5,
            width=1.0,
        )

    fig.tight_layout()

    pdf = out_dir / "fig7_UT_phase_interp_prb_v2.pdf"
    png = out_dir / "fig7_UT_phase_interp_prb_v2.png"
    csv = out_dir / "fig7_UT_phase_summary_v2.csv"

    fig.savefig(pdf, dpi=1200, bbox_inches="tight")
    fig.savefig(png, dpi=400, bbox_inches="tight")

    with open(csv, "w") as f:
        f.write("U,beta,T,J_over_U,J,Z1,Z2,n1,n2,converged,last_delta,file\n")
        for r in sorted(rows, key=lambda x: (x["beta"], x["U"])):
            f.write(
                f"{r['U']},{r['beta']},{r['T']},{r['Jr']},{r['J']},"
                f"{r['Z1']},{r['Z2']},{r['n1']},{r['n2']},"
                f"{r['conv']},{r['delta']},{r['file']}\n"
            )

    print(f"[saved] {pdf}")
    print(f"[saved] {png}")
    print(f"[saved] {csv}")
    print(f"Loaded {len(rows)} records from {args.data}")
    print(f"T range plotted: {T_min:.5f} to {T_max:.5f}")


if __name__ == "__main__":
    main()
