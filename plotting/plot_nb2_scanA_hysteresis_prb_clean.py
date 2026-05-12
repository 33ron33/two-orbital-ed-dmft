#!/usr/bin/env python3
"""
PRB-style plotting for Nb=2 Scan A hysteresis.

This clean version:
  - excludes the trivial U=0 point by default
  - makes the hysteresis Z/D plot more PRB-like
  - makes the diagnostics plot more PRB-like
  - keeps U=0 available only if --include_U0 is passed

Input example:
    data/nb2_scanA_hysteresis/*.npz

Outputs:
    1. fig_scanA_hysteresis_Z_D_prb_clean.pdf
    2. fig_scanA_diagnostics_prb_clean.pdf
    3. scanA_hysteresis_summary_clean.csv
"""

import argparse
import glob
import os
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator, FormatStrFormatter, LogLocator, NullFormatter
from matplotlib.lines import Line2D


# ============================================================
# PRB STYLE
# ============================================================
def set_prb_style(use_tex=True):
    mpl.rcParams.update({
        "text.usetex": use_tex,
        "font.family": "serif",
        "axes.labelsize": 20,
        "font.size": 18,
        "legend.fontsize": 13.5,
        "xtick.labelsize": 16.5,
        "ytick.labelsize": 16.5,
        "axes.linewidth": 1.6,
        "xtick.major.width": 1.35,
        "ytick.major.width": 1.35,
        "xtick.minor.width": 1.0,
        "ytick.minor.width": 1.0,
        "xtick.major.size": 7.5,
        "ytick.major.size": 7.5,
        "xtick.minor.size": 4.2,
        "ytick.minor.size": 4.2,
        "savefig.dpi": 1200,
    })

    if use_tex:
        mpl.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}"

    mpl.rcParams["axes.xmargin"] = 0.02
    mpl.rcParams["axes.ymargin"] = 0.04


def scalar(x):
    arr = np.asarray(x)
    if arr.shape == ():
        return float(arr)
    return float(arr.ravel()[0])


def decode_seed_value(x):
    """
    Robustly decode seed from npz.
    Handles normal strings, byte strings, and 0-d numpy arrays.
    """
    try:
        arr = np.asarray(x)
        if arr.shape == ():
            val = arr.item()
        else:
            val = arr.ravel()[0]

        if isinstance(val, bytes):
            return val.decode("utf-8").lower()

        return str(val).lower()
    except Exception:
        return ""


def infer_seed(filename, d):
    if "seed" in d.files:
        seed = decode_seed_value(d["seed"])
        if "metal" in seed or seed == "met":
            return "metal"
        if "ins" in seed:
            return "insulator"

    base = os.path.basename(filename).lower()

    if "ins" in base:
        return "insulator"
    if "met" in base:
        return "metal"

    return "unknown"


def load_records(
    data_dir,
    beta=None,
    nb=None,
    t2=None,
    j_abs_tol=1e-8,
    include_U0=False,
    u0_tol=1e-12,
):
    files = sorted(glob.glob(os.path.join(data_dir, "*.npz")))
    rows = []
    skipped_U0 = 0

    for f in files:
        try:
            d = np.load(f, allow_pickle=True)
        except Exception:
            continue

        U = scalar(d["U"])
        J = scalar(d["J"])
        b = scalar(d["beta"]) if "beta" in d.files else None
        Nb = int(scalar(d["N_bath"])) if "N_bath" in d.files else None

        # ------------------------------------------------------------
        # Critical fix:
        # omit U=0 by default for hysteresis plots.
        #
        # At U=0 there is only one physical noninteracting solution,
        # so it should not be part of the interacting hysteresis comparison.
        # ------------------------------------------------------------
        if (not include_U0) and U <= u0_tol:
            skipped_U0 += 1
            continue

        # Optional filters
        if beta is not None and b is not None and abs(b - beta) > 1e-8:
            continue

        if nb is not None and Nb is not None and Nb != nb:
            continue

        if abs(J) > j_abs_tol:
            continue

        # t2 may or may not be stored in older files
        if t2 is not None and "t2" in d.files:
            t2_file = scalar(d["t2"])
            if abs(t2_file - t2) > 1e-8:
                continue

        seed = infer_seed(f, d)

        Z1 = scalar(d["Z1"])
        Z2 = scalar(d["Z2"])
        D1 = scalar(d["D1"])
        D2 = scalar(d["D2"])
        n1 = scalar(d["n1"])
        n2 = scalar(d["n2"])

        rows.append(dict(
            file=f,
            U=U,
            J=J,
            beta=b,
            Nb=Nb,
            seed=seed,
            Z1=Z1,
            Z2=Z2,
            Zavg=0.5 * (Z1 + Z2),
            D1=D1,
            D2=D2,
            Davg=0.5 * (D1 + D2),
            n1=n1,
            n2=n2,
            navg=0.5 * (n1 + n2),
            converged=bool(d["converged"]) if "converged" in d.files else True,
            last_delta=scalar(d["last_delta"]) if "last_delta" in d.files else np.nan,
            n_iter=int(scalar(d["n_iter"])) if "n_iter" in d.files else -1,
        ))

    if not rows:
        raise RuntimeError(f"No matching records found in {data_dir}")

    rows = sorted(rows, key=lambda r: (r["seed"], r["U"]))

    print("=" * 90)
    print(f"Loaded records after filtering : {len(rows)}")
    print(f"Skipped U=0 files             : {skipped_U0}")
    print("=" * 90)

    return rows


def branch(rows, seed):
    rr = [r for r in rows if r["seed"] == seed]
    return sorted(rr, key=lambda r: r["U"])


def common_axis_style(ax):
    for spine in ax.spines.values():
        spine.set_linewidth(1.6)

    ax.tick_params(
        axis="both",
        which="major",
        direction="in",
        top=True,
        right=True,
        length=7.5,
        width=1.35,
        pad=7,
    )
    ax.tick_params(
        axis="both",
        which="minor",
        direction="in",
        top=True,
        right=True,
        length=4.2,
        width=1.0,
    )


def add_panel_label(ax, label):
    ax.text(
        0.025,
        0.955,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )


def plot_hysteresis(rows, out_dir, z_thresh=0.05, beta=25, nb=2):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    met = branch(rows, "metal")
    ins = branch(rows, "insulator")

    if not met and not ins:
        raise RuntimeError("Could not find metal/insulator seeded branches.")

    styles = {
        "metal": dict(color="#0072B2", marker="o", label=r"metal seed"),
        "insulator": dict(color="#D55E00", marker="s", label=r"insulator seed"),
    }

    fig, (axZ, axD) = plt.subplots(
        1,
        2,
        figsize=(11.2, 4.55),
        dpi=300,
        gridspec_kw={"wspace": 0.26},
    )

    for seed, rr in [("metal", met), ("insulator", ins)]:
        if not rr:
            continue

        U = np.array([r["U"] for r in rr])
        Z = np.array([r["Zavg"] for r in rr])
        D = np.array([r["Davg"] for r in rr])

        st = styles[seed]

        axZ.plot(
            U, Z,
            lw=2.45,
            color=st["color"],
            marker=st["marker"],
            ms=7.2,
            mec=st["color"],
            mfc=st["color"],
            mew=1.1,
            label=st["label"],
            zorder=4,
        )

        axD.plot(
            U, D,
            lw=2.45,
            color=st["color"],
            marker=st["marker"],
            ms=7.2,
            mec=st["color"],
            mfc=st["color"],
            mew=1.1,
            label=st["label"],
            zorder=4,
        )

    # Z threshold
    axZ.axhline(
        z_thresh,
        color="0.45",
        lw=1.65,
        ls="--",
        zorder=1,
        label=rf"$Z_c={z_thresh:.2f}$",
    )

    # Panel labels
    add_panel_label(axZ, r"(a)")
    add_panel_label(axD, r"(b)")

    # Titles
    axZ.set_title(r"$J=0$, degenerate bands", pad=8)
    axD.set_title(r"Double occupancy", pad=8)

    # Labels
    axZ.set_xlabel(r"Interaction strength $(U/t)$")
    axD.set_xlabel(r"Interaction strength $(U/t)$")

    axZ.set_ylabel(r"Quasiparticle weight $Z$")
    axD.set_ylabel(r"Double occupancy $D$", labelpad=7)

    allU = np.array([r["U"] for r in rows])
    xmin = np.min(allU) - 0.15
    xmax = np.max(allU) + 0.15

    axZ.set_xlim(xmin, xmax)
    axD.set_xlim(xmin, xmax)

    axZ.set_ylim(-0.04, 1.05)

    # Tight but clean D limits
    allD = np.array([r["Davg"] for r in rows])
    ymin = max(0.0, np.min(allD) - 0.008)
    ymax = min(0.27, np.max(allD) + 0.012)
    axD.set_ylim(ymin, ymax)

    for ax in (axZ, axD):
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        common_axis_style(ax)

    axZ.yaxis.set_major_locator(MultipleLocator(0.2))
    axZ.yaxis.set_minor_locator(MultipleLocator(0.1))
    axZ.yaxis.set_major_formatter(FormatStrFormatter("%.1f"))

    axD.yaxis.set_major_locator(MultipleLocator(0.02))
    axD.yaxis.set_minor_locator(MultipleLocator(0.01))
    axD.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    # One compact legend
    handles = [
        Line2D(
            [0], [0],
            color=styles["metal"]["color"],
            marker="o",
            lw=2.45,
            ms=7.2,
            label=r"metal seed",
        ),
        Line2D(
            [0], [0],
            color=styles["insulator"]["color"],
            marker="s",
            lw=2.45,
            ms=7.2,
            label=r"insulator seed",
        ),
        Line2D(
            [0], [0],
            color="0.45",
            lw=1.65,
            ls="--",
            label=rf"$Z_c={z_thresh:.2f}$",
        ),
    ]

    axZ.legend(
        handles=handles,
        frameon=False,
        loc="upper right",
        handlelength=2.35,
        borderpad=0.2,
        labelspacing=0.38,
    )

    fig.suptitle(
        rf"$N_b={nb}$ ED-DMFT hysteresis, $\beta={beta:.0f}$",
        y=1.02,
        fontsize=20,
    )

    fig.tight_layout()

    pdf = out_dir / "fig_scanA_hysteresis_Z_D_prb_clean.pdf"
    png = out_dir / "fig_scanA_hysteresis_Z_D_prb_clean.png"

    fig.savefig(pdf, dpi=1200, bbox_inches="tight")
    fig.savefig(png, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(f"[saved] {pdf}")
    print(f"[saved] {png}")


def plot_diagnostics(rows, out_dir, beta=25, nb=2):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    met = branch(rows, "metal")
    ins = branch(rows, "insulator")

    styles = {
        "metal": dict(color="#0072B2", marker="o", label=r"metal seed"),
        "insulator": dict(color="#D55E00", marker="s", label=r"insulator seed"),
    }

    fig, (axC, axN) = plt.subplots(
        1,
        2,
        figsize=(11.2, 4.55),
        dpi=300,
        gridspec_kw={"wspace": 0.28},
    )

    for seed, rr in [("metal", met), ("insulator", ins)]:
        if not rr:
            continue

        U = np.array([r["U"] for r in rr])
        delta = np.array([r["last_delta"] for r in rr])
        delta = np.where(delta > 0, delta, np.nan)

        ndev = np.array([r["navg"] - 1.0 for r in rr])
        n1dev = np.array([r["n1"] - 1.0 for r in rr])
        n2dev = np.array([r["n2"] - 1.0 for r in rr])

        st = styles[seed]

        axC.plot(
            U, delta,
            lw=2.35,
            color=st["color"],
            marker=st["marker"],
            ms=7.0,
            mec=st["color"],
            mfc=st["color"],
            label=st["label"],
        )

        axN.plot(
            U, ndev,
            lw=2.35,
            color=st["color"],
            marker=st["marker"],
            ms=7.0,
            mec=st["color"],
            mfc=st["color"],
            label=st["label"],
            zorder=4,
        )

        axN.plot(
            U, n1dev,
            lw=1.25,
            color=st["color"],
            ls=":",
            alpha=0.45,
            zorder=2,
        )
        axN.plot(
            U, n2dev,
            lw=1.25,
            color=st["color"],
            ls="--",
            alpha=0.45,
            zorder=2,
        )

    axC.set_yscale("log")

    add_panel_label(axC, r"(a)")
    add_panel_label(axN, r"(b)")

    axC.set_title(r"Convergence quality", pad=8)
    axN.set_title(r"Half-filling check", pad=8)

    axC.set_xlabel(r"Interaction strength $(U/t)$")
    axN.set_xlabel(r"Interaction strength $(U/t)$")

    axC.set_ylabel(r"Final RMS $\Delta\Sigma$")
    axN.set_ylabel(r"Half-filling deviation $\bar n-1$")

    axN.axhline(0.0, color="0.45", lw=1.35, zorder=1)

    allU = np.array([r["U"] for r in rows])
    xmin = np.min(allU) - 0.15
    xmax = np.max(allU) + 0.15

    for ax in (axC, axN):
        ax.set_xlim(xmin, xmax)
        ax.xaxis.set_major_locator(MultipleLocator(1.0))
        ax.xaxis.set_minor_locator(MultipleLocator(0.25))
        ax.xaxis.set_major_formatter(FormatStrFormatter("%.1f"))
        common_axis_style(ax)

    axC.yaxis.set_major_locator(LogLocator(base=10.0, numticks=8))
    axC.yaxis.set_minor_formatter(NullFormatter())

    # Symmetric filling axis
    all_dev = np.array([r["navg"] - 1.0 for r in rows])
    lim = max(1e-4, np.nanmax(np.abs(all_dev)) * 1.25)
    lim = min(max(lim, 0.01), 0.12)
    axN.set_ylim(-lim, lim)

    axN.yaxis.set_major_locator(MultipleLocator(0.02))
    axN.yaxis.set_minor_locator(MultipleLocator(0.01))
    axN.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    axC.legend(
        frameon=False,
        loc="upper right",
        handlelength=2.35,
        borderpad=0.2,
        labelspacing=0.38,
    )

    handles2 = [
        Line2D([0], [0], color="0.25", lw=2.25, ls="-", label=r"$\bar n-1$"),
        Line2D([0], [0], color="0.25", lw=1.25, ls=":", label=r"$n_1-1$"),
        Line2D([0], [0], color="0.25", lw=1.25, ls="--", label=r"$n_2-1$"),
    ]

    axN.legend(
        handles=handles2,
        frameon=False,
        loc="upper right",
        handlelength=2.35,
        borderpad=0.2,
        labelspacing=0.32,
    )

    fig.suptitle(
        rf"$N_b={nb}$ ED-DMFT diagnostics, $\beta={beta:.0f}$",
        y=1.02,
        fontsize=20,
    )

    fig.tight_layout()

    pdf = out_dir / "fig_scanA_diagnostics_prb_clean.pdf"
    png = out_dir / "fig_scanA_diagnostics_prb_clean.png"

    fig.savefig(pdf, dpi=1200, bbox_inches="tight")
    fig.savefig(png, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print(f"[saved] {pdf}")
    print(f"[saved] {png}")


def write_summary(rows, out_dir):
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    csv = out_dir / "scanA_hysteresis_summary_clean.csv"

    with open(csv, "w") as f:
        f.write("U,seed,Z1,Z2,Zavg,D1,D2,Davg,n1,n2,navg,converged,last_delta,n_iter,file\n")
        for r in sorted(rows, key=lambda x: (x["seed"], x["U"])):
            f.write(
                f"{r['U']},{r['seed']},{r['Z1']},{r['Z2']},{r['Zavg']},"
                f"{r['D1']},{r['D2']},{r['Davg']},"
                f"{r['n1']},{r['n2']},{r['navg']},"
                f"{r['converged']},{r['last_delta']},{r['n_iter']},{r['file']}\n"
            )

    print(f"[saved] {csv}")

    flagged = [
        r for r in rows
        if (not r["converged"])
        or (np.isfinite(r["last_delta"]) and r["last_delta"] > 2e-2)
        or abs(r["navg"] - 1.0) > 5e-2
    ]

    print("=" * 90)
    print(f"Loaded records : {len(rows)}")
    print(f"Flagged points : {len(flagged)}")
    print("=" * 90)

    if flagged:
        print("U      seed       conv   delta       Zavg     Davg     navg")
        print("-" * 90)
        for r in flagged:
            print(
                f"{r['U']:5.2f}  {r['seed']:9s}  {str(r['converged']):5s}  "
                f"{r['last_delta']:9.2e}  {r['Zavg']:7.4f}  "
                f"{r['Davg']:7.4f}  {r['navg']:7.4f}"
            )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/nb2_scanA_hysteresis")
    ap.add_argument("--out", default="figures_nb2_scanA_hysteresis_prb_clean")
    ap.add_argument("--nb", type=int, default=2)
    ap.add_argument("--beta", type=float, default=25.0)
    ap.add_argument("--t2", type=float, default=1.0)
    ap.add_argument("--z_thresh", type=float, default=0.05)
    ap.add_argument(
        "--include_U0",
        action="store_true",
        help="Include the trivial U=0 point. Default is to exclude it.",
    )
    ap.add_argument("--no_tex", action="store_true")
    args = ap.parse_args()

    set_prb_style(use_tex=not args.no_tex)

    rows = load_records(
        args.data,
        beta=args.beta,
        nb=args.nb,
        t2=args.t2,
        include_U0=args.include_U0,
    )

    write_summary(rows, args.out)
    plot_hysteresis(rows, args.out, z_thresh=args.z_thresh, beta=args.beta, nb=args.nb)
    plot_diagnostics(rows, args.out, beta=args.beta, nb=args.nb)


if __name__ == "__main__":
    main()
