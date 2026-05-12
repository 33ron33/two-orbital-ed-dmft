#!/usr/bin/env python3
"""
PRB-style Figure 1 validation plot for Nb=2 ED-DMFT.

Purpose:
    Validate U=J=0 limit.

Left panel:
    -Im G_m(iw_n) from ED compared with analytic Bethe Green's function.

Right panel:
    |Sigma_m(iw_n)| showing that self-energy vanishes at U=J=0.

This script searches recursively, so it can find U=0 files even if they were
archived inside data/nb2_scanA_hysteresis/_archived_U0_files.
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


# ============================================================
# PRB STYLE
# ============================================================
def set_prb_style(use_tex=True):
    mpl.rcParams.update({
        "text.usetex": use_tex,
        "font.family": "serif",
        "axes.labelsize": 21,
        "font.size": 19,
        "legend.fontsize": 14.5,
        "xtick.labelsize": 17.5,
        "ytick.labelsize": 17.5,
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


def get_key(d, names):
    for name in names:
        if name in d.files:
            return np.asarray(d[name])
    raise KeyError(f"None of these keys found: {names}")


def get_optional_scalar(d, names, default=None):
    for name in names:
        if name in d.files:
            return scalar(d[name])
    return default


def bethe_green_iw(omega, t):
    """
    Noninteracting Bethe lattice local Green's function.

    Self-consistency at U=0:
        G(iw) = 1 / (iw - t^2 G(iw))

    Stable analytic branch:
        G(z) = [z - sqrt(z^2 - 4 t^2)] / (2 t^2), z = i omega
    """
    z = 1j * np.asarray(omega, dtype=float)
    return (z - np.sqrt(z * z - 4.0 * t * t)) / (2.0 * t * t)


def matsubara_grid(beta, n):
    return (2 * np.arange(n) + 1) * np.pi / beta


def find_validation_file(data_dir, beta_target, nb_target, t2_target=None, recursive=True):
    pattern = "**/*.npz" if recursive else "*.npz"
    files = sorted(glob.glob(str(Path(data_dir) / pattern), recursive=recursive))

    candidates = []

    for f in files:
        try:
            d = np.load(f, allow_pickle=True)
        except Exception:
            continue

        if "U" not in d.files or "J" not in d.files:
            continue

        U = scalar(d["U"])
        J = scalar(d["J"])

        beta = get_optional_scalar(d, ["beta"], default=beta_target)
        Nb = get_optional_scalar(d, ["N_bath", "Nb", "N_b"], default=nb_target)
        t2 = get_optional_scalar(d, ["t2"], default=t2_target)

        if abs(U) > 1e-10:
            continue
        if abs(J) > 1e-10:
            continue
        if beta is not None and abs(beta - beta_target) > 1e-8:
            continue
        if Nb is not None and int(round(Nb)) != nb_target:
            continue
        if t2_target is not None and t2 is not None and abs(t2 - t2_target) > 1e-8:
            continue

        # Need G and Sigma arrays
        has_g = ("G1_iw" in d.files and "G2_iw" in d.files)
        has_s = (
            ("S1_iw" in d.files and "S2_iw" in d.files)
            or ("Sigma1_iw" in d.files and "Sigma2_iw" in d.files)
            or ("Sig1_iw" in d.files and "Sig2_iw" in d.files)
        )

        if has_g and has_s:
            candidates.append(f)

    if not candidates:
        raise RuntimeError(
            f"No U=J=0 validation file found in {data_dir}. "
            f"Try pointing --data to your validation folder or archived U0 folder."
        )

    # Prefer newest candidate
    candidates = sorted(candidates, key=lambda x: Path(x).stat().st_mtime)
    return Path(candidates[-1])


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
        0.035,
        0.94,
        label,
        transform=ax.transAxes,
        ha="left",
        va="top",
        fontsize=17,
        fontweight="bold",
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--data", default="data/nb2_scanA_hysteresis")
    ap.add_argument("--out", default="figures_nb2_fig1_validation_prb")
    ap.add_argument("--nb", type=int, default=2)
    ap.add_argument("--beta", type=float, default=25.0)
    ap.add_argument("--t1", type=float, default=1.0)
    ap.add_argument("--t2", type=float, default=1.0)
    ap.add_argument("--wmax", type=float, default=20.0)
    ap.add_argument("--no_tex", action="store_true")
    args = ap.parse_args()

    set_prb_style(use_tex=not args.no_tex)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    f = find_validation_file(
        args.data,
        beta_target=args.beta,
        nb_target=args.nb,
        t2_target=args.t2,
        recursive=True,
    )

    d = np.load(f, allow_pickle=True)

    G1 = get_key(d, ["G1_iw"]).astype(complex)
    G2 = get_key(d, ["G2_iw"]).astype(complex)

    S1 = get_key(d, ["S1_iw", "Sigma1_iw", "Sig1_iw"]).astype(complex)
    S2 = get_key(d, ["S2_iw", "Sigma2_iw", "Sig2_iw"]).astype(complex)

    n = min(len(G1), len(G2), len(S1), len(S2))
    G1 = G1[:n]
    G2 = G2[:n]
    S1 = S1[:n]
    S2 = S2[:n]

    omega = matsubara_grid(args.beta, n)
    mask = omega <= args.wmax

    omega = omega[mask]
    G1 = G1[mask]
    G2 = G2[mask]
    S1 = S1[mask]
    S2 = S2[mask]

    Gbethe1 = bethe_green_iw(omega, args.t1)
    Gbethe2 = bethe_green_iw(omega, args.t2)

    fig, (axG, axS) = plt.subplots(
        1,
        2,
        figsize=(11.3, 4.55),
        dpi=300,
        gridspec_kw={"wspace": 0.27},
    )

    col1 = "#0072B2"
    col2 = "#009E73"

    # ========================================================
    # Left: Green's function validation
    # ========================================================
    axG.plot(
        omega,
        -np.imag(Gbethe1),
        color=col1,
        lw=2.4,
        alpha=0.55,
        label=r"Bethe orb. 1",
    )
    axG.plot(
        omega,
        -np.imag(Gbethe2),
        color=col2,
        lw=2.4,
        ls="--",
        alpha=0.55,
        label=r"Bethe orb. 2",
    )

    axG.plot(
        omega,
        -np.imag(G1),
        linestyle="none",
        marker="o",
        ms=5.8,
        color=col1,
        label=r"ED orb. 1",
    )
    axG.plot(
        omega,
        -np.imag(G2),
        linestyle="none",
        marker="s",
        ms=5.8,
        color=col2,
        label=r"ED orb. 2",
    )

    axG.set_title(r"$U=J=0$ Green's function", pad=8)
    axG.set_xlabel(r"Matsubara frequency $\omega_n$")
    axG.set_ylabel(r"$-\mathrm{Im}\,G_m(i\omega_n)$")

    axG.set_xlim(0.0, args.wmax)
    ymax = max(np.max(-np.imag(G1)), np.max(-np.imag(G2)), np.max(-np.imag(Gbethe1)), np.max(-np.imag(Gbethe2)))
    axG.set_ylim(0.0, 1.08 * ymax)

    axG.xaxis.set_major_locator(MultipleLocator(5.0))
    axG.xaxis.set_minor_locator(MultipleLocator(1.0))
    axG.yaxis.set_major_locator(MultipleLocator(0.25))
    axG.yaxis.set_minor_locator(MultipleLocator(0.125))
    axG.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))
    axG.yaxis.set_major_formatter(FormatStrFormatter("%.2f"))

    axG.legend(
        frameon=False,
        loc="upper right",
        handlelength=2.4,
        borderpad=0.2,
        labelspacing=0.35,
    )

    # ========================================================
    # Right: self-energy vanishes
    # ========================================================
    eps = 1e-18
    absS1 = np.maximum(np.abs(S1), eps)
    absS2 = np.maximum(np.abs(S2), eps)

    axS.plot(
        omega,
        absS1,
        color=col1,
        marker="o",
        ms=5.4,
        lw=2.2,
        label=r"$|\Sigma_1|$",
    )
    axS.plot(
        omega,
        absS2,
        color=col2,
        marker="s",
        ms=5.4,
        lw=2.2,
        label=r"$|\Sigma_2|$",
    )

    axS.set_yscale("log")
    axS.set_title(r"Self-energy vanishes", pad=8)
    axS.set_xlabel(r"Matsubara frequency $\omega_n$")
    axS.set_ylabel(r"$|\Sigma_m(i\omega_n)|$")

    axS.set_xlim(0.0, args.wmax)

    finite = np.concatenate([absS1[np.isfinite(absS1)], absS2[np.isfinite(absS2)]])
    ymin = max(1e-16, 0.5 * np.nanmin(finite))
    ymax_s = min(1e-9, 2.0 * np.nanmax(finite))
    if ymax_s <= ymin:
        ymax_s = ymin * 100.0

    axS.set_ylim(ymin, ymax_s)

    axS.xaxis.set_major_locator(MultipleLocator(5.0))
    axS.xaxis.set_minor_locator(MultipleLocator(1.0))
    axS.xaxis.set_major_formatter(FormatStrFormatter("%.0f"))

    axS.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
    axS.yaxis.set_minor_formatter(NullFormatter())

    axS.legend(
        frameon=False,
        loc="upper left",
        handlelength=2.4,
        borderpad=0.2,
        labelspacing=0.35,
    )

    add_panel_label(axG, r"(a)")
    add_panel_label(axS, r"(b)")

    for ax in (axG, axS):
        common_axis_style(ax)

    fig.suptitle(
        rf"$N_b={args.nb}$ ED-DMFT validation, $\beta={args.beta:.0f}$",
        y=1.02,
        fontsize=20,
    )

    fig.tight_layout()

    pdf = out_dir / "fig1_validation_green_sigma_prb.pdf"
    png = out_dir / "fig1_validation_green_sigma_prb.png"

    fig.savefig(pdf, dpi=1200, bbox_inches="tight")
    fig.savefig(png, dpi=400, bbox_inches="tight")
    plt.close(fig)

    print("=" * 90)
    print(f"Validation file used: {f}")
    print(f"[saved] {pdf}")
    print(f"[saved] {png}")
    print("=" * 90)


if __name__ == "__main__":
    main()
