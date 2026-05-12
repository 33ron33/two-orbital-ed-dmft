#!/usr/bin/env python3
"""
plot_ed_dmft_all.py

Comprehensive plotting/diagnostic script for two-orbital Kanamori ED-DMFT NPZ files.

Works for two use cases:

  1) Smoke test / single-point diagnostics
     - Loads one or more ED_*.npz files from --data/--smoke_dir.
     - Makes quick dashboard plots: G(iw), Sigma(iw), observables, bath parameters.
     - Optionally parses a log file and plots convergence history.

  2) Production paper figures
     - Figure 1: U=J=0 Green's-function validation against analytic Bethe result.
     - Figure 2: Scan A, J=0 Mott transition: Z and D vs U, metal/insulator branches.
     - Figure 3: Scan C, Z1/Z2 vs U for several J/U values.
     - Figure 4: Self-energy inside OSMT, default U=3, J/U=0.10.
     - Figure 5: Scan B, D1/D2/D12 vs J/U at fixed U.
     - Figure 6: Scan C, (U,J/U) phase diagram.
     - Figure 7: Scan D, (U,T) phase diagram plus Z(U) at lowest T.

The script is intentionally robust: if a figure's data are missing, it skips that
figure instead of crashing. It also writes all loaded scalar observables to CSV.

Example smoke test:

  python plot_ed_dmft_all.py \
      --mode smoke \
      --data smoke_ED \
      --logs logs/smoke_ED.log \
      --out figs_smoke \
      --nb 2 --beta 8 --t2 1.0

Example production:

  python plot_ed_dmft_all.py \
      --mode production \
      --scanA /lustre/.../data/scanA \
      --scanB /lustre/.../data/scanB \
      --scanC /lustre/.../data/scanC \
      --scanD /lustre/.../data/scanD \
      --out figures_prod \
      --nb 2 --beta 25 --t2A 1.0 --t2B 1.0 --t2C 0.5 --t2D 0.5 \
      --z_thresh 0.05
"""

from __future__ import annotations

import argparse
import csv
import glob
import math
import os
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, BoundaryNorm
from matplotlib.lines import Line2D


# -----------------------------------------------------------------------------
# Style
# -----------------------------------------------------------------------------

COL_ORB1 = "#378ADD"   # blue  — orbital 1 / wide band
COL_ORB2 = "#1D9E75"   # green — orbital 2 / narrow band
COL_QMC  = "#03A2FE"   # cyan  — reference
COL_INS  = "#E24B4A"   # red   — insulator branch / Mott
COL_OSMT = "#F0A202"   # orange — OSMT
COL_GREY = "0.35"

PHASE_LABELS = {
    0: "both metal",
    1: "OSMT",
    2: "both Mott",
    3: "missing/other",
}


def setup_matplotlib(usetex: bool = False) -> None:
    """PRB-like matplotlib defaults. Disable usetex if LaTeX is unavailable."""
    plt.rcParams.update({
        "text.usetex"        : bool(usetex),
        "font.family"        : "serif",
        "font.serif"         : ["Computer Modern Roman", "DejaVu Serif"],
        "mathtext.fontset"   : "cm",
        "axes.labelsize"     : 10,
        "axes.titlesize"     : 10,
        "xtick.labelsize"    : 9,
        "ytick.labelsize"    : 9,
        "legend.fontsize"    : 8,
        "xtick.direction"    : "in",
        "ytick.direction"    : "in",
        "xtick.top"          : True,
        "ytick.right"        : True,
        "xtick.minor.visible": True,
        "ytick.minor.visible": True,
        "axes.linewidth"     : 0.8,
        "lines.linewidth"    : 1.5,
        "savefig.bbox"       : "tight",
        "pdf.fonttype"       : 42,
        "ps.fonttype"        : 42,
    })


def savefig(fig: plt.Figure, outdir: Path, name: str, dpi: int = 600, png: bool = True) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    pdf = outdir / f"{name}.pdf"
    fig.savefig(pdf, dpi=dpi)
    if png:
        fig.savefig(outdir / f"{name}.png", dpi=min(dpi, 300))
    print(f"[saved] {pdf}")
    plt.close(fig)


# -----------------------------------------------------------------------------
# Data model / loading
# -----------------------------------------------------------------------------

@dataclass
class EDRecord:
    path: Path
    branch: str
    U: float
    J: float
    beta: float
    t1: float
    t2: float
    nb: int
    converged: bool
    n_iter: int
    last_delta: float
    Z1: float
    Z2: float
    D1: float
    D2: float
    D12: float
    n1: float
    n2: float
    iw_vals: np.ndarray
    G1: np.ndarray
    G2: np.ndarray
    S1: np.ndarray
    S2: np.ndarray
    eps_bath: Optional[np.ndarray]
    V_bath: Optional[np.ndarray]

    @property
    def Jr(self) -> float:
        return self.J / self.U if abs(self.U) > 1e-14 else np.nan

    @property
    def T(self) -> float:
        return 1.0 / self.beta if self.beta != 0 else np.nan

    @property
    def Zavg(self) -> float:
        return 0.5 * (self.Z1 + self.Z2)

    @property
    def Davg(self) -> float:
        return 0.5 * (self.D1 + self.D2)


def scalar_from_npz(d: np.lib.npyio.NpzFile, key: str, default=np.nan):
    if key not in d:
        return default
    arr = d[key]
    try:
        return arr.item()
    except Exception:
        return np.asarray(arr).reshape(-1)[0]


def parse_filename(path: Path) -> Dict[str, float | str | int]:
    """Parse ED_MET_U5.00_J0.750_t20.50_b25.0_Nb2.npz style names."""
    info: Dict[str, float | str | int] = {}
    m = re.search(
        r"ED_(?P<branch>MET|INS)_U(?P<U>[-+0-9.]+)_J(?P<J>[-+0-9.]+)_t2(?P<t2>[-+0-9.]+)_b(?P<beta>[-+0-9.]+)_Nb(?P<nb>\d+)",
        path.name,
    )
    if m:
        info["branch"] = m.group("branch")
        info["U"] = float(m.group("U"))
        info["J"] = float(m.group("J"))
        info["t2"] = float(m.group("t2"))
        info["beta"] = float(m.group("beta"))
        info["nb"] = int(m.group("nb"))
    return info


def reconstruct_iw_vals(beta: float, n_total: int) -> np.ndarray:
    """
    Reconstruct symmetric fermionic Matsubara mesh if iw_vals is absent.

    TRIQS MeshImFreq with n_iw positive frequencies usually has 2*n_iw points
    ordered from negative to positive. This reconstruction is sufficient for plots.
    """
    half = n_total // 2
    n = np.arange(-half, n_total - half)
    # Fermionic frequencies: omega_n = (2n+1) pi / beta.
    omega = (2 * n + 1) * np.pi / beta
    return 1j * omega


def load_record(path: Path) -> Optional[EDRecord]:
    try:
        d = np.load(path, allow_pickle=False)
    except Exception as exc:
        print(f"[skip] failed to load {path}: {exc}")
        return None

    fn = parse_filename(path)

    def get_float(key: str, default=np.nan):
        val = scalar_from_npz(d, key, fn.get(key, default))
        try:
            return float(val)
        except Exception:
            return float(default)

    def get_int(key: str, default=0):
        val = scalar_from_npz(d, key, fn.get(key, default))
        try:
            return int(val)
        except Exception:
            return int(default)

    branch = str(scalar_from_npz(d, "seed", fn.get("branch", "MET")))
    if branch.lower().startswith("metal"):
        branch = "MET"
    elif branch.lower().startswith("ins"):
        branch = "INS"
    elif branch not in {"MET", "INS"}:
        branch = str(fn.get("branch", "MET"))

    if "G1_iw" not in d or "G2_iw" not in d or "S1_iw" not in d or "S2_iw" not in d:
        print(f"[skip] {path}: missing one of G1_iw/G2_iw/S1_iw/S2_iw")
        return None

    G1 = np.asarray(d["G1_iw"], dtype=np.complex128)
    G2 = np.asarray(d["G2_iw"], dtype=np.complex128)
    S1 = np.asarray(d["S1_iw"], dtype=np.complex128)
    S2 = np.asarray(d["S2_iw"], dtype=np.complex128)

    beta = get_float("beta", float(fn.get("beta", np.nan)))
    if "iw_vals" in d:
        iw_vals = np.asarray(d["iw_vals"], dtype=np.complex128)
    else:
        iw_vals = reconstruct_iw_vals(beta, len(G1))

    eps_bath = np.asarray(d["eps_bath"]) if "eps_bath" in d else None
    V_bath = np.asarray(d["V_bath"]) if "V_bath" in d else None

    return EDRecord(
        path=path,
        branch=branch,
        U=get_float("U", float(fn.get("U", np.nan))),
        J=get_float("J", float(fn.get("J", np.nan))),
        beta=beta,
        t1=get_float("t1", 1.0),
        t2=get_float("t2", float(fn.get("t2", np.nan))),
        nb=get_int("N_bath", int(fn.get("nb", 0))),
        converged=bool(scalar_from_npz(d, "converged", False)),
        n_iter=get_int("n_iter", 0),
        last_delta=get_float("last_delta", np.nan),
        Z1=get_float("Z1", np.nan),
        Z2=get_float("Z2", np.nan),
        D1=get_float("D1", np.nan),
        D2=get_float("D2", np.nan),
        D12=get_float("D12", np.nan),
        n1=get_float("n1", np.nan),
        n2=get_float("n2", np.nan),
        iw_vals=iw_vals,
        G1=G1,
        G2=G2,
        S1=S1,
        S2=S2,
        eps_bath=eps_bath,
        V_bath=V_bath,
    )


def expand_dirs(paths: Sequence[str]) -> List[Path]:
    files: List[Path] = []
    for p in paths:
        if not p:
            continue
        pp = Path(p).expanduser()
        if pp.is_file() and pp.suffix == ".npz":
            files.append(pp)
        elif pp.is_dir():
            files.extend(sorted(pp.glob("*.npz")))
            files.extend(sorted(pp.glob("**/*.npz")))
        else:
            files.extend(Path(x) for x in glob.glob(str(pp)))
    # Deduplicate while preserving order.
    seen = set()
    out = []
    for f in files:
        r = f.resolve()
        if r not in seen and f.name.endswith(".npz"):
            seen.add(r)
            out.append(f)
    return out


def load_records(paths: Sequence[str]) -> List[EDRecord]:
    files = expand_dirs(paths)
    records = []
    for f in files:
        r = load_record(f)
        if r is not None:
            records.append(r)
    records.sort(key=lambda x: (x.beta, x.t2, x.U, x.J, x.branch, str(x.path)))
    print(f"[load] {len(records)} ED records from {len(files)} npz files")
    return records


def write_summary_csv(records: Sequence[EDRecord], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / "ed_dmft_summary.csv"
    fields = [
        "path", "branch", "U", "J", "J_over_U", "beta", "T", "t1", "t2", "Nb",
        "converged", "n_iter", "last_delta", "Z1", "Z2", "D1", "D2", "D12", "n1", "n2",
    ]
    with path.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in records:
            w.writerow({
                "path": str(r.path), "branch": r.branch, "U": r.U, "J": r.J,
                "J_over_U": r.Jr, "beta": r.beta, "T": r.T, "t1": r.t1, "t2": r.t2,
                "Nb": r.nb, "converged": r.converged, "n_iter": r.n_iter,
                "last_delta": r.last_delta, "Z1": r.Z1, "Z2": r.Z2,
                "D1": r.D1, "D2": r.D2, "D12": r.D12, "n1": r.n1, "n2": r.n2,
            })
    print(f"[saved] {path}")


def approx(a: float, b: float, tol: float = 1e-8) -> bool:
    if np.isnan(a) or np.isnan(b):
        return False
    return abs(a - b) <= tol


def filter_records(
    records: Sequence[EDRecord],
    *,
    branch: Optional[str] = None,
    beta: Optional[float] = None,
    t2: Optional[float] = None,
    nb: Optional[int] = None,
    U: Optional[float] = None,
    J: Optional[float] = None,
    Jr: Optional[float] = None,
    tol: float = 1e-6,
) -> List[EDRecord]:
    out: List[EDRecord] = []
    for r in records:
        if branch is not None and r.branch != branch:
            continue
        if beta is not None and not approx(r.beta, beta, tol):
            continue
        if t2 is not None and not approx(r.t2, t2, tol):
            continue
        if nb is not None and r.nb != nb:
            continue
        if U is not None and not approx(r.U, U, tol):
            continue
        if J is not None and not approx(r.J, J, tol):
            continue
        if Jr is not None and not approx(round(r.Jr, 8), Jr, max(tol, 1e-5)):
            continue
        out.append(r)
    return out


def best_by_grid(records: Sequence[EDRecord]) -> List[EDRecord]:
    """If duplicate grid points exist, prefer converged, then lower last_delta, then latest mtime."""
    buckets: Dict[Tuple, List[EDRecord]] = {}
    for r in records:
        key = (r.branch, round(r.U, 8), round(r.J, 8), round(r.beta, 8), round(r.t2, 8), r.nb)
        buckets.setdefault(key, []).append(r)
    out = []
    for group in buckets.values():
        group = sorted(
            group,
            key=lambda r: (
                not r.converged,
                np.inf if np.isnan(r.last_delta) else r.last_delta,
                -r.path.stat().st_mtime if r.path.exists() else 0.0,
            ),
        )
        out.append(group[0])
    out.sort(key=lambda r: (r.beta, r.t2, r.Jr if not np.isnan(r.Jr) else -1, r.U, r.branch))
    return out


def pos_freq(r: EDRecord) -> Tuple[np.ndarray, slice]:
    n = len(r.iw_vals)
    n0 = n // 2
    sl = slice(n0, n)
    omega = np.imag(r.iw_vals[sl])
    # If reconstruction ordering produced negative/zero weirdness, sort positive by value.
    return omega, sl


# -----------------------------------------------------------------------------
# Log parsing for smoke/production diagnostics
# -----------------------------------------------------------------------------

@dataclass
class IterLog:
    it: np.ndarray
    dSigma: np.ndarray
    D1: np.ndarray
    D2: np.ndarray
    Z1: np.ndarray
    Z2: np.ndarray
    n1: np.ndarray
    n2: np.ndarray
    fit1: np.ndarray
    fit2: np.ndarray


def parse_log_file(path: str) -> Optional[IterLog]:
    """Parse solver/ed_solver_dmft.py logs robustly.

    The previous single-regex parser could miss dSigma_rms because the field was
    optional. Here each quantity is extracted independently from the same line.
    """
    p = Path(path).expanduser()
    if not p.exists():
        print(f"[skip] log file not found: {p}")
        return None

    num = r"[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?"

    def find_float(line: str, key: str, default=np.nan) -> float:
        m = re.search(rf"{re.escape(key)}=\s*({num})", line)
        return float(m.group(1)) if m else float(default)

    it, dS, D1, D2, Z1, Z2, n1, n2, fit1, fit2 = [], [], [], [], [], [], [], [], [], []

    for line in p.read_text(errors="ignore").splitlines():
        if "it=" not in line or "D1=" not in line or "Z1=" not in line:
            continue
        m_it = re.search(r"it=\s*(\d+)", line)
        if not m_it:
            continue

        it.append(int(m_it.group(1)))
        dS.append(find_float(line, "dSigma_rms", np.nan))
        D1.append(find_float(line, "D1"))
        D2.append(find_float(line, "D2"))
        Z1.append(find_float(line, "Z1"))
        Z2.append(find_float(line, "Z2"))
        n1.append(find_float(line, "n1"))
        n2.append(find_float(line, "n2"))

        mf = re.search(rf"fit=\(\s*({num})\s*,\s*({num})\s*\)", line)
        fit1.append(float(mf.group(1)) if mf else np.nan)
        fit2.append(float(mf.group(2)) if mf else np.nan)

    if not it:
        print(f"[skip] no iteration lines parsed from {p}")
        return None

    return IterLog(
        it=np.array(it), dSigma=np.array(dS), D1=np.array(D1), D2=np.array(D2),
        Z1=np.array(Z1), Z2=np.array(Z2), n1=np.array(n1), n2=np.array(n2),
        fit1=np.array(fit1), fit2=np.array(fit2),
    )


def plot_log_convergence(log: IterLog, outdir: Path, name: str = "smoke_log_convergence") -> None:
    fig, axes = plt.subplots(2, 2, figsize=(6.75, 4.8))
    ax = axes[0, 0]
    mask = np.isfinite(log.dSigma)
    if mask.any():
        ax.semilogy(log.it[mask], log.dSigma[mask], "o-", color=COL_GREY)
    ax.set_xlabel("DMFT iteration")
    ax.set_ylabel(r"RMS $\Delta\Sigma$")
    ax.set_title("Convergence")

    ax = axes[0, 1]
    ax.plot(log.it, log.Z1, "o-", color=COL_ORB1, label=r"$Z_1$")
    ax.plot(log.it, log.Z2, "s-", color=COL_ORB2, label=r"$Z_2$")
    ax.set_xlabel("DMFT iteration")
    ax.set_ylabel(r"$Z_m$")
    ax.legend(frameon=False)

    ax = axes[1, 0]
    ax.plot(log.it, log.D1, "o-", color=COL_ORB1, label=r"$D_1$")
    ax.plot(log.it, log.D2, "s-", color=COL_ORB2, label=r"$D_2$")
    ax.set_xlabel("DMFT iteration")
    ax.set_ylabel(r"$D_m$")
    ax.legend(frameon=False)

    ax = axes[1, 1]
    if np.isfinite(log.fit1).any() or np.isfinite(log.fit2).any():
        ax.semilogy(log.it, log.fit1, "o-", color=COL_ORB1, label="orb. 1")
        ax.semilogy(log.it, log.fit2, "s-", color=COL_ORB2, label="orb. 2")
        ax.set_ylabel("bath-fit residual")
    else:
        ax.plot(log.it, log.n1, "o-", color=COL_ORB1, label=r"$n_1$")
        ax.plot(log.it, log.n2, "s-", color=COL_ORB2, label=r"$n_2$")
        ax.set_ylabel(r"$n_m$")
    ax.set_xlabel("DMFT iteration")
    ax.legend(frameon=False)

    fig.tight_layout()
    savefig(fig, outdir, name)


# -----------------------------------------------------------------------------
# Analytic non-interacting Bethe Green's function
# -----------------------------------------------------------------------------


def bethe_G0_analytic(z: np.ndarray, t: float) -> np.ndarray:
    """Local Green's function for semicircular DOS with Bethe hopping t."""
    out = np.empty_like(z, dtype=np.complex128)
    for i, zz in enumerate(z):
        if abs(t) < 1e-14:
            out[i] = 1.0 / zz
            continue
        disc = zz**2 - 4.0 * t**2
        root = np.sqrt(disc)
        g1 = (zz - root) / (2.0 * t**2)
        g2 = (zz + root) / (2.0 * t**2)
        # Choose root with correct high-frequency asymptotic 1/z.
        out[i] = g1 if abs(g1 - 1.0 / zz) < abs(g2 - 1.0 / zz) else g2
    return out


# -----------------------------------------------------------------------------
# Smoke/single-run plots
# -----------------------------------------------------------------------------


def plot_smoke_dashboard(records: Sequence[EDRecord], outdir: Path, nfreq: int = 60) -> None:
    if not records:
        return
    # Prefer most recently modified file.
    r = sorted(records, key=lambda x: x.path.stat().st_mtime if x.path.exists() else 0.0)[-1]
    omega, sl_all = pos_freq(r)
    nshow = min(nfreq, len(omega))
    sl = slice(sl_all.start, sl_all.start + nshow)
    w = omega[:nshow]

    # 2x3 layout avoids the confusing twin-axis scientific offset in the bath plot.
    fig, axes = plt.subplots(2, 3, figsize=(8.2, 5.0))

    ax = axes[0, 0]
    ax.plot(w, -r.G1[sl].imag, "o-", ms=3, color=COL_ORB1, label=r"$-\mathrm{Im}\,G_1$")
    ax.plot(w, -r.G2[sl].imag, "s-", ms=3, color=COL_ORB2, label=r"$-\mathrm{Im}\,G_2$")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$-\mathrm{Im}\,G_m(i\omega_n)$")
    ax.set_title("Impurity Green's function")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[0, 1]
    ax.plot(w, -r.S1[sl].imag, "o-", ms=3, color=COL_ORB1, label=r"$-\mathrm{Im}\,\Sigma_1$")
    ax.plot(w, -r.S2[sl].imag, "s-", ms=3, color=COL_ORB2, label=r"$-\mathrm{Im}\,\Sigma_2$")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$-\mathrm{Im}\,\Sigma_m(i\omega_n)$")
    ax.set_title("Self-energy")
    ax.legend(frameon=False, fontsize=7)

    ax = axes[0, 2]
    labels = [r"$Z_1$", r"$Z_2$", r"$D_1$", r"$D_2$", r"$n_1$", r"$n_2$"]
    vals = [r.Z1, r.Z2, r.D1, r.D2, r.n1, r.n2]
    ax.bar(np.arange(len(vals)), vals, color=[COL_ORB1, COL_ORB2, COL_ORB1, COL_ORB2, COL_ORB1, COL_ORB2])
    ax.set_xticks(np.arange(len(vals)), labels)
    ax.set_ylabel("value")
    ax.set_title("Observables")

    ax = axes[1, 0]
    txt = (
        fr"$U={r.U:.2f}$, $J={r.J:.3f}$" + "\n" +
        fr"$\beta={r.beta:g}$, $t_2={r.t2:g}$, $N_b={r.nb}$" + "\n" +
        fr"converged = {r.converged}" + "\n" +
        fr"iterations = {r.n_iter}" + "\n" +
        fr"last $\Delta\Sigma$ = {r.last_delta:.2e}" + "\n" +
        fr"$n_1-1={r.n1-1:.1e}$, $n_2-1={r.n2-1:.1e}$"
    )
    ax.axis("off")
    ax.text(0.02, 0.98, txt, transform=ax.transAxes, va="top", ha="left", fontsize=9)

    ax = axes[1, 1]
    if r.eps_bath is not None:
        k = np.arange(r.eps_bath.shape[1])
        ax.plot(k, r.eps_bath[0], "o-", color=COL_ORB1, label=r"$\epsilon_{1k}$")
        ax.plot(k, r.eps_bath[1], "s-", color=COL_ORB2, label=r"$\epsilon_{2k}$")
        ax.axhline(0, color="0.5", lw=0.8)
        ax.set_xlabel("bath index $k$")
        ax.set_ylabel(r"$\epsilon_{mk}$")
        ax.set_title("Bath energies")
        ax.legend(frameon=False, fontsize=7)
    else:
        ax.axis("off")

    ax = axes[1, 2]
    if r.V_bath is not None:
        k = np.arange(r.V_bath.shape[1])
        ax.plot(k, r.V_bath[0], "o-", color=COL_ORB1, label=r"$V_{1k}$")
        ax.plot(k, r.V_bath[1], "s-", color=COL_ORB2, label=r"$V_{2k}$")
        ax.set_xlabel("bath index $k$")
        ax.set_ylabel(r"$V_{mk}$")
        ax.set_title("Bath hybridizations")
        ax.legend(frameon=False, fontsize=7)
    else:
        ax.axis("off")

    fig.suptitle(f"Smoke dashboard: {r.path.name}", fontsize=10)
    fig.tight_layout()
    savefig(fig, outdir, "smoke_dashboard")


def fig1_green_validation(records: Sequence[EDRecord], outdir: Path, beta: Optional[float], t2: Optional[float], nb: Optional[int]) -> None:
    cand = [r for r in records if approx(r.U, 0.0, 1e-8) and approx(r.J, 0.0, 1e-8)]
    if beta is not None:
        cand = [r for r in cand if approx(r.beta, beta, 1e-6)]
    if t2 is not None:
        cand = [r for r in cand if approx(r.t2, t2, 1e-6)]
    if nb is not None:
        cand = [r for r in cand if r.nb == nb]
    cand = best_by_grid(cand)
    if not cand:
        print("[skip] Figure 1: no U=J=0 file found")
        return
    r = cand[0]
    omega, sl_all = pos_freq(r)
    nshow = min(80, len(omega))
    sl = slice(sl_all.start, sl_all.start + nshow)
    z = r.iw_vals[sl]
    w = omega[:nshow]
    G1_exact = bethe_G0_analytic(z, r.t1)
    G2_exact = bethe_G0_analytic(z, r.t2)

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))
    ax = axes[0]
    ax.plot(w, -r.G1[sl].imag, "o", ms=3, color=COL_ORB1, label="ED orb. 1")
    ax.plot(w, -G1_exact.imag, "-", color=COL_ORB1, alpha=0.65, label="Bethe orb. 1")
    ax.plot(w, -r.G2[sl].imag, "s", ms=3, color=COL_ORB2, label="ED orb. 2")
    ax.plot(w, -G2_exact.imag, "--", color=COL_ORB2, alpha=0.65, label="Bethe orb. 2")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$-\mathrm{Im}\,G(i\omega_n)$")
    ax.set_title(r"$U=J=0$ Green's function")
    ax.legend(frameon=False, ncol=1, fontsize=7)

    ax = axes[1]
    absS1 = np.abs(r.S1[sl])
    absS2 = np.abs(r.S2[sl])
    # Avoid zeros on a log axis.
    floor = 1e-18
    ax.semilogy(w, np.maximum(absS1, floor), "o-", ms=3, color=COL_ORB1, label=r"$|\Sigma_1|$")
    ax.semilogy(w, np.maximum(absS2, floor), "s-", ms=3, color=COL_ORB2, label=r"$|\Sigma_2|$")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$|\Sigma(i\omega_n)|$")
    ax.set_title(r"Self-energy vanishes")
    ax.legend(frameon=False, fontsize=7)
    fig.tight_layout()
    savefig(fig, outdir, "fig1_green_validation")


def fig2_mott_scan_A(records: Sequence[EDRecord], outdir: Path, beta: float, t2: float, nb: Optional[int]) -> None:
    rec = [r for r in records if approx(r.beta, beta, 1e-6) and approx(r.t2, t2, 1e-6) and approx(r.J, 0.0, 1e-6)]
    if nb is not None:
        rec = [r for r in rec if r.nb == nb]
    rec = best_by_grid(rec)
    if len(rec) < 2:
        print("[skip] Figure 2: not enough Scan A J=0 records")
        return

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))
    for branch, marker, color, label in [("MET", "o", COL_ORB1, "metal seed"), ("INS", "s", COL_INS, "insulator seed")]:
        rr = sorted([r for r in rec if r.branch == branch], key=lambda x: x.U)
        if not rr:
            continue
        U = np.array([x.U for x in rr])
        Z = np.array([x.Zavg for x in rr])
        D = np.array([x.Davg for x in rr])
        axes[0].plot(U, Z, marker + "-", color=color, label=label)
        axes[1].plot(U, D, marker + "-", color=color, label=label)
    axes[0].set_xlabel(r"$U/t$")
    axes[0].set_ylabel(r"$Z$")
    axes[0].set_title(r"$J=0$, degenerate")
    axes[0].legend(frameon=False)
    axes[1].set_xlabel(r"$U/t$")
    axes[1].set_ylabel(r"$D=\langle n_\uparrow n_\downarrow\rangle$")
    axes[1].legend(frameon=False)
    fig.tight_layout()
    savefig(fig, outdir, "fig2_scanA_Z_D_vs_U")


def unique_rounded(vals: Iterable[float], ndigits: int = 3) -> List[float]:
    out = sorted(set(round(float(v), ndigits) for v in vals if np.isfinite(v)))
    return out


def fig3_Z_vs_U_scan_C(records: Sequence[EDRecord], outdir: Path, beta: float, t2: float, nb: Optional[int]) -> None:
    rec = [r for r in records if r.branch == "MET" and approx(r.beta, beta, 1e-6) and approx(r.t2, t2, 1e-6)]
    if nb is not None:
        rec = [r for r in rec if r.nb == nb]
    rec = best_by_grid(rec)
    if len(rec) < 3:
        print("[skip] Figure 3: not enough Scan C records")
        return

    jrs = unique_rounded([r.Jr for r in rec], 3)
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8), sharey=True)
    cmap = plt.get_cmap("viridis")
    denom = max(1, len(jrs) - 1)
    for i, jr in enumerate(jrs):
        rr = sorted([r for r in rec if approx(round(r.Jr, 3), jr, 1e-8)], key=lambda x: x.U)
        if len(rr) < 2:
            continue
        color = cmap(i / denom)
        U = np.array([x.U for x in rr])
        axes[0].plot(U, [x.Z1 for x in rr], "o-", color=color, label=fr"$J/U={jr:.2f}$")
        axes[1].plot(U, [x.Z2 for x in rr], "s-", color=color, label=fr"$J/U={jr:.2f}$")
    axes[0].set_xlabel(r"$U/t_1$")
    axes[1].set_xlabel(r"$U/t_1$")
    axes[0].set_ylabel(r"$Z_m$")
    axes[0].set_title(r"wide orbital, $Z_1$")
    axes[1].set_title(r"narrow orbital, $Z_2$")
    axes[1].legend(frameon=False, fontsize=7)
    fig.tight_layout()
    savefig(fig, outdir, "fig3_scanC_Z_vs_U_by_Jratio")


def nearest_record(records: Sequence[EDRecord], U: float, Jr: float, beta: float, t2: float, nb: Optional[int], branch: str = "MET") -> Optional[EDRecord]:
    cand = [r for r in records if r.branch == branch and approx(r.beta, beta, 1e-6) and approx(r.t2, t2, 1e-6)]
    if nb is not None:
        cand = [r for r in cand if r.nb == nb]
    if not cand:
        return None
    return min(cand, key=lambda r: abs(r.U - U) + 5.0 * abs(r.Jr - Jr))


def fig4_self_energy_osmt(records: Sequence[EDRecord], outdir: Path, beta: float, t2: float, nb: Optional[int], U_osmt: float, Jr_osmt: float) -> None:
    r = nearest_record(records, U_osmt, Jr_osmt, beta, t2, nb, branch="MET")
    if r is None:
        print("[skip] Figure 4: no OSMT candidate record")
        return
    omega, sl_all = pos_freq(r)
    nshow = min(80, len(omega))
    sl = slice(sl_all.start, sl_all.start + nshow)
    w = omega[:nshow]

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))
    ax = axes[0]
    ax.plot(w, -r.S1[sl].imag, "o-", ms=3, color=COL_ORB1, label=fr"orb. 1, $Z_1={r.Z1:.3g}$")
    ax.plot(w, -r.S2[sl].imag, "s-", ms=3, color=COL_ORB2, label=fr"orb. 2, $Z_2={r.Z2:.3g}$")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$-\mathrm{Im}\,\Sigma_m(i\omega_n)$")
    ax.set_title(fr"Self-energy, $U={r.U:.2f}$, $J/U={r.Jr:.2f}$")
    ax.legend(frameon=False)

    ax = axes[1]
    ax.plot(w, -r.G1[sl].imag, "o-", ms=3, color=COL_ORB1, label=r"$-\mathrm{Im}\,G_1$")
    ax.plot(w, -r.G2[sl].imag, "s-", ms=3, color=COL_ORB2, label=r"$-\mathrm{Im}\,G_2$")
    ax.set_xlabel(r"$\omega_n$")
    ax.set_ylabel(r"$-\mathrm{Im}\,G_m(i\omega_n)$")
    ax.legend(frameon=False)
    fig.tight_layout()
    savefig(fig, outdir, "fig4_self_energy_osmt")


def fig5_hund_docc_scan_B(records: Sequence[EDRecord], outdir: Path, beta: float, t2: float, nb: Optional[int], U_fixed: float) -> None:
    rec = [r for r in records if r.branch == "MET" and approx(r.beta, beta, 1e-6) and approx(r.t2, t2, 1e-6) and approx(r.U, U_fixed, 1e-6)]
    if nb is not None:
        rec = [r for r in rec if r.nb == nb]
    rec = best_by_grid(rec)
    rec = sorted(rec, key=lambda x: x.Jr)
    if len(rec) < 2:
        print("[skip] Figure 5: not enough Scan B records")
        return

    Jr = np.array([r.Jr for r in rec])
    fig, ax = plt.subplots(figsize=(3.375, 2.8))
    ax.plot(Jr, [r.D1 for r in rec], "o-", color=COL_ORB1, label=r"$D_1$")
    ax.plot(Jr, [r.D2 for r in rec], "s-", color=COL_ORB2, label=r"$D_2$")
    ax.set_xlabel(r"$J/U$")
    ax.set_ylabel(r"intra-orbital $D_m$")

    ax2 = ax.twinx()
    ax2.plot(Jr, [r.D12 for r in rec], "^-", color=COL_INS, label=r"$D_{12}$")
    ax2.set_ylabel(r"inter-orbital $D_{12}$")

    ax.text(
        0.05, 0.08,
        r"Hund drives $D_{12}\to 1$" + "\n" + r"one electron per orbital",
        transform=ax.transAxes,
        fontsize=8,
        va="bottom",
    )
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, loc="best")
    fig.tight_layout()
    savefig(fig, outdir, "fig5_scanB_D_vs_Jratio")


def classify_phase(Z1: float, Z2: float, z_thresh: float) -> int:
    if not np.isfinite(Z1) or not np.isfinite(Z2):
        return 3
    m1 = Z1 > z_thresh
    m2 = Z2 > z_thresh
    if m1 and m2:
        return 0
    if m1 != m2:
        return 1
    return 2


def phase_cmap_norm():
    cmap = ListedColormap([COL_ORB1, COL_OSMT, COL_INS, "0.8"])
    norm = BoundaryNorm([-0.5, 0.5, 1.5, 2.5, 3.5], cmap.N)
    return cmap, norm


def fig6_phase_diagram_U_J(records: Sequence[EDRecord], outdir: Path, beta: float, t2: float, nb: Optional[int], z_thresh: float) -> None:
    rec = [r for r in records if r.branch == "MET" and approx(r.beta, beta, 1e-6) and approx(r.t2, t2, 1e-6)]
    if nb is not None:
        rec = [r for r in rec if r.nb == nb]
    rec = best_by_grid(rec)
    if len(rec) < 3:
        print("[skip] Figure 6: not enough Scan C records")
        return

    U = np.array([r.U for r in rec])
    Jr = np.array([r.Jr for r in rec])
    phase = np.array([classify_phase(r.Z1, r.Z2, z_thresh) for r in rec])
    cmap, norm = phase_cmap_norm()

    fig, ax = plt.subplots(figsize=(3.75, 3.15))
    sc = ax.scatter(U, Jr, c=phase, cmap=cmap, norm=norm, s=150, marker="s", edgecolor="k", linewidth=0.4)
    # Also show labels with (Z1,Z2) for sparse debugging.
    for r, ph in zip(rec, phase):
        ax.text(r.U, r.Jr, f"{r.Z1:.2f}\n{r.Z2:.2f}", ha="center", va="center", fontsize=5.5, color="white" if ph in [0, 2] else "black")
    ax.set_xlabel(r"$U/t_1$")
    ax.set_ylabel(r"$J/U$")
    ax.set_title(fr"Phase map, $Z_c={z_thresh:g}$")
    handles = [
        Line2D([0], [0], marker="s", color="w", label=PHASE_LABELS[i], markerfacecolor=cmap(norm(i)), markeredgecolor="k", markersize=8)
        for i in [0, 1, 2]
    ]
    ax.legend(handles=handles, frameon=False, loc="upper left")
    fig.tight_layout()
    savefig(fig, outdir, "fig6_scanC_U_J_phase_diagram")


def fig7_phase_diagram_U_T(records: Sequence[EDRecord], outdir: Path, t2: float, nb: Optional[int], Jr_target: float, z_thresh: float) -> None:
    rec = [r for r in records if r.branch == "MET" and approx(r.t2, t2, 1e-6) and approx(round(r.Jr, 3), Jr_target, 2e-3)]
    if nb is not None:
        rec = [r for r in rec if r.nb == nb]
    rec = best_by_grid(rec)
    if len(rec) < 3:
        print("[skip] Figure 7: not enough Scan D records")
        return

    phase = np.array([classify_phase(r.Z1, r.Z2, z_thresh) for r in rec])
    cmap, norm = phase_cmap_norm()

    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.9))
    ax = axes[0]
    ax.scatter([r.U for r in rec], [r.T for r in rec], c=phase, cmap=cmap, norm=norm, s=145, marker="s", edgecolor="k", linewidth=0.4)
    ax.set_xlabel(r"$U/t_1$")
    ax.set_ylabel(r"$T/t_1=1/\beta$")
    ax.set_title(fr"$J/U={Jr_target:.2f}$, $Z_c={z_thresh:g}$")
    handles = [
        Line2D([0], [0], marker="s", color="w", label=PHASE_LABELS[i], markerfacecolor=cmap(norm(i)), markeredgecolor="k", markersize=8)
        for i in [0, 1, 2]
    ]
    ax.legend(handles=handles, frameon=False, loc="best")

    # Right panel: lowest temperature / highest beta cut.
    max_beta = max(r.beta for r in rec)
    rr = sorted([r for r in rec if approx(r.beta, max_beta, 1e-6)], key=lambda x: x.U)
    ax = axes[1]
    if rr:
        ax.plot([r.U for r in rr], [r.Z1 for r in rr], "o-", color=COL_ORB1, label=r"$Z_1$")
        ax.plot([r.U for r in rr], [r.Z2 for r in rr], "s-", color=COL_ORB2, label=r"$Z_2$")
        ax.axhline(z_thresh, color="0.5", ls="--", lw=0.8, label=r"$Z_c$")
    ax.set_xlabel(r"$U/t_1$")
    ax.set_ylabel(r"$Z_m$")
    ax.set_title(fr"lowest $T$: $\beta={max_beta:g}$")
    ax.legend(frameon=False)
    fig.tight_layout()
    savefig(fig, outdir, "fig7_scanD_U_T_phase_diagram")


# -----------------------------------------------------------------------------
# Extra QA plots useful for production scans
# -----------------------------------------------------------------------------


def plot_convergence_quality(records: Sequence[EDRecord], outdir: Path) -> None:
    if not records:
        return
    rec = best_by_grid(records)
    fig, axes = plt.subplots(1, 2, figsize=(6.75, 2.8))

    ax = axes[0]
    vals = np.array([r.last_delta for r in rec])
    x = np.arange(len(rec))
    mask = np.isfinite(vals)
    if mask.any():
        ax.semilogy(x[mask], vals[mask], "o", color=COL_GREY)
    ax.set_xlabel("loaded point index")
    ax.set_ylabel(r"final RMS $\Delta\Sigma$")
    ax.set_title("Convergence quality")
    if len(rec) == 1:
        ax.set_xlim(-0.5, 0.5)

    ax = axes[1]
    Uvals = np.array([r.U for r in rec])
    n1dev = np.array([r.n1 - 1.0 for r in rec])
    n2dev = np.array([r.n2 - 1.0 for r in rec])
    ax.plot(Uvals, n1dev, "o", color=COL_ORB1, label=r"$n_1-1$")
    ax.plot(Uvals, n2dev, "s", color=COL_ORB2, label=r"$n_2-1$")
    ax.axhline(0.0, color="0.5", lw=0.8)
    ax.set_xlabel(r"$U/t_1$")
    ax.set_ylabel(r"$n_m-1$")
    ax.set_title("Half-filling check")
    if len(rec) == 1:
        u = float(Uvals[0])
        ax.set_xlim(u - 0.5, u + 0.5)
        ymax = max(1e-12, float(np.nanmax(np.abs([n1dev[0], n2dev[0]]))) * 2.5)
        ax.set_ylim(-ymax, ymax)
    ax.legend(frameon=False)
    fig.tight_layout()
    savefig(fig, outdir, "diagnostic_convergence_and_filling")


def main() -> int:
    p = argparse.ArgumentParser(description="Plot ED-DMFT smoke diagnostics and production figures.")
    p.add_argument("--mode", choices=["smoke", "production", "all"], default="all")
    p.add_argument("--data", nargs="*", default=[], help="Generic data directories/files/globs; used for smoke or all.")
    p.add_argument("--scanA", nargs="*", default=[], help="Scan A dirs/files: J=0 degenerate U sweep.")
    p.add_argument("--scanB", nargs="*", default=[], help="Scan B dirs/files: fixed U, J/U sweep.")
    p.add_argument("--scanC", nargs="*", default=[], help="Scan C dirs/files: U,J phase map, t2<1.")
    p.add_argument("--scanD", nargs="*", default=[], help="Scan D dirs/files: U,T phase map.")
    p.add_argument("--logs", nargs="*", default=[], help="Optional log files for convergence-history plots.")
    p.add_argument("--out", default="figures_ed_dmft", help="Output figure directory.")
    p.add_argument("--usetex", action="store_true", help="Use LaTeX rendering. Disable if cluster lacks LaTeX.")
    p.add_argument("--beta", type=float, default=25.0, help="Default beta for Scan A/B/C figures.")
    p.add_argument("--nb", type=int, default=None, help="Filter by bath levels. Omit to use all.")
    p.add_argument("--t2", type=float, default=None, help="Generic t2 filter for smoke/Fig.1 when needed.")
    p.add_argument("--t2A", type=float, default=1.0)
    p.add_argument("--t2B", type=float, default=1.0)
    p.add_argument("--t2C", type=float, default=0.5)
    p.add_argument("--t2D", type=float, default=0.5)
    p.add_argument("--U_hund", type=float, default=5.0, help="Fixed U for Fig.5 Scan B.")
    p.add_argument("--U_osmt", type=float, default=3.0, help="Target U for Fig.4 OSMT self-energy.")
    p.add_argument("--Jr_osmt", type=float, default=0.10, help="Target J/U for Fig.4 OSMT self-energy.")
    p.add_argument("--Jr_D", type=float, default=0.10, help="Target J/U for Fig.7 U-T map.")
    p.add_argument("--z_thresh", type=float, default=0.05, help="Z threshold for phase classification.")
    args = p.parse_args()

    setup_matplotlib(args.usetex)
    outdir = Path(args.out).expanduser()
    outdir.mkdir(parents=True, exist_ok=True)

    # Load data. Keep scan pools separate so figures do not contaminate each other.
    generic_records = load_records(args.data)
    scanA = load_records(args.scanA) if args.scanA else []
    scanB = load_records(args.scanB) if args.scanB else []
    scanC = load_records(args.scanC) if args.scanC else []
    scanD = load_records(args.scanD) if args.scanD else []

    all_records = best_by_grid(generic_records + scanA + scanB + scanC + scanD)
    write_summary_csv(all_records, outdir)
    plot_convergence_quality(all_records, outdir)

    if args.mode in {"smoke", "all"}:
        smoke_records = generic_records or all_records
        if args.nb is not None:
            smoke_records = [r for r in smoke_records if r.nb == args.nb]
        if args.t2 is not None:
            smoke_records = [r for r in smoke_records if approx(r.t2, args.t2, 1e-6)]
        plot_smoke_dashboard(smoke_records, outdir)
        for i, log_path in enumerate(args.logs):
            log = parse_log_file(log_path)
            if log is not None:
                plot_log_convergence(log, outdir, name=f"log_convergence_{i+1}_{Path(log_path).stem}")

    if args.mode in {"production", "all"}:
        # If specific scan dirs were not given, fall back to all generic records.
        A = scanA if scanA else all_records
        B = scanB if scanB else all_records
        C = scanC if scanC else all_records
        D = scanD if scanD else all_records

        fig1_green_validation(all_records, outdir, beta=args.beta if args.beta else None, t2=args.t2, nb=args.nb)
        fig2_mott_scan_A(A, outdir, beta=args.beta, t2=args.t2A, nb=args.nb)
        fig3_Z_vs_U_scan_C(C, outdir, beta=args.beta, t2=args.t2C, nb=args.nb)
        fig4_self_energy_osmt(C, outdir, beta=args.beta, t2=args.t2C, nb=args.nb, U_osmt=args.U_osmt, Jr_osmt=args.Jr_osmt)
        fig5_hund_docc_scan_B(B, outdir, beta=args.beta, t2=args.t2B, nb=args.nb, U_fixed=args.U_hund)
        fig6_phase_diagram_U_J(C, outdir, beta=args.beta, t2=args.t2C, nb=args.nb, z_thresh=args.z_thresh)
        fig7_phase_diagram_U_T(D, outdir, t2=args.t2D, nb=args.nb, Jr_target=args.Jr_D, z_thresh=args.z_thresh)

    print("\nDone. Generated figures in:", outdir.resolve())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
