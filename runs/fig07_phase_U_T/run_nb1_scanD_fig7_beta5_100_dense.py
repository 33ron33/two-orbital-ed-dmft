#!/usr/bin/env python3
import os
import glob
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np

OUTDIR = Path("data/nb1_scanD_fig7_beta5_100_dense")
LOGDIR = Path("logs")
OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

PYTHON = os.environ.get("ED_PYTHON", "python")

# ============================================================
# Dense Figure 7 scan
# ============================================================
BETAS = [float(b) for b in range(5, 101, 5)]          # 5, 10, ..., 100
U_GRID = [round(float(x), 4) for x in np.arange(0.0, 6.0001, 0.125)]
JR = 0.10
T2 = 0.3

def solver_mesh_settings(beta):
    """
    Use larger Matsubara grids and more iterations at lower temperatures.
    """
    if beta <= 20:
        return dict(N_IW=160, n_iw_fit=70, conv_n_iw=30, n_iter=180, tol="5e-3", mix="0.08")
    elif beta <= 50:
        return dict(N_IW=240, n_iw_fit=100, conv_n_iw=40, n_iter=240, tol="4e-3", mix="0.06")
    else:
        return dict(N_IW=360, n_iw_fit=140, conv_n_iw=60, n_iter=350, tol="3e-3", mix="0.045")

def load_meta(npzfile):
    d = np.load(npzfile)
    return {
        "U": float(d["U"]),
        "J": float(d["J"]),
        "beta": float(d["beta"]),
        "converged": bool(d["converged"]),
        "last_delta": float(d["last_delta"]),
        "Z1": float(d["Z1"]),
        "Z2": float(d["Z2"]),
        "n1": float(d["n1"]),
        "n2": float(d["n2"]),
    }

def find_point(U_target, beta_target, J_target, tol_U=1e-3, tol_beta=1e-8, tol_J=1e-3):
    files = sorted(glob.glob(str(OUTDIR / "*Nb1.npz")))
    matches = []

    for f in files:
        try:
            m = load_meta(f)
        except Exception:
            continue

        if (
            abs(m["U"] - U_target) < tol_U
            and abs(m["beta"] - beta_target) < tol_beta
            and abs(m["J"] - J_target) < tol_J
        ):
            matches.append(f)

    if not matches:
        return None

    matches = sorted(matches, key=lambda x: Path(x).stat().st_mtime)
    return Path(matches[-1])

def find_previous_seed(U_target, beta_target, jr=JR, tol_jr=1.2e-3):
    """
    Use nearest lower-U point at same beta and same J/U as warm start.
    """
    files = sorted(glob.glob(str(OUTDIR / "*Nb1.npz")))
    candidates = []

    for f in files:
        try:
            m = load_meta(f)
            U = m["U"]
            J = m["J"]
            beta = m["beta"]
            Jr_file = J / U if U != 0 else 0.0
        except Exception:
            continue

        if U < U_target and abs(beta - beta_target) < 1e-8 and abs(Jr_file - jr) < tol_jr:
            candidates.append((U, f))

    if not candidates:
        return None

    candidates = sorted(candidates, key=lambda x: x[0])
    return Path(candidates[-1][1])

def run_beta_branch(beta):
    prev = None
    mesh = solver_mesh_settings(beta)

    print("=" * 96, flush=True)
    print(f"Starting beta branch beta={beta:.1f}, T={1.0/beta:.6f}", flush=True)
    print(f"settings: {mesh}", flush=True)
    print("=" * 96, flush=True)

    for U in U_GRID:
        J = JR * U

        existing = find_point(U, beta, J)
        if existing is not None:
            print(f"[skip exists] beta={beta:.1f}, U={U:.4f}, J={J:.6f} -> {existing}", flush=True)
            prev = existing
            continue

        if prev is None:
            prev = find_previous_seed(U, beta)

        cmd = [
            PYTHON, "solver/ed_solver_dmft.py",
            "--beta", f"{beta:.1f}",
            "--U", f"{U:.6f}",
            "--J", f"{J:.6f}",
            "--t1", "1.0",
            "--t2", f"{T2:.1f}",
            "--Nb", "1",
            "--N_IW", str(mesh["N_IW"]),
            "--n_iw_fit", str(mesh["n_iw_fit"]),
            "--conv_n_iw", str(mesh["conv_n_iw"]),
            "--n_iter", str(mesh["n_iter"]),
            "--tol", mesh["tol"],
            "--mix", mesh["mix"],
            "--fit_restarts", "3",
            "--seed", "metal",
            "--outdir", str(OUTDIR),
        ]

        if prev is not None:
            cmd += ["--init_from", str(prev)]

        tag = f"nb1_fig7_b{beta:.1f}_U{U:.4f}_JR{JR:.2f}".replace(".", "p")
        logfile = LOGDIR / f"{tag}.log"

        print("-" * 96, flush=True)
        print(f"beta={beta:.1f}, T={1.0/beta:.6f}, U={U:.4f}, J={J:.6f}", flush=True)
        print(f"prev = {prev}", flush=True)
        print(" ".join(cmd), flush=True)
        print(f"log = {logfile}", flush=True)
        print("-" * 96, flush=True)

        with open(logfile, "w") as f:
            ret = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

        if ret.returncode != 0:
            raise RuntimeError(f"Solver failed beta={beta:.1f}, U={U:.4f}. Check {logfile}")

        out = find_point(U, beta, J)
        if out is None:
            raise FileNotFoundError(
                f"Finished run but output not found for beta={beta:.1f}, U={U:.4f}, J={J:.6f}"
            )

        prev = out
        print(f"[saved] {out}", flush=True)

    print("=" * 96, flush=True)
    print(f"Finished beta branch beta={beta:.1f}", flush=True)
    print("=" * 96, flush=True)

def main():
    parallel_beta = int(os.environ.get("PARALLEL_BETA", "4"))

    print("=" * 96)
    print("Dense Nb=1 Figure 7 scan: beta=5..100, U=0..6")
    print(f"betas      : {BETAS}")
    print(f"U points   : {len(U_GRID)}")
    print(f"U range    : {U_GRID[0]} ... {U_GRID[-1]}")
    print(f"J/U        : {JR:.2f}")
    print(f"t2/t1      : {T2:.1f}")
    print(f"target pts : {len(BETAS) * len(U_GRID)}")
    print(f"outdir     : {OUTDIR}")
    print("=" * 96)

    failures = []
    with ThreadPoolExecutor(max_workers=parallel_beta) as ex:
        futs = {ex.submit(run_beta_branch, beta): beta for beta in BETAS}
        for fut in as_completed(futs):
            beta = futs[fut]
            try:
                fut.result()
            except Exception as e:
                failures.append((beta, str(e)))
                print(f"[FAIL] beta={beta:.1f}: {e}", flush=True)

    if failures:
        print("\nFailures:")
        for beta, msg in failures:
            print(f"  beta={beta:.1f}: {msg}")
        raise SystemExit(1)

    print("\nAll dense Figure 7 beta branches completed successfully.")

if __name__ == "__main__":
    main()
