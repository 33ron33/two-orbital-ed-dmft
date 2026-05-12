#!/usr/bin/env python3
import os
import subprocess
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed

OUTDIR = Path("data/nb2_scanC_t2p30_coarse")
LOGDIR = Path("logs")
OUTDIR.mkdir(parents=True, exist_ok=True)
LOGDIR.mkdir(parents=True, exist_ok=True)

PYTHON = os.environ.get("ED_PYTHON", "python")

# Coarse grid: not too expensive, enough to see orbital differentiation.
U_GRID = [0.50, 1.00, 1.50, 2.00, 2.50, 3.00, 3.50, 4.00, 4.50, 5.00, 5.50, 6.00]
JR_GRID = [0.10, 0.20, 0.30]

BASE = [
    PYTHON, "cluster_dmft_ed_updated.py",
    "--beta", "25",
    "--t1", "1.0",
    "--t2", "0.3",
    "--Nb", "2",
    "--N_IW", "160",
    "--n_iw_fit", "70",
    "--conv_n_iw", "30",
    "--n_iter", "100",
    "--tol", "5e-3",
    "--mix", "0.10",
    "--fit_restarts", "2",
    "--seed", "metal",
    "--ph_sym_bath",
    "--outdir", str(OUTDIR),
]

def expected_file(U, J):
    return OUTDIR / f"ED_MET_U{U:.2f}_J{J:.3f}_t20.30_b25.0_Nb2.npz"

def run_jr_branch(jr):
    prev = None

    print("=" * 90, flush=True)
    print(f"Starting Nb=2, t2=0.3 branch J/U={jr:.2f}", flush=True)
    print("=" * 90, flush=True)

    for U in U_GRID:
        J = U * jr
        Utag = f"{U:.2f}"
        Jtag = f"{J:.6f}"
        logtag = f"nb2_t2p30_JR{jr:.2f}_U{U:.2f}".replace(".", "p")
        logfile = LOGDIR / f"{logtag}.log"

        cmd = BASE + ["--U", Utag, "--J", Jtag]

        # Metallic continuation: increasing U for each fixed J/U.
        if prev is not None:
            cmd += ["--init_from", str(prev)]

        print("-" * 90, flush=True)
        print(f"J/U={jr:.2f}, U={Utag}, J={J:.6f}", flush=True)
        print(f"prev = {prev}", flush=True)
        print(" ".join(cmd), flush=True)
        print(f"log = {logfile}", flush=True)
        print("-" * 90, flush=True)

        with open(logfile, "w") as f:
            ret = subprocess.run(cmd, stdout=f, stderr=subprocess.STDOUT)

        if ret.returncode != 0:
            raise RuntimeError(f"Failed J/U={jr:.2f}, U={Utag}. Check {logfile}")

        out = expected_file(U, J)
        if not out.exists():
            raise FileNotFoundError(f"Expected output not found: {out}")

        prev = out

    print("=" * 90, flush=True)
    print(f"Finished branch J/U={jr:.2f}", flush=True)
    print("=" * 90, flush=True)

def main():
    parallel_jr = int(os.environ.get("PARALLEL_JR", "2"))

    failures = []
    with ThreadPoolExecutor(max_workers=parallel_jr) as ex:
        futs = {ex.submit(run_jr_branch, jr): jr for jr in JR_GRID}
        for fut in as_completed(futs):
            jr = futs[fut]
            try:
                fut.result()
            except Exception as e:
                failures.append((jr, str(e)))
                print(f"[FAIL] J/U={jr:.2f}: {e}", flush=True)

    if failures:
        print("Failures:")
        for jr, msg in failures:
            print(f"  J/U={jr:.2f}: {msg}")
        raise SystemExit(1)

    print("All Nb=2 t2=0.3 coarse Scan C branches completed.")

if __name__ == "__main__":
    main()
