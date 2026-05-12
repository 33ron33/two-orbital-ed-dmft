#!/usr/bin/env python3
"""
cluster_dmft_ed_updated.py

Single-point two-orbital Kanamori ED-DMFT for a Bethe lattice.
Designed for cluster/SLURM sweeps where each job runs one (U, J, seed) point.

Main fixes compared with the earlier script:
  1. Correct half-filled Kanamori chemical potential in Weiss/Dyson:
        G0^{-1}(iw) = iw + mu_H - Delta(iw)
        G0^{-1,new}(iw) = iw + mu_H - t_m^2 G_m(iw)
     where mu_H = (3U - 5J)/2.
  2. Real DMFT stabilization by mixing fitted bath parameters, not only Sigma.
  3. Optional particle-hole-symmetric bath fitting for half-filled ED.
  4. Better convergence metric: RMS change of low-frequency Sigma.
  5. Default exit code is 0 even if not converged, so SLURM scans do not mark
     saved-but-unconverged points as failed. Use --fail_on_nonconvergence if needed.

Usage example:
  python cluster_dmft_ed_updated.py --U 5.0 --J 0.75 --beta 25 \
      --t1 1.0 --t2 0.5 --Nb 2 --seed metal --ph_sym_bath \
      --outdir data/results
"""

import argparse
import datetime
import os
import sys
import time
import warnings

import numpy as np
from scipy.optimize import minimize

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
LOCAL_PYED = os.path.join(PROJECT_DIR, "pyed")
if os.path.isdir(os.path.join(LOCAL_PYED, "pyed")) and LOCAL_PYED not in sys.path:
    sys.path.insert(0, LOCAL_PYED)

from triqs.operators import c, c_dag, Operator
from triqs.gf import Gf, MeshImFreq
from pyed.TriqsExactDiagonalization import TriqsExactDiagonalization

# Different xargs/SLURM workers should not get identical bath-fit random starts.
np.random.seed(os.getpid() % (2**31))

UP, DN = 0, 1


# -----------------------------------------------------------------------------
# Basic helpers
# -----------------------------------------------------------------------------

def kanamori_mu_half_filling(U, J):
    """Two-orbital Kanamori chemical potential at half filling."""
    return 0.5 * (3.0 * U - 5.0 * J)


def ed_size_summary(N_b):
    """
    ED Hilbert-space size estimate.

    Two impurity orbitals, two spins = 4 impurity modes.
    Each orbital has N_b bath sites, each with two spins = 4*N_b bath modes.
    """
    n_modes = 4 + 4 * N_b
    n_states = 2 ** n_modes
    dense_complex_gib = (
        n_states * n_states * np.dtype(np.complex128).itemsize / 1024**3
    )
    return n_modes, n_states, dense_complex_gib


def iw_array_from_mesh(mesh):
    return np.array([complex(w) for w in mesh], dtype=np.complex128)


def positive_iw_slice(mesh, n_keep=None):
    n_total = len(list(mesh))
    n0 = n_total // 2
    if n_keep is None:
        return slice(n0, n_total)
    return slice(n0, min(n_total, n0 + int(n_keep)))


def make_gf_like(mesh):
    return Gf(mesh=mesh, target_shape=[1, 1])


# -----------------------------------------------------------------------------
# Hamiltonian and ED impurity solve
# -----------------------------------------------------------------------------

def build_kanamori_aim(U, J, eps_bath, V_bath):
    """
    Build two-orbital Kanamori Anderson impurity model.

    Orbitals 0 and 1 are impurity orbitals.
    Bath orbital index convention:
        bath for impurity orbital m starts at bi = 2 + m*N_b.

    Half-filling is imposed through the impurity chemical potential
        mu_H = (3U - 5J)/2.
    """
    N_b = eps_bath.shape[1]
    Uprime = U - 2.0 * J
    mu_H = kanamori_mu_half_filling(U, J)

    H = Operator()

    # Intra-orbital Hubbard U.
    for m in [0, 1]:
        n_up = c_dag(UP, m) * c(UP, m)
        n_dn = c_dag(DN, m) * c(DN, m)
        H += U * n_up * n_dn

    # Inter-orbital density-density terms.
    n0u = c_dag(UP, 0) * c(UP, 0)
    n0d = c_dag(DN, 0) * c(DN, 0)
    n1u = c_dag(UP, 1) * c(UP, 1)
    n1d = c_dag(DN, 1) * c(DN, 1)

    # Opposite spin inter-orbital: U'.
    H += Uprime * n0u * n1d
    H += Uprime * n1u * n0d

    # Same spin inter-orbital: U' - J.
    H += (Uprime - J) * n0u * n1u
    H += (Uprime - J) * n0d * n1d

    # Spin-flip terms.
    H += -J * (c_dag(UP, 0) * c(DN, 0) * c_dag(DN, 1) * c(UP, 1))
    H += -J * (c_dag(UP, 1) * c(DN, 1) * c_dag(DN, 0) * c(UP, 0))

    # Pair-hopping terms. Keep the same sign convention as the original script.
    H += J * (c_dag(UP, 0) * c_dag(DN, 0) * c(DN, 1) * c(UP, 1))
    H += J * (c_dag(UP, 1) * c_dag(DN, 1) * c(DN, 0) * c(UP, 0))

    # Impurity chemical potential at half filling.
    for m in [0, 1]:
        H -= mu_H * (c_dag(UP, m) * c(UP, m) + c_dag(DN, m) * c(DN, m))

    # Bath levels and hybridization.
    for m in [0, 1]:
        for k in range(N_b):
            bi = 2 + m * N_b + k
            for sp in [UP, DN]:
                H += eps_bath[m, k] * c_dag(sp, bi) * c(sp, bi)
                H += V_bath[m, k] * (
                    c_dag(sp, m) * c(sp, bi) + c_dag(sp, bi) * c(sp, m)
                )

    # Fundamental operators required by pyed.
    fund_ops = [c(UP, m) for m in [0, 1]] + [c(DN, m) for m in [0, 1]]
    for m in [0, 1]:
        for k in range(N_b):
            bi = 2 + m * N_b + k
            fund_ops += [c(UP, bi), c(DN, bi)]

    return H, fund_ops


def ed_solve_2orb(U, J, eps_bath, V_bath, beta, n_iw):
    """Solve one ED impurity problem and return G, Sigma, G0 and observables."""
    H, fund_ops = build_kanamori_aim(U, J, eps_bath, V_bath)
    ed = TriqsExactDiagonalization(H, fund_ops, beta)

    mesh = MeshImFreq(beta, 'Fermion', n_iw)
    iw = iw_array_from_mesh(mesh)
    mu_H = kanamori_mu_half_filling(U, J)

    G_iw = []
    for m in [0, 1]:
        g = make_gf_like(mesh)
        ed.set_g2_iwn(g[0, 0], c(UP, m), c_dag(UP, m))
        G_iw.append(g)

    n0u = c_dag(UP, 0) * c(UP, 0)
    n0d = c_dag(DN, 0) * c(DN, 0)
    n1u = c_dag(UP, 1) * c(UP, 1)
    n1d = c_dag(DN, 1) * c(DN, 1)

    D1 = float(ed.get_expectation_value(n0u * n0d).real)
    D2 = float(ed.get_expectation_value(n1u * n1d).real)
    D12 = float(ed.get_expectation_value((n0u + n0d) * (n1u + n1d)).real)
    occ1 = float(ed.get_expectation_value(n0u + n0d).real)
    occ2 = float(ed.get_expectation_value(n1u + n1d).real)

    G0_iw, Sig_iw = [], []
    for m in [0, 1]:
        hyb = np.sum(
            V_bath[m] ** 2 / (iw[:, None] - eps_bath[m][None, :]), axis=1
        )

        # Important correction: because H contains -mu_H*n_imp,
        # the noninteracting Weiss inverse is iw + mu_H - Delta(iw).
        G0inv = iw + mu_H - hyb

        g0 = make_gf_like(mesh)
        g0.data[:, 0, 0] = 1.0 / G0inv

        sig = make_gf_like(mesh)
        sig.data[:, 0, 0] = G0inv - 1.0 / G_iw[m].data[:, 0, 0]

        G0_iw.append(g0)
        Sig_iw.append(sig)

    return G_iw, Sig_iw, G0_iw, D1, D2, D12, occ1, occ2


# -----------------------------------------------------------------------------
# ED bath fitting and Bethe update
# -----------------------------------------------------------------------------

def _hyb_from_bath(iw, eps, V):
    return np.sum(V**2 / (iw[:, None] - eps[None, :]), axis=1)


def _bath_weiss(iw, eps, V, mu):
    return 1.0 / (iw + mu - _hyb_from_bath(iw, eps, V))


def _expand_ph_bath(x, N_b):
    """
    Expand PH-symmetric variables into full bath arrays.

    For N_b=2:
        x = [e1, v1]
        eps = [-e1, +e1]
        V   = [ v1,  v1]

    For N_b=4:
        x = [e1, e2, v1, v2]
        eps = [-e1, +e1, -e2, +e2]
        V   = [ v1,  v1,  v2,  v2]
    """
    if N_b % 2 != 0:
        raise ValueError("PH-symmetric bath fitting requires even N_b.")

    M = N_b // 2
    e = np.abs(np.asarray(x[:M], dtype=float))
    v = np.abs(np.asarray(x[M:], dtype=float))

    eps = np.empty(N_b, dtype=float)
    V = np.empty(N_b, dtype=float)
    eps[0::2] = -e
    eps[1::2] = +e
    V[0::2] = v
    V[1::2] = v
    return eps, V


def _compress_to_ph_guess(eps, V, N_b, t):
    """Convert an arbitrary bath into a reasonable PH-symmetric initial guess."""
    M = N_b // 2
    e_abs = np.sort(np.abs(np.asarray(eps, dtype=float)))
    v_abs = np.abs(np.asarray(V, dtype=float))

    if len(e_abs) >= N_b:
        e_guess = e_abs[::2][:M]
    else:
        e_guess = np.linspace(0.25 * t, 1.5 * t, M)

    if len(v_abs) >= N_b:
        v_guess = []
        for i in range(M):
            v_guess.append(np.mean(v_abs[2*i:2*i+2]))
        v_guess = np.array(v_guess)
    else:
        v_guess = np.full(M, t / np.sqrt(max(N_b, 1)))

    e_guess = np.clip(e_guess, 1e-4, 5.0 * t)
    v_guess = np.clip(v_guess, 1e-4, 5.0 * t)
    return np.r_[e_guess, v_guess]


def fit_bath(
    G0_target,
    N_b,
    beta,
    n_iw_fit=80,
    t=1.0,
    mu=0.0,
    ph_sym=False,
    prev_eps=None,
    prev_V=None,
    n_restart=6,
):
    """
    Fit ED bath parameters to a target Weiss field.

    If ph_sym=True and N_b is even, fit only PH-symmetric bath pairs.
    The previous bath is used as the first optimizer start when available.
    """
    iw_all = iw_array_from_mesh(G0_target.mesh)
    pos = positive_iw_slice(G0_target.mesh, n_iw_fit)
    iw_fit = iw_all[pos]
    G0_tgt = G0_target.data[pos, 0, 0]

    # Low-frequency weighted residual. ED bath fitting should care most about
    # the low Matsubara region for DMFT stability and Z/self-energy behavior.
    j = np.arange(len(iw_fit), dtype=float)
    weights = 1.0 / (1.0 + (j / 20.0) ** 2)

    best_x = None
    best_f = np.inf

    if ph_sym:
        if N_b % 2 != 0:
            raise ValueError("ph_sym=True requires even N_b")
        M = N_b // 2

        def unpack(x):
            return _expand_ph_bath(x, N_b)

        def residual(x):
            eps, V = unpack(x)
            d = _bath_weiss(iw_fit, eps, V, mu) - G0_tgt
            return float(np.sum(weights * (d.real**2 + d.imag**2)))

        bounds = [(1e-4, 5.0 * t)] * M + [(1e-4, 5.0 * t)] * M

        starts = []
        if prev_eps is not None and prev_V is not None:
            starts.append(_compress_to_ph_guess(prev_eps, prev_V, N_b, t))

        # Deterministic starting guesses are useful for reproducible scans.
        starts.append(np.r_[np.linspace(0.25 * t, 1.50 * t, M),
                            np.full(M, t / np.sqrt(N_b))])

        for _ in range(max(0, n_restart - len(starts))):
            e0 = np.random.uniform(0.05 * t, 2.5 * t, M)
            v0 = np.abs(np.random.normal(t / np.sqrt(N_b), 0.2 * t, M)) + 1e-3
            starts.append(np.r_[e0, v0])

        for x0 in starts:
            res = minimize(
                residual,
                x0,
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': 3000, 'ftol': 1e-11, 'gtol': 1e-8},
            )
            if res.fun < best_f:
                best_f = float(res.fun)
                best_x = np.array(res.x, dtype=float)

        eps_best, V_best = unpack(best_x)
        return eps_best, np.abs(V_best), best_f

    # General, non-PH-constrained bath.
    def residual(x):
        eps = x[:N_b]
        V = np.abs(x[N_b:])
        d = _bath_weiss(iw_fit, eps, V, mu) - G0_tgt
        return float(np.sum(weights * (d.real**2 + d.imag**2)))

    bounds = [(-5.0 * t, 5.0 * t)] * N_b + [(1e-4, 5.0 * t)] * N_b

    starts = []
    if prev_eps is not None and prev_V is not None:
        starts.append(np.r_[prev_eps, np.abs(prev_V)])

    starts.append(np.r_[np.zeros(N_b), np.full(N_b, t / np.sqrt(N_b))])

    for _ in range(max(0, n_restart - len(starts))):
        e0 = np.random.uniform(-2.0 * t, 2.0 * t, N_b)
        v0 = np.abs(np.random.normal(t / np.sqrt(N_b), 0.2 * t, N_b)) + 1e-3
        starts.append(np.r_[e0, v0])

    for x0 in starts:
        res = minimize(
            residual,
            x0,
            method='L-BFGS-B',
            bounds=bounds,
            options={'maxiter': 3000, 'ftol': 1e-11, 'gtol': 1e-8},
        )
        if res.fun < best_f:
            best_f = float(res.fun)
            best_x = np.array(res.x, dtype=float)

    eps_best = best_x[:N_b]
    V_best = np.abs(best_x[N_b:])
    return eps_best, V_best, best_f


def bethe_update_G0(G_iw, U, J, t_orb=1.0):
    """
    Bethe lattice self-consistency for one orbital:
        G0^{-1}(iw) = iw + mu_H - t_orb^2 G(iw)
    """
    iw = iw_array_from_mesh(G_iw.mesh)
    mu_H = kanamori_mu_half_filling(U, J)
    G0inv = iw + mu_H - t_orb**2 * G_iw.data[:, 0, 0]

    g0 = make_gf_like(G_iw.mesh)
    g0.data[:, 0, 0] = 1.0 / G0inv
    return g0


# -----------------------------------------------------------------------------
# Observables and convergence
# -----------------------------------------------------------------------------

def quasiparticle_Z(Sigma_iw, beta):
    """Eliashberg estimate Z = 1/[1 - Im Sigma(iw0)/w0]."""
    n0 = len(list(Sigma_iw.mesh)) // 2
    ImS0 = float(Sigma_iw.data[n0, 0, 0].imag)
    omega0 = np.pi / beta
    denom = 1.0 - ImS0 / omega0
    return float(1.0 / denom) if abs(denom) > 1e-12 else np.nan


def sigma_rms_delta(Sig_new, Sig_old, n_low=40):
    """RMS low-frequency self-energy change, max over two orbitals."""
    vals = []
    for m in [0, 1]:
        pos = positive_iw_slice(Sig_new[m].mesh, n_low)
        d = Sig_new[m].data[pos, 0, 0] - Sig_old[m].data[pos, 0, 0]
        vals.append(float(np.sqrt(np.mean(np.abs(d) ** 2))))
    return max(vals)


def initialize_bath(U, t1, t2, N_b, seed='metal', ph_sym=False):
    """Initial bath guess for metal/insulator branch."""
    t_orbs = [t1, t2]
    eps_bath = np.zeros((2, N_b), dtype=float)
    V_bath = np.zeros((2, N_b), dtype=float)

    for m, tm in enumerate(t_orbs):
        if seed == 'insulator':
            if ph_sym and N_b % 2 == 0:
                M = N_b // 2
                e = np.linspace(max(0.25 * tm, U / 4.0), max(0.5 * tm, U / 2.0), M)
                eps_bath[m, 0::2] = -e
                eps_bath[m, 1::2] = +e
                V_bath[m, :] = 0.5 * tm / np.sqrt(N_b)
            else:
                half = max(N_b // 2, 1)
                eps_row = np.concatenate([
                    np.linspace(-U / 2.0, -U / 4.0, N_b - half),
                    np.linspace(+U / 4.0, +U / 2.0, half),
                ])
                eps_bath[m, :] = eps_row
                V_bath[m, :] = 0.5 * tm / np.sqrt(N_b)
        else:
            if ph_sym and N_b % 2 == 0:
                M = N_b // 2
                e = np.linspace(0.25 * tm, 1.25 * tm, M)
                eps_bath[m, 0::2] = -e
                eps_bath[m, 1::2] = +e
            else:
                eps_bath[m, :] = 0.0
            V_bath[m, :] = tm / np.sqrt(N_b)

    return eps_bath, V_bath


# -----------------------------------------------------------------------------
# DMFT loop
# -----------------------------------------------------------------------------

def run_dmft_2orb(
    U,
    J,
    beta,
    t1=1.0,
    t2=1.0,
    N_b=2,
    n_iw=200,
    n_iter=60,
    tol=1e-4,
    mix=0.5,
    seed='metal',
    enforce_sym=False,
    ph_sym_bath=False,
    n_iw_fit=80,
    conv_n_iw=40,
    fit_restarts=6,
    verbose=True,
    eps_init=None,
    V_init=None,
):
    """Run one two-orbital ED-DMFT point."""
    if ph_sym_bath and N_b % 2 != 0:
        raise ValueError("--ph_sym_bath requires an even --Nb, e.g. --Nb 2.")

    t_orbs = [t1, t2]
    mu_H = kanamori_mu_half_filling(U, J)

    if eps_init is not None and V_init is not None:
        eps_bath = np.array(eps_init, dtype=float).copy()
        V_bath = np.array(V_init, dtype=float).copy()
        if eps_bath.shape != (2, N_b) or V_bath.shape != (2, N_b):
            raise ValueError(
                f"Warm-start bath shape mismatch. Expected {(2, N_b)}, "
                f"got eps {eps_bath.shape}, V {V_bath.shape}."
            )
    else:
        eps_bath, V_bath = initialize_bath(U, t1, t2, N_b, seed=seed, ph_sym=ph_sym_bath)

    Sig_old = None
    G_iw = Sig_iw = None
    D1 = D2 = D12 = occ1 = occ2 = np.nan
    converged = False
    last_delta = np.inf
    fit_err = [np.nan, np.nan]

    for it in range(n_iter):
        t_it = time.time()

        G_iw, Sig_iw, _, D1, D2, D12, occ1, occ2 = ed_solve_2orb(
            U, J, eps_bath, V_bath, beta, n_iw
        )

        Z1 = quasiparticle_Z(Sig_iw[0], beta)
        Z2 = quasiparticle_Z(Sig_iw[1], beta)

        if Sig_old is not None:
            last_delta = sigma_rms_delta(Sig_iw, Sig_old, n_low=conv_n_iw)
            if verbose:
                print(
                    f"  it={it:3d}  dSigma_rms={last_delta:.4e}  "
                    f"D1={D1:.6f} D2={D2:.6f}  "
                    f"Z1={Z1:.4f} Z2={Z2:.4f}  "
                    f"n1={occ1:.4f} n2={occ2:.4f}  "
                    f"fit=({fit_err[0]:.2e},{fit_err[1]:.2e})  "
                    f"({time.time() - t_it:.0f}s)",
                    flush=True,
                )
            if last_delta < tol:
                converged = True
                break
        else:
            if verbose:
                print(
                    f"  it={it:3d}  init  "
                    f"D1={D1:.6f} D2={D2:.6f}  "
                    f"Z1={Z1:.4f} Z2={Z2:.4f}  "
                    f"n1={occ1:.4f} n2={occ2:.4f}  "
                    f"({time.time() - t_it:.0f}s)",
                    flush=True,
                )

        # Bethe self-consistency builds the new target Weiss field.
        G0_new = [bethe_update_G0(G_iw[m], U, J, t_orbs[m]) for m in [0, 1]]

        # Only use orbital symmetry when t1 == t2 and the physics truly has equivalent orbitals.
        if enforce_sym:
            avg = 0.5 * (G0_new[0].data[:, 0, 0] + G0_new[1].data[:, 0, 0])
            for m in [0, 1]:
                G0_new[m].data[:, 0, 0] = avg

        # Fit new ED bath and mix bath parameters directly.
        eps_next = np.zeros_like(eps_bath)
        V_next = np.zeros_like(V_bath)
        fit_err = [np.nan, np.nan]

        for m in [0, 1]:
            eps_fit, V_fit, ferr = fit_bath(
                G0_new[m],
                N_b=N_b,
                beta=beta,
                n_iw_fit=n_iw_fit,
                t=t_orbs[m],
                mu=mu_H,
                ph_sym=ph_sym_bath,
                prev_eps=eps_bath[m],
                prev_V=V_bath[m],
                n_restart=fit_restarts,
            )
            fit_err[m] = ferr

            if it > 0 and mix < 1.0:
                eps_next[m] = mix * eps_fit + (1.0 - mix) * eps_bath[m]
                V_next[m] = mix * V_fit + (1.0 - mix) * V_bath[m]
            else:
                eps_next[m] = eps_fit
                V_next[m] = V_fit

            # Re-enforce exact PH pair structure after linear mixing.
            if ph_sym_bath:
                eps_next[m], V_next[m] = _expand_ph_bath(
                    _compress_to_ph_guess(eps_next[m], V_next[m], N_b, t_orbs[m]), N_b
                )

        eps_bath = eps_next
        V_bath = np.abs(V_next)
        Sig_old = Sig_iw

    Z1 = quasiparticle_Z(Sig_iw[0], beta)
    Z2 = quasiparticle_Z(Sig_iw[1], beta)

    return {
        'G_iw': G_iw,
        'Sig_iw': Sig_iw,
        'D1': D1,
        'D2': D2,
        'D12': D12,
        'n1': occ1,
        'n2': occ2,
        'Z1': Z1,
        'Z2': Z2,
        'eps_bath': eps_bath,
        'V_bath': V_bath,
        'converged': converged,
        'n_iter_done': it + 1,
        'last_delta': last_delta,
        'last_fit_err': np.array(fit_err, dtype=float),
    }


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description='Two-orbital Kanamori ED-DMFT, single cluster point.'
    )
    parser.add_argument('--U', type=float, required=True, help='Intra-orbital Hubbard U')
    parser.add_argument('--J', type=float, required=True, help='Hund coupling J')
    parser.add_argument('--beta', type=float, default=25.0, help='Inverse temperature')
    parser.add_argument('--t1', type=float, default=1.0, help='Hopping orbital 1')
    parser.add_argument('--t2', type=float, default=1.0, help='Hopping orbital 2')
    parser.add_argument('--Nb', type=int, default=2, help='Bath levels per orbital')
    parser.add_argument('--N_IW', type=int, default=200, help='Number of positive Matsubara frequencies')
    parser.add_argument('--n_iw_fit', type=int, default=80, help='Positive Matsubara points used in bath fit')
    parser.add_argument('--conv_n_iw', type=int, default=40, help='Positive Matsubara points used for convergence')
    parser.add_argument('--n_iter', type=int, default=60, help='Max DMFT iterations')
    parser.add_argument('--tol', type=float, default=1e-4, help='Convergence tolerance')
    parser.add_argument('--mix', type=float, default=0.5, help='Bath mixing parameter')
    parser.add_argument('--fit_restarts', type=int, default=6, help='Number of bath fit starts per orbital')
    parser.add_argument('--seed', type=str, default='metal', choices=['metal', 'insulator'])
    parser.add_argument('--outdir', type=str, default='data/results')
    parser.add_argument('--enforce_sym', action='store_true', help='Force orbital symmetry; use only when t1 == t2')
    parser.add_argument('--ph_sym_bath', action='store_true', help='Use particle-hole symmetric ED bath fit; requires even Nb')
    parser.add_argument('--init_from', type=str, default=None, help='Warm-start from previous .npz file')
    parser.add_argument('--fail_on_nonconvergence', action='store_true', help='Return exit code 1 if not converged')
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    if not (0.0 < args.mix <= 1.0):
        raise ValueError('--mix must satisfy 0 < mix <= 1')

    if args.ph_sym_bath and args.Nb % 2 != 0:
        raise ValueError('--ph_sym_bath requires even --Nb, e.g. --Nb 2.')

    branch = 'MET' if args.seed == 'metal' else 'INS'
    enforce_sym = args.enforce_sym or (abs(args.t1 - args.t2) < 1e-12)

    if args.enforce_sym and abs(args.t1 - args.t2) > 1e-12:
        warnings.warn(
            '--enforce_sym was passed even though t1 != t2. This removes orbital differentiation.',
            RuntimeWarning,
        )

    eps_init = V_init = None
    if args.init_from:
        if not os.path.exists(args.init_from):
            print(f'WARNING: --init_from file not found: {args.init_from}. Using seed instead.', flush=True)
        else:
            prev = np.load(args.init_from)
            eps_init = prev['eps_bath']
            V_init = prev['V_bath']
            print(f'Warm-start from: {args.init_from}', flush=True)

    mu_H = kanamori_mu_half_filling(args.U, args.J)
    n_modes, n_states, dense_complex_gib = ed_size_summary(args.Nb)

    print('=' * 72, flush=True)
    print(f'Started: {datetime.datetime.now()}', flush=True)
    print(
        f'U={args.U}  J={args.J}  mu_H={mu_H:.8f}  beta={args.beta}  '
        f't1={args.t1}  t2={args.t2}  N_b={args.Nb}  seed={args.seed}',
        flush=True,
    )
    print(
        f'ED size: {n_modes} fermionic modes, {n_states} states, '
        f'one dense complex matrix={dense_complex_gib:.2f} GiB',
        flush=True,
    )
    print(
        f'enforce_sym={enforce_sym}  ph_sym_bath={args.ph_sym_bath}  '
        f'n_iter={args.n_iter}  tol={args.tol}  mix={args.mix}',
        flush=True,
    )
    print(
        f'N_IW={args.N_IW}  n_iw_fit={args.n_iw_fit}  '
        f'conv_n_iw={args.conv_n_iw}  fit_restarts={args.fit_restarts}',
        flush=True,
    )
    print('=' * 72, flush=True)

    t_start = time.time()

    res = run_dmft_2orb(
        U=args.U,
        J=args.J,
        beta=args.beta,
        t1=args.t1,
        t2=args.t2,
        N_b=args.Nb,
        n_iw=args.N_IW,
        n_iter=args.n_iter,
        tol=args.tol,
        mix=args.mix,
        seed=args.seed,
        enforce_sym=enforce_sym,
        ph_sym_bath=args.ph_sym_bath,
        n_iw_fit=args.n_iw_fit,
        conv_n_iw=args.conv_n_iw,
        fit_restarts=args.fit_restarts,
        verbose=True,
        eps_init=eps_init,
        V_init=V_init,
    )

    elapsed = time.time() - t_start

    G_iw = res['G_iw']
    Sig_iw = res['Sig_iw']
    iw_vals = iw_array_from_mesh(G_iw[0].mesh)

    fname = (
        f'{args.outdir}/ED_{branch}_U{args.U:.2f}_J{args.J:.3f}'
        f'_t2{args.t2:.2f}_b{args.beta:.1f}_Nb{args.Nb}.npz'
    )

    np.savez(
        fname,
        iw_vals=iw_vals,
        G1_iw=G_iw[0].data[:, 0, 0],
        G2_iw=G_iw[1].data[:, 0, 0],
        S1_iw=Sig_iw[0].data[:, 0, 0],
        S2_iw=Sig_iw[1].data[:, 0, 0],
        Z1=res['Z1'],
        Z2=res['Z2'],
        D1=res['D1'],
        D2=res['D2'],
        D12=res['D12'],
        n1=res['n1'],
        n2=res['n2'],
        eps_bath=res['eps_bath'],
        V_bath=res['V_bath'],
        converged=res['converged'],
        n_iter=res['n_iter_done'],
        last_delta=res['last_delta'],
        last_fit_err=res['last_fit_err'],
        U=args.U,
        J=args.J,
        mu_H=mu_H,
        beta=args.beta,
        t1=args.t1,
        t2=args.t2,
        N_bath=args.Nb,
        seed=args.seed,
        enforce_sym=enforce_sym,
        ph_sym_bath=args.ph_sym_bath,
        tol=args.tol,
        mix=args.mix,
    )

    print('=' * 72, flush=True)
    print(f'Finished: {datetime.datetime.now()}', flush=True)
    print(
        f'Elapsed: {elapsed/60:.1f} min  '
        f'({res["n_iter_done"]} iters, converged={res["converged"]}, '
        f'last_delta={res["last_delta"]:.4e})',
        flush=True,
    )
    print(f'D1={res["D1"]:.6f}  D2={res["D2"]:.6f}  D12={res["D12"]:.6f}', flush=True)
    print(f'Z1={res["Z1"]:.6f}  Z2={res["Z2"]:.6f}', flush=True)
    print(f'n1={res["n1"]:.4f}  n2={res["n2"]:.4f}', flush=True)
    print('Final eps_bath:', res['eps_bath'], flush=True)
    print('Final V_bath  :', res['V_bath'], flush=True)
    print(f'Saved: {fname}', flush=True)
    print('=' * 72, flush=True)

    if args.fail_on_nonconvergence and not res['converged']:
        sys.exit(1)
    sys.exit(0)


if __name__ == '__main__':
    main()
