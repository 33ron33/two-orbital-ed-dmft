# Two-Orbital Kanamori ED-DMFT

This repository contains a two-orbital exact-diagonalization dynamical mean-field theory implementation for the Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Main solver

```bash
python solver/cluster_dmft_ed_updated.py \
  --U 4.0 \
  --J 0.0 \
  --beta 25 \
  --t1 1.0 \
  --t2 1.0 \
  --Nb 2 \
  --N_IW 160 \
  --n_iw_fit 70 \
  --conv_n_iw 30 \
  --n_iter 100 \
  --tol 5e-3 \
  --mix 0.10 \
  --fit_restarts 2 \
  --ph_sym_bath \
  --seed metal \
  --outdir data/test
```

## Environment

The working ED-DMFT environment is exported from the original `NN_Solver` conda environment and renamed as:

```bash
ed_pyed
```

Create it with:

```bash
conda env create -f environment.yml
conda activate ed_pyed
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

The ED solver requires TRIQS and pyED. The local `pyed/` folder is included for the exact ED interface used in this project.

## Directory layout

- `solver/`: main ED-DMFT solver.
- `runs/`: scripts used to generate data for individual figures.
- `plotting/`: plotting scripts for PRB-style figures.
- `slurm/`: cluster submission scripts.
- `joblists/`: generated command lists for local or SLURM sweeps.
- `data/`: raw ED-DMFT outputs, ignored by git by default.
- `figures/`: generated figures, ignored by git by default except selected files in `figures/final/`.
- `docs/`: notes, figure manifest, and project documentation.

## Main observables

The solver saves the main observables used throughout the analysis:

- orbital quasiparticle weights `Z1` and `Z2`,
- intra-orbital double occupancies `D1` and `D2`,
- inter-orbital occupancy correlation `D12`,
- orbital fillings `n1` and `n2`,
- Matsubara Green's functions `G1_iw` and `G2_iw`,
- Matsubara self-energies `S1_iw` and `S2_iw`,
- bath parameters `eps_bath` and `V_bath`,
- convergence metadata such as `converged`, `last_delta`, and `n_iter`.

## Figure workflow

The project is organized around the following figure workflow.

### Figure 1: Optional validation

Noninteracting `U=J=0` validation. This checks that the ED-DMFT solver reproduces the noninteracting Bethe Green's function and that the self-energy is numerically zero.

Recommended use: appendix or validation section.

### Figure 2: Degenerate-band Mott transition

Scan A with:

- `J=0`
- `t1=t2=1`
- typically `Nb=2`
- typically `beta=25`

This figure shows the metal-seeded and insulator-seeded DMFT branches through the Mott transition using `Z` and double occupancy.

### Figure 3: Orbital differentiation

Scan C with:

- unequal bandwidths, usually `t2/t1=0.3`
- several Hund ratios `J/U`
- quasiparticle weights `Z1` and `Z2`

This figure shows that the narrow orbital loses coherence before the wide orbital.

### Figure 4: OSMT self-energy

A representative orbital-selective Mott point is selected from Scan C. The figure compares low-frequency self-energies and Green's functions for the wide and narrow orbitals.

### Figure 5: Hund-coupling scan

Scan B at fixed interaction strength. This figure tracks `D1`, `D2`, and `D12` as a function of `J/U`, showing how Hund coupling suppresses intra-orbital charge fluctuations and enhances inter-orbital correlations.

### Figure 6: U-J/U phase diagram

Scan C phase classification using a quasiparticle-weight threshold `Zc`. The phases are classified as:

- both orbitals metallic,
- orbital-selective Mott regime,
- both orbitals localized.

### Figure 7: U-T phase diagram

Scan D phase classification in the `U/t1` versus `T/t1=1/beta` plane at fixed `J/U`.

## Reproducing a small smoke test

After creating the environment, run:

```bash
mkdir -p data/smoke_env_test logs

python solver/cluster_dmft_ed_updated.py \
  --U 0.0 \
  --J 0.0 \
  --beta 25 \
  --t1 1.0 \
  --t2 1.0 \
  --Nb 1 \
  --N_IW 60 \
  --n_iw_fit 20 \
  --conv_n_iw 10 \
  --n_iter 3 \
  --tol 1e-5 \
  --mix 0.5 \
  --seed metal \
  --enforce_sym \
  --outdir data/smoke_env_test \
  2>&1 | tee logs/smoke_env_test.log
```

A successful run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Notes on data

Raw `.npz` output files and logs are ignored by git by default because they can become large. The scripts in `runs/` and `plotting/` are used to regenerate data and figures.

Selected final figures can be placed in:

```bash
figures/final/
```

## License

The original project scripts are released under the license included in this repository.

The local `pyed/` directory contains its own license and copying files. Those terms apply to that directory.
