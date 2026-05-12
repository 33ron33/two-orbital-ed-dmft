# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.# Two-Orbital Kanamori ED-DMFT

This repository contains an exact-diagonalization dynamical mean-field theory implementation for the two-orbital Hubbard-Kanamori model on the Bethe lattice.

The project focuses on:

- degenerate-band Mott transition,
- Hund-coupling effects,
- orbital differentiation,
- orbital-selective Mott physics,
- U-J/U and U-T phase diagrams.

## Repository layout

- `solver/`: ED-DMFT solver.
- `plotting/`: plotting and diagnostic workflow.
- `runs/`: scripts used to generate parameter scans.
- `pyed/`: local pyED dependency used by the solver.
- `figures/final/`: selected final figures.
- `environment.yml`: conda environment for running the solver and plots.
- `requirements_plotting.txt`: lightweight plotting-only requirements.

Raw `.npz` data, logs, generated job lists, and temporary figures are ignored by git.

## Environment

Create the conda environment with:

```bash
conda env create -f environment.yml
conda activate ed_dmft
pip install -e ./pyed
```

For plotting-only usage:

```bash
pip install -r requirements_plotting.txt
```

## Solver example

```bash
python solver/ed_solver_dmft.py \
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

## Small smoke test

```bash
mkdir -p data/smoke_env_test logs

python solver/ed_solver_dmft.py \
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

A successful noninteracting run should give approximately:

- `Z1 = Z2 = 1`
- `D1 = D2 = 0.25`
- `n1 = n2 = 1`
- `converged=True`

## Plotting workflow

The main plotting script is:

```bash
python plotting/ed_dmft_figures.py --help
```

Example smoke plotting command:

```bash
python plotting/ed_dmft_figures.py \
  --mode smoke \
  --data data/smoke_env_test \
  --out figures/test_plot_workflow \
  --nb 1 \
  --t2 1.0
```

Example production-style command:

```bash
python plotting/ed_dmft_figures.py \
  --mode production \
  --scanA data/scanA \
  --scanB data/scanB \
  --scanC data/scanC \
  --scanD data/scanD \
  --out figures/production \
  --nb 2 \
  --beta 25 \
  --t2A 1.0 \
  --t2B 1.0 \
  --t2C 0.3 \
  --t2D 0.3 \
  --z_thresh 0.15
```

## Main observables

The solver writes NumPy archives containing:

- quasiparticle weights `Z1`, `Z2`,
- double occupancies `D1`, `D2`,
- inter-orbital correlation `D12`,
- fillings `n1`, `n2`,
- Green's functions `G1_iw`, `G2_iw`,
- self-energies `S1_iw`, `S2_iw`,
- final bath parameters `eps_bath`, `V_bath`,
- convergence metadata.

## License

This repository is distributed under the GNU General Public License v3.0.

The included `pyed/` directory carries its own GPL license files. Those terms apply to that directory.
