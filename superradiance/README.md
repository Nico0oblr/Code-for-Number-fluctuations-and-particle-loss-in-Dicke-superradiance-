# Intensity profiles and particle loss

Python code and data for Fig. 4(b-e) and Fig. 7(a-d) of *Number
fluctuations and particle loss in Dicke superradiance*. Figure 4 compares
emission profiles for different initial particle-number distributions.
Figure 7 shows the effect of atom loss on the emission and ensemble
dynamics.

| Figure | Plot script | Saved input |
| --- | --- | --- |
| Fig. 4(b-e) | [fig4/plot.py](fig4/plot.py) | Fixed-N intensity kernels and time grids, `data/NumPopDensMatrix-*.npy` |
| Fig. 7(a,c,d) | [fig7/plot.py](fig7/plot.py) | Monte Carlo ensemble-averaged intensities and time grids, `data/MC-*.npy` |
| Fig. 7(b) | [fig7/plot.py](fig7/plot.py) | Mean-trajectory summaries, `data/fig7b-*.npz` |

## Setup and plotting

Use Python 3.11 or later with the dependencies in
[requirements.txt](requirements.txt). Run the following from this directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python fig4/plot.py
python fig7/plot.py
```

The scripts save `fig4/output/fig4_numeric.pdf` and
`fig7/output/fig7_numeric.pdf`. You can pass a different output filename
as an argument. Both scripts also work from another working directory.
For a server without a display, set `MPLBACKEND=Agg`.

The Fig. 4 output is a two-by-two grid corresponding to manuscript panels
(b-e). The Fig. 7 output contains panels (a-d) and their shared colorbar.

## Data and numerical conventions

The `data/` folder contains fixed-N emission kernels, averaged Monte Carlo
intensity traces, and mean trajectories. The [data README](data/README.md)
describes the arrays and a known issue in some older population files.
[manifest.csv](data/manifest.csv) lists checksums for the saved arrays.

In the code, `Gs` is the collective emission rate Gamma, `Ge` and `Gg` are
the excited- and ground-state loss rates, `N` is the total emitter number,
and `m` is the excited-emitter number. The equal-loss parameter is
`gr = Ge*log(N0)/(Gs*N0)`. The plots use
`t_SR = log(N0)/(Gamma*N0)` and `I_SR = Gamma*N0**2/5`.

## Regeneration and random seeds

The scripts use saved `.npy` files when available. If a required dataset
is missing, they run the simulation and save the result. To rerun a
dataset, move both its intensity or kernel file and its time file aside.
Monte Carlo simulations at large N can take a long time.

The Monte Carlo generator starts from
[`np.random.SeedSequence(42)`](python/MonteCarlo.py#L88) and derives a
separate random seed for each trajectory. The original saved traces do
not record their seeds or trajectory counts. The current simulation
settings are:

| Panel | Time intervals per trace | Monte Carlo trajectories |
| --- | ---: | ---: |
| Fig. 7(a) | 100 | 10,000 |
| Fig. 7(c) | 100 | 1,000 |
| Fig. 7(d), N = 10^4 and 10^5 | 1,000 | 10,000 |
| Fig. 7(d), N = 10^6 and 10^7 | 100 | 1,000 |

Filenames record the rates, N, and time resolution, but not the seed or
trajectory count. Move existing results aside before changing those
settings. Plotting the saved traces does not generate new random samples.

## Fig. 7(b): trajectory data

The `.npz` files contain the mean total and excited populations, intensity,
and time grid for N0 = 500 and gr = 0.2, 1, 3. These are all the data
needed to plot panel (b); the full population histories are not included.

If a summary is missing, the script can calculate it from a full
population history. To read `Rho-...Rho.npy` files from another directory:

```sh
python fig7/plot.py --raw-data-dir /path/to/original/data
```

Existing summaries are loaded first. To rebuild one, move its
`data/fig7b-N-500-gr-*.npz` file aside before running the command above.

If you have neither the summary nor the full population history, run:

```sh
python fig7/plot.py --regenerate-trajectories --raw-data-dir /path/to/scratch/data
```

This solves the population rate equations, writes the full histories to
the chosen directory, and saves the summaries in `data/`. Without
`--raw-data-dir`, the full histories are also written to `data/` and are
ignored by Git. This calculation does not use random sampling.
