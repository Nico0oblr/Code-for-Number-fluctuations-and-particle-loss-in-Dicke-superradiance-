# Intensity profiles and particle loss

Numerical code and data for Fig. 4(b-e) and Fig. 7(a-d) of *Number
fluctuations and particle loss in Dicke superradiance*, adapted from
Christian Niehof's `Thesis-DickeVariableN.py` and its numerical modules.

| Figure | Plot script | Saved input |
| --- | --- | --- |
| Fig. 4(b-e) | [fig4/plot.py](fig4/plot.py) | Fixed-N intensity kernels and time grids, `data/NumPopDensMatrix-*.npy` |
| Fig. 7(a,c,d) | [fig7/plot.py](fig7/plot.py) | Monte Carlo ensemble-averaged intensities and time grids, `data/MC-*.npy` |
| Fig. 7(b) | [fig7/plot.py](fig7/plot.py) | Mean-trajectory summaries, `data/fig7b-*.npz` |

## Setup and plotting

Use Python 3.11 or later. The scripts were checked with Python 3.14.7 and
the dependency versions pinned in [requirements.txt](requirements.txt).
From this directory:

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install -r requirements.txt
python fig4/plot.py
python fig7/plot.py
```

The outputs are `fig4/output/fig4_numeric.pdf` and
`fig7/output/fig7_numeric.pdf`. An optional positional argument selects a
different output PDF. Both scripts resolve their input paths relative to
the package and can also run from another working directory. For a server
without a display, set `MPLBACKEND=Agg`.

The Fig. 4 output is a two-by-two grid corresponding to manuscript panels
(b-e). The Fig. 7 output contains panels (a-d) and their shared colorbar.

## Data and numerical conventions

The `data/` directory contains only the inputs requested by these scripts:
four fixed-N kernel/time files, 5,000 Monte Carlo intensity/time files, and
three small trajectory summaries. See [data/README.md](data/README.md) for
array definitions and provenance, and [data/manifest.csv](data/manifest.csv)
for SHA-256 hashes of the supplied arrays.

Historical full-population caches with an incomplete m = 0 sector are
identified in [data/known_population_errors.json](data/known_population_errors.json).
Their m = 0 probabilities and mean total emitter number are unreliable;
intensity and mean excited population are unaffected by this omission.
These labels concern legacy files and archived summaries. The supplied
Fig. 7(b) summaries use regenerated, normalized populations and have no
error label. See the [data notes](data/README.md#historical-population-cache-labels).

The numerical routines retain Christian's emission rates, integration
methods, sampling, distribution weights, and analytical peak formulas.
Unused studies and alternative solvers have been removed. Changes to the
retained code concern portable paths, explicit SciPy imports, command-line
entry points, and the Fig. 7(b) summary cache. Full populations are read by
memory mapping, and the intensity reduction is evaluated in row blocks to
limit temporary memory.

In the code, `Gs` is the collective emission rate Gamma, `Ge` and `Gg` are
the excited- and ground-state loss rates, `N` is the total emitter number,
and `m` is the excited-emitter number. The equal-loss parameter is
`gr = Ge*log(N0)/(Gs*N0)`. The plots use
`t_SR = log(N0)/(Gamma*N0)` and `I_SR = Gamma*N0**2/5`.
The decimal constants in the analytical curves are retained as supplied.

## Regeneration and random seeds

The kernel and Monte Carlo constructors load existing `.npy` files first.
If a required pair is absent, the relevant numerical calculation runs and
saves it under the same name. Move both the intensity/kernel file and its
time file aside to regenerate a dataset. Large-N Monte Carlo calculations
can be expensive.

The Monte Carlo code initializes its random streams with
[`np.random.SeedSequence(42)`](python/MonteCarlo.py#L88), spawns one child
sequence per trajectory, and passes the generated integer seeds to
`np.random.default_rng`. This is the seed prescription for regeneration;
the provided files do not store their historical seeds or trajectory
counts, so those cannot be independently recovered from the arrays.

The retained plotting calls request:

| Panel | Time intervals per trace | Monte Carlo trajectories |
| --- | ---: | ---: |
| Fig. 7(a) | 100 | 10,000 |
| Fig. 7(c) | 100 | 1,000 |
| Fig. 7(d), N = 10^4 and 10^5 | 1,000 | 10,000 |
| Fig. 7(d), N = 10^6 and 10^7 | 100 | 1,000 |

These are the settings requested by the current code, not an assertion of
the settings that produced each historical cache. Cache filenames include
rates, N, and time resolution, but omit seed and trajectory count. Some
panels share a cache at the same rates, N, and resolution. For new runs
with different sampling settings, use a separate copy of the package or
move the relevant cached pair aside. Replotting the supplied data does
not generate random samples.

## Fig. 7(b): small trajectory summaries

The three `.npz` files (91,978 bytes in total) store time, mean total and
excited populations, intensity, the peak index, and parameter metadata.
The trajectories are plotted using the stored values. They replace three
full population histories totaling about 3.02 GB. The full histories are
excluded from version control and are unnecessary for ordinary plotting.

These histories were regenerated with the corrected population generator
for N0 = 500 and gr = 0.2, 1, 3. At all 1,001 stored times, the full
populations sum to 1 within 1.03e-13. Generation settings, source hashes,
and numerical checks are recorded in
[data/fig7b_generation.json](data/fig7b_generation.json). Both mean
trajectories are direct contractions of those populations.

On a summary cache miss, the plot script reduces the corresponding full
population history and saves the summary. To use an external directory
containing the normalized `Rho-...Rho.npy` population files:

```sh
python fig7/plot.py --raw-data-dir /path/to/original/data
```

Existing summaries are always preferred. To rebuild one, first move its
`data/fig7b-N-500-gr-*.npz` file aside, then run the command above.

If neither a summary nor its source is available, the script stops with
an explanation. To explicitly allow regeneration of the full populations:

```sh
python fig7/plot.py --regenerate-trajectories --raw-data-dir /path/to/scratch/data
```

This uses the retained sparse rate-equation solver, writes the large
population files to the selected directory, and saves the small summaries
in this package. Without `--raw-data-dir`, the population files are written
to `data/` and ignored by Git. No random seeds are involved in this solver.
