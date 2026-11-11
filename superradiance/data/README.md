# Saved numerical data

This folder contains the data used by the Python plotting scripts.
[manifest.csv](manifest.csv) lists the arrays, their checksums, and the
figure panels that use them.

## Fixed-N intensity kernels

`NumPopDensMatrix-...-In.npy` has shape `(Nmax+1, nt+1)`, with row N
containing the emission kernel for a fully excited fixed-N ensemble.
The matching `...-T.npy` is its time grid. Gamma is set to 1 in these data.
The selected banks have `(Nmax, n) = (100, 20)` and `(1000, 100)`;
`n` sets the default time window `3*log(n)/(Gamma*n)`.

## Monte Carlo intensity traces

`MC-...I.npy` stores the ensemble-averaged intensity
`Gamma*<m*(N-m+1)>` at `nt+1` sample times. `MC-...T.npy` is the matching
time grid. Rates, the initial N, and nt appear in the filename. These
files contain the averages, without individual trajectories or sampling
uncertainties.

See the [package README](../README.md#regeneration-and-random-seeds) for
simulation settings and random seeds.

## Fig. 7(b) trajectory summaries

`fig7b-N-500-gr-{0.2,1,3}.npz` contains:

| Key | Meaning |
| --- | --- |
| `T` | Physical sample times with Gamma = 1 |
| `N_mean` | Mean total emitter number at each time |
| `m_mean` | Mean excited-emitter number at each time |
| `I` | Physical intensity `Gamma*<m*(N-m+1)>` |
| `peak_index` | Index of the largest sampled intensity |
| `N0`, `Gs`, `Ge`, `Gg`, `gr`, `tnum`, `Tmax` | Initial population, rates, scaled loss, and time-grid settings |
| `source_file` | Full-population filename used to calculate the summary |

Each file contains 1,001 sample times. The plotting script loads these
summaries directly, so it does not need the full Dicke-triangle
populations. The solver settings are in
[fig7b_generation.json](fig7b_generation.json).

## Older population files

Some older full-population files have an incomplete m = 0 sector. Their
m = 0 probabilities and mean total emitter number are unreliable.
Intensity, mean excited population, and the intensity peak index are
unaffected by this omission. Do not renormalize these arrays, as that
would incorrectly change those observables.

Affected files have a `<filename>.metadata.json` file beside them with
the label `data_error = "missing_m0_ground_loss"`. They are also listed in
[known_population_errors.json](known_population_errors.json). Check both
the filename and checksum when using older data, since different versions
can share a filename.
