# Saved numerical data

The `.npy` files were selected from Christian Niehof's delivered `NB & CN/data`
directory and copied without changing their contents. `manifest.csv` lists
each supplied array, its SHA-256 hash, and the figure panels that use it.

## Fixed-N intensity kernels

`NumPopDensMatrix-...-In.npy` has shape `(Nmax+1, nt+1)`, with row N
containing the emission kernel for a fully excited fixed-N ensemble.
The matching `...-T.npy` is its time grid. Gamma is set to 1 in these data.
The selected banks have `(Nmax, n) = (100, 20)` and `(1000, 100)`;
`n` sets the default time window `3*log(n)/(Gamma*n)`. The four files
occupy 8,841,344 bytes in total.

## Monte Carlo intensity traces

`MC-...I.npy` stores the ensemble-averaged intensity
`Gamma*<m*(N-m+1)>` at `nt+1` sample times. `MC-...T.npy` is the matching
time grid. Rates, the initial N, and nt appear in the filename. The 5,000
selected files occupy 33,480,000 bytes in total. Individual stochastic
paths and sampling uncertainties are not stored in these files.

The seed construction is given in
[MonteCarlo.py](../python/MonteCarlo.py#L88): a root `SeedSequence(42)`
spawns the per-trajectory seeds. The original arrays do not include
historical seeds or trajectory counts. The requested regeneration settings
and the limitation of these cache keys are described in the
[package README](../README.md#regeneration-and-random-seeds).

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
| `source_file` | Original full-population filename |

The full population histories for N0 = 500 and gr = 0.2, 1, 3 were
regenerated on 2026-11-11 with the retained solver, including ground-state
loss from m = 0. The summaries use the original `get_trajectory`
contractions and per-time intensity sums. The data retain 1,001 sample
times per curve. Full populations have
shape `(tnum+1, (N0+1)*(N0+2)//2)` with flattened index
`N*(N+1)//2 + m`; they are unnecessary for plotting the summaries.

The maximum absolute normalization errors across all 1,001 stored times
are 1.03e-13, 1.12e-15, and 5.31e-14 for gr = 0.2, 1, and 3,
respectively. The total-population contractions agree with the equal-loss
identity `N_mean = N0*exp(-Gg*T)` within 2.85e-11 emitters. Intensity and
mean excited population agree with the earlier histories to numerical
precision; the intensity peak indices are unchanged. No analytical
replacement or error label is applied to the summaries.

See [fig7b_generation.json](fig7b_generation.json) for the solver and
library versions, source-population hashes, and per-history checks.

The summary cache and its regeneration path are implemented in
[`load_trajectory`](../fig7/plot.py#L16).

## Historical population-cache labels

Some older full-population files omit ground-loss gain terms within the
m = 0 sector. Their m = 0 probabilities are incomplete, so their total
probability and mean total emitter number are unreliable. This omission
does not feed back into m > 0: intensity, mean excited population, and
the intensity peak index are unaffected. Do not renormalize the arrays,
because doing so would incorrectly change those observables.

Affected historical files are annotated with a neighboring
`<filename>.metadata.json` sidecar carrying
`data_error = "missing_m0_ground_loss"`. The original array contents and
filenames are preserved. [known_population_errors.json](known_population_errors.json)
records the label, SHA-256 hash, normalization checks, and affected
quantities for seven remaining legacy population files, three superseded
panel (b) population files, and their three archived summaries. In the
archived summaries, the affected saved quantity is `N_mean`.

Match both filename and hash when identifying an affected file: the
regenerated panel (b) populations and summaries reuse the original
filenames. The three currently supplied Fig. 7(b) summaries and their
regenerated full-population sources are normalized and carry no error
label. The registry describes historical data, without applying any
automatic correction to plotting.
