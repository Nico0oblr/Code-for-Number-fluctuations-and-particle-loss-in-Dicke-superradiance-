# Particle-number reconstruction

| Figure | Study | Saved input |
| --- | --- | --- |
| Fig. 5(b–e) | Reconstruction versus preparation count, detection efficiency, and mean occupation | [fig5_results.jld2](fig5/data/fig5_results.jld2) |
| Fig. S4 | Sensitivity to smoothing strength and correlation order | [figS4_results.jld2](figS4/data/figS4_results.jld2) |

Each figure directory contains `plot.jl` and `generate.jl`. The JLD2 files
store the study settings, binned observations, intensity kernels, and fits.
The shared Julia code defines the forward model, photon sampling, fitting,
and plotting style.

## Setup

Run from this directory with Julia 1.13.0, the version used to generate
the data and the manifest. The Julia plot scripts use PyPlot, Matplotlib
3.10.8, and LaTeX for axis labels.

```sh
python3 -m venv .venv
. .venv/bin/activate
python -m pip install matplotlib==3.10.8
PYTHON="$(command -v python)" julia --project=. -e 'using Pkg; Pkg.instantiate(); Pkg.build("PyCall")'
```

## Plot the saved data

```sh
julia --project=. fig5/plot.jl
julia --project=. figS4/plot.jl
```

The scripts write `fig5/output/fig5_numeric_source.pdf` and
`figS4/output/figS4.pdf`. Each accepts an optional input JLD2 path and output
PDF path. The Fig. 5 script reproduces the numerical plots used for panels
(b–e), before their placement in the manuscript figure.

## Regenerate data and fits

```sh
julia --project=. --threads=auto fig5/generate.jl
julia --project=. --threads=auto figS4/generate.jl
```

The generators use the settings and seeds at the top of each script. They
save the binned observations before fitting, and write to each figure's
`output/` directory by default. An optional output filename overrides
that location. Existing files are not overwritten.

The acquisition windows recorded in the saved data are `Gamma*T = 5`
for Fig. 5(b–e) and `Gamma*T = 20` for Fig. S4.
