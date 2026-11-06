# Plot only: adapted from smoothing_lambda_plots.ipynb.
using JLD2
include(joinpath(@__DIR__, "..", "plot_style.jl"))

length(ARGS) <= 2 || error("Usage: julia --project=. figS4/plot.jl [input.jld2] [output.pdf]")
data_file = isempty(ARGS) ? joinpath(@__DIR__, "data", "figS4_results.jld2") : abspath(ARGS[1])
output = length(ARGS) < 2 ? joinpath(@__DIR__, "output", "figS4.pdf") : abspath(ARGS[2])
mkpath(dirname(output))

saved = JLD2.load(data_file)
haskey(saved, "fits") || error("Observations are saved, but fits have not completed.")
observations = saved["observations"]
fits = saved["fits"]
settings = fits.settings
(; mu, nmax, T, dt, tau, eta, b, R, n_runs, solver_tol,
   lambda_values, families, moment_orders) = settings
records = observations.records
results = fits.results
nothing

# Display choices only; these select summaries already calculated by the script.
moment_lambda = settings.moment_lambda
distribution_lambdas = settings.distribution_lambdas
distribution_n_display = settings.distribution_n_display

moment_results = Dict{String,Any}()
distribution_results = Dict{String,Any}()
for family in families
    scan = results[family]
    index = findfirst(r -> isapprox(r.lambda, moment_lambda; rtol=1e-10, atol=0), scan)
    isnothing(index) && error("moment_lambda was not saved in this dataset.")
    moment_results[family] = scan[index].moments
    selected = Any[]
    for smoothing in distribution_lambdas
        index = findfirst(r -> isapprox(r.lambda, smoothing; rtol=1e-10, atol=0), scan)
        isnothing(index) && error("lambda=$smoothing was not saved in this dataset.")
        push!(selected, scan[index].distribution)
    end
    distribution_results[family] = selected
end
nothing

set_figure_font!(18)

fig_combined, axes_combined = subplots(3, 4, figsize=(17, 10.5))

panel_colors = Dict(
    "coherent" => "#0072B2", "thermal" => "#D55E00",
    "fock" => "#009E73", "even-cat" => "#CC79A7",
)
panel_titles = Dict(
    "coherent" => "Coherent", "thermal" => "Thermal",
    "fock" => "Fock", "even-cat" => "Even cat",
)
fit_colors = ("#0072B2", "#D55E00", "#009E73", "#CC79A7")
fit_markers = ("o", "s", "^", "D")

for (column, family) in enumerate(families)
    color = panel_colors[family]

    # First row: smoothing scan.
    ax = axes_combined[1, column]
    scan = results[family]
    x = [result.lambda for result in scan]
    average = [result.average for result in scan]
    spread = [result.spread for result in scan]
    clean = [result.noiseless for result in scan]

    ax.fill_between(x, average .- spread, average .+ spread,
        color=color, alpha=0.18, linewidth=0)
    ax.plot(x, average, "o-", color=color,
        linewidth=1.6, markersize=6,
        markerfacecolor="white", markeredgewidth=1.2,
        markevery=3, label=raw"Reconstructed: mean $\pm$ SD")
    ax.plot(x, clean, ":", color="black",
        linewidth=1.7, label="Noiseless fit")
    ax.axhline(records[family].truth, color="black",
        linestyle="--", linewidth=1.2, label="True")

    ax.set_xscale("symlog", linthresh=1e-7, linscale=0.6)
    ax.set_xlim(0, maximum(x))
    ax.set_xticks([0.0, 1e-6, 1e-4, 1e-2])
    ax.set_xlabel(raw"$\lambda$")
    ax.set_title(panel_titles[family], pad=12)

    if family == "thermal"
        ax.set_yticks([1.98, 2.0])
    elseif family == "fock"
        ax.set_yticks([0.8, 0.9])
    end

    # Second row: reconstructed moments.
    ax = axes_combined[2, column]
    result = moment_results[family]

    ax.errorbar(moment_orders, result.average,
        yerr=result.spread, fmt="o", color=color,
        markersize=7, markerfacecolor="white",
        markeredgewidth=1.3, elinewidth=1.2,
        capsize=4, capthick=1.2,
        label=raw"Reconstructed: mean $\pm$ SD")

    for (i, k) in enumerate(moment_orders)
        ax.hlines(result.truth[i], k-0.22, k+0.22,
            color="black", linestyle="--", linewidth=1.6,
            label=i == 1 ? "True" : nothing)
    end

    ax.set_xticks(moment_orders)
    ax.set_xlim(1.5, 4.5)
    ax.set_xlabel(raw"$k$")

    if family in ("coherent", "even-cat")
        ax.set_yticks([0.0, 2.5, 5.0])
    end

    # Third row: full distributions at the selected smoothing strengths.
    ax = axes_combined[3, column]
    numbers = collect(0:length(records[family].p_true)-1)
    n_display = min(distribution_n_display, last(numbers))

    ax.step(numbers, records[family].p_true,
        where="mid", color="black",
        linewidth=1.8, linestyle="-", label="True")

    for (i, result) in enumerate(distribution_results[family])
        style = mod1(i, length(fit_colors))
        ax.errorbar(numbers, result.average,
            yerr=result.spread,
            fmt=fit_markers[style], color=fit_colors[style],
            markersize=5.5, markerfacecolor="white",
            markeredgewidth=1.2, elinewidth=1.0,
            capsize=2, capthick=1.0,
            label=string(raw"$\lambda=", result.lambda, raw"$"))
    end

    ax.set_xlim(-0.5, n_display + 0.5)
    ax.set_xticks(collect(0:5:n_display))
    ax.set_xlabel(raw"$n$")
    ax.set_yscale("log")
    ax.set_ylim(1e-6, 2e0)
end

# One axis label per row; vertical scales remain independent.
axes_combined[1, 1].set_ylabel(raw"$g_N^{(2)}$")
axes_combined[2, 1].set_ylabel(raw"$g_N^{(k)}$")
axes_combined[3, 1].set_ylabel(raw"$P_n$")

# Shared legend for the first two rows.
handles, labels = axes_combined[1, 1].get_legend_handles_labels()
fig_combined.legend(handles, labels,
    loc="upper center", bbox_to_anchor=(0.5, 1.0),
    ncol=3, frameon=false, fontsize=14)

# Smoothing-strength legend for the distribution row.
handles, labels = axes_combined[3, 1].get_legend_handles_labels()
fig_combined.legend(handles, labels,
    loc="lower center", bbox_to_anchor=(0.5, 0.0),
    ncol=length(handles), frameon=false, fontsize=14)

# Consistent panel labels: (a)–(d), (e)–(h), (i)–(l).
for row in 1:3, column in 1:4
    ax = axes_combined[row, column]
    letter = Char(Int('a') + 4 * (row - 1) + column - 1)
    ax.text(-0.18, 0.9, "($letter)",
        transform=ax.transAxes, fontsize=18, va="bottom")

    for spine in ("left", "right", "top", "bottom")
        ax.spines[spine].set_visible(true)
    end
    ax.tick_params(direction="out", length=4, width=0.9)
    ax.grid(false, which="both")
end

for column in 1:4
    axes_combined[3, column].set_yticks([1e-6, 1e-4, 1e-2])
end

fig_combined.tight_layout(
    rect=(0, 0.05, 1, 0.95),
    h_pad=1.5, w_pad=1.5,
)

fig_combined.savefig(output; bbox_inches="tight")
println("Saved $output")
