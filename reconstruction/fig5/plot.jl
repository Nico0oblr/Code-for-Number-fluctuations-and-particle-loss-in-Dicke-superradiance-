# Plot only: adapted from reconstruction_plots.ipynb.
using JLD2
include(joinpath(@__DIR__, "..", "plot_style.jl"))

length(ARGS) <= 2 || error("Usage: julia --project=. fig5/plot.jl [input.jld2] [output.pdf]")
data_file = isempty(ARGS) ? joinpath(@__DIR__, "data", "fig5_results.jld2") : abspath(ARGS[1])
output = length(ARGS) < 2 ? joinpath(@__DIR__, "output", "fig5_numeric_source.pdf") : abspath(ARGS[2])
mkpath(dirname(output))

saved = JLD2.load(data_file)
haskey(saved, "fits") || error("Observations are saved, but fits have not completed.")
observations = saved["observations"]
fits = saved["fits"]
settings = fits.settings
(; mu, nmax, T, dt, tau, eta, b, lambda, solver_tol, n_runs,
   R_fixed, R_values, eta_values, mean_values, families) = settings
cases = observations.cases
summary_results = fits.summary_results
distribution_results = fits.distribution_results
nothing

# APS/PRA RevTeX 4.2: a4paper, twocolumn, 10pt -> columnwidth = 246 TeX pt.
# Include the exported figure at width=\columnwidth (one column for all four panels).
paper_columnwidth_pt = 246.0
paper_canvas_width_in = 7.0
paper_reduction = (paper_columnwidth_pt / 72.27) / paper_canvas_width_in
paper_label_font = 10.0 / paper_reduction
paper_tick_font = 9.0 / paper_reduction
paper_legend_font = 9.0 / paper_reduction
paper_schematic_font = 8.0 / paper_reduction
paper_line_width = 1.3 / paper_reduction
paper_reference_width = 1.0 / paper_reduction
paper_errorbar_width = 0.8 / paper_reduction
paper_frame_width = 0.8 / paper_reduction
paper_marker_size = 4.0 / paper_reduction
paper_marker_edge = 0.7 / paper_reduction
set_figure_font!(paper_label_font)
PyPlot.matplotlib.rc("xtick"; labelsize=paper_tick_font)
PyPlot.matplotlib.rc("ytick"; labelsize=paper_tick_font)
fig_uncertainty, panel_axes = subplots(2, 2, figsize=(paper_canvas_width_in, paper_canvas_width_in))
fig_uncertainty.subplots_adjust(left=0.12, right=0.98, bottom=0.19,
    top=0.88, wspace=0.40, hspace=0.46)
colors = Dict("coherent" => "#0072B2", "thermal" => "#D55E00")

ax = panel_axes[1, 1]
p_peak = maximum(maximum(cases[(f, mu)].p_true) for f in families)
for family in families
    case = cases[(family, mu)]
    result = distribution_results[family]
    color = colors[family]
    ax.fill_between(case.n, result.lower, result.upper,
        step="mid", color=color, alpha=0.18, linewidth=0)
    ax.step(case.n, case.p_true, where="mid", color=color, linewidth=paper_line_width)
    ax.plot(case.n, result.average, "o", color=color, markersize=paper_marker_size,
        markerfacecolor="white", markeredgewidth=paper_marker_edge)
end
ax.set_xlim(-0.5, 15.5)
ax.set_xticks([0, 5, 10, 15])
ax.set_ylim(0, 1.7p_peak)
ax.set_yticks([0.0, 0.2])
ax.set_xlabel(raw"$n$")
ax.set_ylabel(raw"$p_n$")
node_positions = (0.08, 0.36, 0.64, 0.92)
node_labels = (raw"$p_n$", raw"$I(t)$", "data", raw"$p_n^{\mathrm{rec}}$")
for (x, label) in zip(node_positions, node_labels)
    ax.text(x, 0.88, label, transform=ax.transAxes,
        ha="center", va="center", fontsize=paper_schematic_font,
        bbox=Dict("boxstyle" => "round,pad=0.25", "fc" => "#F5F5F5", "ec" => "#BBBBBB"))
end
for i in 1:3
    ax.annotate("", xy=(node_positions[i+1]-0.09, 0.88),
        xytext=(node_positions[i]+0.09, 0.88), xycoords="axes fraction",
        arrowprops=Dict("arrowstyle" => "->", "color" => "#555555", "lw" => paper_frame_width))
end

panels = (
    (ax=panel_axes[1,2], scan=:R, x=R_values, scale="log", xlabel=raw"$R$"),
    (ax=panel_axes[2,1], scan=:eta, x=eta_values, scale="linear", xlabel=raw"$\eta$"),
    (ax=panel_axes[2,2], scan=:N, x=mean_values, scale="log", xlabel=raw"$\langle \hat{N}\rangle$"),
)
for panel in panels
    local ax = panel.ax
    for family in families
        result = summary_results[(panel.scan, family, :g2)]
        ax.errorbar(panel.x, result.center, yerr=result.spread, fmt="o-",
            color=colors[family], linewidth=paper_line_width, markersize=paper_marker_size,
            markerfacecolor="white", markeredgewidth=paper_marker_edge,
            elinewidth=paper_errorbar_width, capsize=1.5/paper_reduction, capthick=paper_errorbar_width)
        ax.axhline(cases[(family, mu)].truth.g2, color=colors[family],
            linestyle="--", linewidth=paper_reference_width, alpha=0.8)
    end
    ax.set_xscale(panel.scale)
    # Include complete error bars; no hard upper limit concealing large scatter.
    lower = minimum(minimum(summary_results[(panel.scan, f, :g2)].lower) for f in families)
    upper = maximum(maximum(summary_results[(panel.scan, f, :g2)].upper) for f in families)
    ax.set_ylim(min(0.8, lower-0.1), max(2.2, upper+0.1))
    ax.set_xlabel(panel.xlabel)
    ax.set_ylabel(raw"$g_N^{(2)}$")
end
panel_axes[2,1].set_xlim(0, 1.02)
panel_axes[1,2].xaxis.set_major_locator(PyPlot.matplotlib.ticker.LogLocator(base=10, numticks=5))
panel_axes[2,2].set_xticks([1, 3, 10, 20])
panel_axes[2,2].set_xticklabels(["1", "3", "10", "20"])
panel_axes[2,2].set_yticks([1.0, 1.5, 2.0])
panel_axes[1,2].set_yticks([1.0, 1.5, 2.0])
panel_axes[2,1].set_yticks([1.0, 1.5, 2.0])


for (ax, letter) in zip((panel_axes[1,1], panel_axes[1,2], panel_axes[2,1], panel_axes[2,2]), ("a", "b", "c", "d"))
    ax.text(-0.18, 0.9, "($letter)", transform=ax.transAxes, fontsize=paper_label_font, va="bottom")
    for spine in ("left", "right", "top", "bottom")
        ax.spines[spine].set_visible(true)
        ax.spines[spine].set_linewidth(paper_frame_width)
    end
    ax.tick_params(direction="out", length=2.5/paper_reduction, width=paper_frame_width)
    ax.grid(false, which="both")
end

Line2D = PyPlot.matplotlib.lines.Line2D
state_handles = [Line2D([], [], color=colors[f], linewidth=paper_line_width) for f in families]
fig_uncertainty.legend(state_handles, ["Coherent", "Thermal"],
    loc="upper center", bbox_to_anchor=(0.5, 0.95), ncol=2,
    frameon=false, fontsize=paper_legend_font, handlelength=2, columnspacing=2)
style_handles = [
    Line2D([], [], color="black", linewidth=paper_line_width),
    Line2D([], [], color="black", linestyle="--", linewidth=paper_reference_width),
    Line2D([], [], color="black", linestyle="-", marker="o", markersize=paper_marker_size,
        markerfacecolor="white", markeredgewidth=paper_marker_edge),
]
fig_uncertainty.legend(style_handles,
    [raw"True $p_n$", raw"True $g_N^{(2)}$", "Reconstructed"],
    loc="lower center", bbox_to_anchor=(0.5, 0.005), ncol=3,
    frameon=false, fontsize=paper_legend_font, handlelength=1.8, columnspacing=1.3)
nothing

fig_uncertainty.savefig(output; bbox_inches=fig_uncertainty.bbox_inches, pad_inches=0)
println("Saved $output")
