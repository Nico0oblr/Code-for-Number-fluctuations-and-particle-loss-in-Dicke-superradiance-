get!(ENV, "MPLBACKEND", "Agg")
using PyPlot

function set_figure_font!(fontsize)
    mpl = PyPlot.matplotlib
    mpl.rc("text"; usetex=true)
    mpl.rc("text.latex"; preamble=raw"\usepackage{amsmath}\newcommand{\widebar}[1]{\overline{#1}}")
    mpl.rc("font"; family="serif", size=fontsize)
    mpl.rc("axes"; titlesize=fontsize, labelsize=fontsize)
    mpl.rc("xtick"; labelsize=0.88fontsize)
    mpl.rc("ytick"; labelsize=0.88fontsize)
    mpl.rc("legend"; fontsize=0.84fontsize)
    mpl.rc("figure"; titlesize=fontsize)
    return fontsize
end
