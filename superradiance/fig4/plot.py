"""Numerical plots adapted from Christian Niehof's Thesis-DickeVariableN.py."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np

from python.NumPopDensMatrix import NumPopDensMatrix

def pick_rho(n):
    if n == 20:
        return NumPopDensMatrix({"nmax": 100, "Gs": 1, "nt": 1000, "n": 20})
    if n == 100:
        return NumPopDensMatrix({"nmax": 1000, "Gs": 1, "nt": 1000, "n": 100})
    raise ValueError("Fig. 4 uses mean occupations 20 and 100.")

def fig4(output):
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "mathtext.fontset": "cm",
    })
    # 2x2 Raster
    fig, axes = plt.subplots(
        2, 2,
        figsize=(246/72, 246/72 * 3/4),
        constrained_layout=True
    )
    cbars=[]

    # (a)
    ax = axes[0, 0]
    n = 20
    rho = pick_rho(n)
    T = rho.T * n / np.log(n)

    ax.plot(T, rho.In[n,:]/(n**2/5), color='black')
    ax.plot(T, [(n+(n**2-2*n)*t)/(n**2/5) for t in rho.T], color='black', ls='--')
    ax.set_xlim(0, 3)
    ax.set_xticks([0, 1, 2, 3])
    ax.tick_params(axis='x', labelbottom=False)
    ax.set_yticks([0, 1])
    ax.set_ylabel(r"$I(t)/I_\mathrm{SR}$")
    ax.set_ylim(0, 1.1)

    # (b)
    ax = axes[0,1]
    n = 20
    rho = pick_rho(n)
    T = rho.T * n / np.log(n)

    sn=[0.5, 1, 1.5, 2, 3, 5]

    cmap = plt.cm.turbo
    norm = mpl.colors.Normalize(vmin=0, vmax=max(sn))
    ax.plot(T, rho.In[n,:]/(n**2/5), color=cmap(norm(0)))
    for s in sn:
        ax.plot(T, rho.intensity({"statistics": "Gaussian", "n": n, "sigma": s*n})/(n**2/5), color=cmap(norm(s)))
    ax.set_xticks([0, 1, 2, 3])
    ax.tick_params(axis='x', labelbottom=False)
    ax.set_yticks([0, 1, 2, 3])

    # ax.set_xlabel(r"$t/t_\mathrm{SR}$")
    # ax.set_ylabel(r"$I(t)/I_\mathrm{SR}$")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_ticks([0, 1, 2, 3, 4, 5])
    cbar.set_label("$\\sigma/\\langle\\hat N\\rangle$")
    cbars.append(cbar)

    # (c)
    ax = axes[1,0]
    n = 20
    rho = pick_rho(n)
    T = rho.T * n / np.log(n)

    ax.plot(T, rho.intensity({"statistics": "Coherent", "n": n})/(n**2/5), color='red')
    ax.plot(T, rho.intensity({"statistics": "Cat", "n": n, "normalized": True})/(n**2/5), color='green', ls=(0,(5,5)))
    ax.plot(T, rho.intensity({"statistics": "Squeezed", "n": n})/(n**2/5), color='orange')
    ax.plot(T, rho.intensity({"statistics": "Thermal", "n": n}) / (n ** 2 / 5), color='blue', ls='--')

    ax.set_xlabel(r"$t/t_\mathrm{SR}$")
    ax.set_ylabel(r"$I(t)/I_\mathrm{SR}$")

    # (d)
    ax = axes[1,1]
    n = 100
    rho = pick_rho(n)
    T = rho.T * n / np.log(n)

    G = [0.2, 0.5, 1, 2]

    cmap = plt.cm.turbo
    norm = mpl.colors.Normalize(vmin=0, vmax=max(G))
    ax.plot(T, rho.In[n,:]/(n**2/5), color=cmap(norm(0)))
    for g in G:
        ax.plot(T, rho.intensity({"statistics": "Lorentz", "n": n, "gamma": g})/(n**2/5), color=cmap(norm(g)))
    ax.set_yticks([0, 1, 2])
    ax.set_xlabel(r"$t/t_\mathrm{SR}$")

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax)
    cbar.set_ticks([0, 1, 2])
    cbar.set_label("$\\gamma/\\langle\\hat N\\rangle$")
    cbars.append(cbar)

    for ax in axes.flat:
        ax.set_xlim(0,3)
        ax.set_xticks([0,1,2,3])

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output)
    plt.close(fig)
    print(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot the numerical panels (b-e) of Fig. 4.")
    parser.add_argument("output", nargs="?", type=Path,
                        default=Path(__file__).resolve().parent / "output/fig4_numeric.pdf")
    args = parser.parse_args()
    fig4(args.output)
