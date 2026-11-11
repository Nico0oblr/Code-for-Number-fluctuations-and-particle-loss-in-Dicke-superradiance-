"""Numerical plots adapted from Christian Niehof's Thesis-DickeVariableN.py."""
import argparse
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
from concurrent.futures import ProcessPoolExecutor as PPE
from tqdm import tqdm
from python.MonteCarlo import MonteCarlo as MC

def load_trajectory(n, gr, raw_data_dir=None, regenerate=False):
    """Prefer the small summary; reduce the full populations on a cache miss."""
    summary_file = ROOT / "data" / f"fig7b-N-{n}-gr-{gr}.npz"
    if summary_file.is_file():
        with np.load(summary_file, allow_pickle=False) as saved:
            return {name: saved[name] for name in saved.files}

    from python.NumPopDensLosses import NumPopDensLosses
    data_dir = ROOT / "data" if raw_data_dir is None else Path(raw_data_dir)
    rate = gr * n / np.log(n)
    tmax = 3 * np.log(n) / n
    source_file = (f"Rho-Gs-1-Gg-{rate}-Ge-{rate}-tnum-1000-n-{n}"
                   f"-Tmax-{tmax}Rho.npy")
    if not (data_dir / source_file).is_file() and not regenerate:
        raise FileNotFoundError(
            f"Missing {summary_file.name}. Supply the original populations with "
            "--raw-data-dir PATH, or use --regenerate-trajectories to run the loss solver.")
    rho = NumPopDensLosses({"n": n, "Gs": 1, "Gg": rate, "tnum": 1000,
                           "data_dir": data_dir})
    N_mean, m_mean = rho.get_trajectory()
    weights = np.concatenate([np.arange(k + 1) * (k - np.arange(k + 1) + 1)
                              for k in range(n + 1)])
    # Same per-row sum as the original intensity(), with bounded temporary memory.
    intensity = np.empty(len(rho.T))
    for start in range(0, len(rho.T), 32):
        stop = start + 32
        intensity[start:stop] = np.sum(rho.rho[start:stop] * weights, axis=1)
    saved = dict(T=rho.T, N_mean=N_mean, m_mean=m_mean, I=intensity,
                 peak_index=np.argmax(intensity), N0=n, Gs=1, Ge=rate, Gg=rate,
                 gr=gr, tnum=1000, Tmax=tmax, source_file=source_file)
    np.savez_compressed(summary_file, **saved)
    return saved

def calculate_I(params):
    Gs = 1
    if len(params) == 2:
        N, gr = params
        nt = 10**3 if N==10**4 or N==10**5 else 10**2
        nMC = 10**4 if N==10**4 or N==10**5 else 10**3
    elif len(params) == 3:
        N, gr, b = params
        nt = 10**2
        nMC = 10**3
    return MC({"Gs": Gs, "Ge": gr * Gs * N / np.log(N), "nt": nt, "N": N, "nMC": nMC, "sequential": True}).getIMax()

def plot_I_gr(GR, I, N, cl, fig):
    GR = np.asarray(GR)
    I = np.asarray(I)
    # gr=1
    fig.semilogy(1, 0.3891/np.log(N)**2 *5, marker='*', color=cl)

    # gr>1
    mask = GR >= 1
    imax = [4*gr**2/np.e**2 * N**(1/gr-1)/np.log(N)**2 *5 for gr in GR[mask]]
    fig.semilogy(GR[mask], imax, ls=':', color=cl)


    fig.semilogy(GR, I, color=cl)

def fig7(output, raw_data_dir=None, regenerate_trajectories=False):
    plt.rcParams.update({
        "font.size": 10,
        "axes.labelsize": 10,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "mathtext.fontset": "cm",
    })
    fig = plt.figure(
        figsize=(246 / 72, 246 / 72 * 3 / 4 * 1.1),
        constrained_layout=True
    )
    gs = fig.add_gridspec(3, 2, height_ratios=[1, 1, 0.05])
    axes = [
        fig.add_subplot(gs[0, 0]),
        fig.add_subplot(gs[0, 1]),
        fig.add_subplot(gs[1, 0]),
        fig.add_subplot(gs[1, 1]),
        fig.add_subplot(gs[2, :])
    ]

    cmap = plt.cm.turbo
    norm = mpl.colors.Normalize(vmin=0, vmax=5)

    # (a)
    ax = axes[0]
    N0 = int(10 ** 7)
    GE = [2914097.1729438975, 2538084.6344995233, 2068068.9614440564, 1598053.288388589, 1128037.6153331215,
          658021.9422776542, 407346.9166480717, 0.0]

    for ge in GE:
        mc = MC({"Gs": 1, "Ge": ge, "nt": 100, "N": N0, "sequential": False})
        ax.semilogy(mc.getT(), mc.getI(), color=cmap(norm(ge*np.log(N0)/N0)))

    ax.set_xlabel(r"$t/t_\mathrm{SR}$")
    ax.set_xticks([0, 1, 2, 3])

    ax.set_yticks([1, 10 ** -7, 10 ** -14], labels=["$N^0$", "$N^{-1}$", "$N^{-2}$"])
    ax.set_ylabel(r"$I(t)/I_\mathrm{SR}$")

    # (b)
    ax = axes[1]
    n = 500
    GR = [0.2, 1, 3]
    ax.fill([0, n, n], [0, 0, n], facecolor="lightgray", edgecolor="black")
    for k, gr in enumerate(tqdm(GR)):
        saved = load_trajectory(n, gr, raw_data_dir, regenerate_trajectories)
        tmax = int(saved["peak_index"])
        trajectory = saved["N_mean"], saved["m_mean"]
        ax.plot(*trajectory, color=cmap(norm(gr)))
        ax.plot(*tuple(x[tmax] for x in trajectory), "+", color=cmap(norm(gr)), markersize=10)
    ax.set_xticks([0, n], labels=["$0$", "$N$"])
    ax.set_xlabel("$\\langle\\hat N\\rangle$")
    ax.set_yticks([0, n], labels=["$0$", "$N$"])
    ax.set_ylabel("$\\langle\\hat m\\rangle$")

    # (c)
    ax = axes[2]
    with PPE() as executor:
        GR = [0.2, 1., 3.]
        N = [int(10 ** logn) for logn in np.linspace(4, 7, 100)]
        for gr in GR:
            params = [(n, gr, True) for n in N]
            I = list(tqdm(executor.map(calculate_I, params),
                          total=len(params),
                          desc=f"gr = {gr}"
                          ))
            ax.loglog(N, I, color=cmap(norm(gr)))

            if gr < 1:
                I = [0.1957 * (1 - gr) ** 2 * 5 for n in N]
                ax.loglog(N, I, color=cmap(norm(gr)), linestyle='--')
            elif gr > 1:
                I = [4 * gr ** 2 / np.e ** 2 * n ** (1 / gr - 1) / np.log(n) ** 2 * 5 for n in N]
                ax.loglog(N, I, color=cmap(norm(gr)), linestyle='--')
            else:
                I = [0.3891 / np.log(n) ** 2 * 5 for n in N]
                ax.loglog(N, I, color=cmap(norm(gr)), linestyle='--')

    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$I_\mathrm{max}/I_\mathrm{SR}$")

    # (d)
    ax = axes[3]
    with PPE() as executor:
        for gr in [0.2, 1, 3]:
            ax.axvline(x=gr, color=cmap(norm(gr)), alpha=0.7)

        GR = np.linspace(0, 5, 1000)
        for n in [int(10 ** k) for k in [4, 5]]:
            params = [(n, gr) for gr in GR]
            I = list(tqdm(executor.map(
                        calculate_I,
                        params
                    ),total=len(params),desc=f"n = 10**{int(np.log10(n)+.5)}"))
            plot_I_gr(GR, I, n, 'brown' if n==10**4 else 'magenta', ax)

        GR = np.linspace(0, 5, 100)
        for n in [int(10 ** k) for k in range(6, 8)]:
            params = [(n, gr) for gr in GR]
            I = list(tqdm(executor.map(
                        calculate_I,
                        params
                    ),total=len(params),desc=f"n = 10**{int(np.log10(n)+.5)}"))
            plot_I_gr(GR, I, n, 'purple' if n == 10 ** 6 else 'olive', ax)

    X = np.linspace(0, 1, 1000)
    imax = [0.1957 * (1 - gr) ** 2 *5 for gr in X]
    ax.semilogy(X, imax, ls=':', color='black')
    ax.set_xlabel("$g_r$")
    ax.set_xticks([0, 1, 2, 3, 4, 5])
    ax.set_yticks([1, 10**-6])
    ax.set_ylabel(r"$I_\mathrm{max}/I_\mathrm{SR}$")

    # colorbar
    ax = axes[4]

    sm = mpl.cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = plt.colorbar(sm, cax=ax, orientation="horizontal")
    cbar.set_label("$g_r$")

    output.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output, bbox_inches="tight")
    plt.close(fig)
    print(output)

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Plot numerical panels (a-d) of Fig. 7.")
    parser.add_argument("output", nargs="?", type=Path,
                        default=Path(__file__).resolve().parent / "output/fig7_numeric.pdf")
    parser.add_argument("--raw-data-dir", type=Path,
                        help="Directory with the original full population histories, used only on a summary cache miss.")
    parser.add_argument("--regenerate-trajectories", action="store_true",
                        help="Allow the expensive population simulation if a summary and its source are both missing.")
    args = parser.parse_args()
    fig7(args.output, args.raw_data_dir, args.regenerate_trajectories)
