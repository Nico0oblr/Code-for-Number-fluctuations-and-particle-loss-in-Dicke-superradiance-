import numpy as np
from pathlib import Path
from scipy.sparse import csr_matrix
from scipy.sparse.linalg import expm_multiply

def i(n, m):
    return n * (n + 1) // 2 + m

def nm(idx):
    n = int((np.sqrt(8 * idx + 1) - 1) / 2)
    m = idx - n * (n + 1) // 2
    return n, m

class NumPopDensLosses:
    """Full Dicke-triangle populations for the fixed initial N in Fig. 7(b)."""
    def __init__(self, args):
        self.n = args["n"] if "n" in args else 100
        Gs = args["Gs"] if "Gs" in args else 1
        Gg = args["Gg"] if "Gg" in args else Gs
        Ge = args["Ge"] if "Ge" in args else Gg
        tnum = args["tnum"] if "tnum" in args else 1000

        Tmax = args["Tmax"] if "Tmax" in args else 3 * np.log(self.n) / (self.n*Gs)
        self.T = np.array(np.linspace(0, Tmax, tnum+1))

        data_dir = Path(args.get("data_dir", Path(__file__).resolve().parents[1] / "data"))
        data_dir.mkdir(parents=True, exist_ok=True)
        file = str(data_dir / "Rho")
        file += "-Gs-" + str(Gs)
        file += "-Gg-" + str(Gg)
        file +=  "-Ge-" + str(Ge)
        file += "-tnum-" + str(tnum)
        file += "-n-" + str(self.n)
        file += "-Tmax-" + str(Tmax)

        if Path(file + "Rho.npy").is_file():
            self.rho = np.load(file + "Rho.npy", mmap_mode="r")
        else:
            N = i(self.n, self.n)+1
            rho0 = np.zeros(N)
            rho0[-1] = 1.0

            rows = []
            cols = []
            values = []
            trace = 0
            for n in range(self.n + 1):
                for m in range(n + 1):
                    # Lindbladian matrix M
                    rows.append(i(n, m))
                    cols.append(i(n, m))
                    values.append(-Gs * m * (n - m + 1) - Gg * (n - m) - Ge * m)
                    trace += values[-1]

                    if n > m >= 0:  # Knock out ground state
                        rows.append(i(n - 1, m))
                        cols.append(i(n, m))
                        values.append(Gg * (n - m))

                    if m > 0:  # Knock out excited & Dicke
                        rows.append(i(n - 1, m - 1))
                        cols.append(i(n, m))
                        values.append(Ge * m)

                        rows.append(i(n, m - 1))
                        cols.append(i(n, m))
                        values.append(Gs * m * (n - m + 1))

            M = csr_matrix((values, (rows, cols)), shape=(N, N))

            self.rho = expm_multiply(M, rho0, start=0, stop=Tmax, num=tnum+1, endpoint=True, traceA=trace)
            np.save(Path(file + "Rho.npy"), self.rho)


    def get_trajectory(self):
        N = np.zeros_like(self.rho[0])
        M = np.zeros_like(self.rho[0])

        for idx in range(len(self.rho[0])):
            n, m = nm(idx)
            N[idx] = n
            M[idx] = m

        return self.rho @ N, self.rho @ M
