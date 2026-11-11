import numpy as np
from concurrent.futures import ProcessPoolExecutor as PPE
from tqdm import tqdm
from pathlib import Path as Path

class MonteCarlo:
    def MC(self, seed):
        N = self.nmax
        m = N

        res = []
        t = 0
        tj = 0

        rng = np.random.default_rng(seed)
        k = 0
        while k <= self.nt:
            r1 = rng.random()
            G = self.Gs * m * (N - m + 1) + self.Ge * m + self.Gg * (N - m)
            tj += -np.log(r1) / G if G > 0 else self.tf

            while t < tj and k <= self.nt:
                res.append(self.Gs * m * (N - m + 1))
                t += self.dt
                k += 1

            r2 = rng.random()
            if N == 0 or G==0:
                tj += self.tf
            elif r2 < self.Gs * m * (N - m + 1) / G:
                m += -1
            elif r2 < (self.Gs * m * (N - m + 1) + self.Ge * m) / G:
                N += -1
                m += -1
            else:
                N += -1
        return res

    def MC_Dicke(self, seed):
        N = self.nmax
        m = N

        res = []
        t = 0
        tj = 0

        rng = np.random.default_rng(seed)
        k = 0
        while k <= self.nt:
            r1 = rng.random()
            G = self.Gs * m * (N - m + 1)
            tj += -np.log(r1) / G if G > 0 else self.tf

            while t < tj and k <= self.nt:
                res.append(self.Gs * m * (N - m + 1))
                t += self.dt
                k += 1

            if m==0 or G==0:
                tj += self.tf
            else:
                m += -1
        return res

    def __init__(self, args):
        self.Gs = args["Gs"] if "Gs" in args else 1
        self.Ge = args["Ge"] if "Ge" in args else self.Gs
        self.Gg = args["Gg"] if "Gg" in args else self.Ge

        self.nt = args["nt"] if "nt" in args else 1000
        self.tf = 3 * np.log(args["N"]) / (args["N"] * self.Gs)
        self.dt = self.tf / self.nt

        self.nmax = args["N"]

        file = str(Path(__file__).resolve().parents[1] / "data/MC")
        file += "-Gs-" + str(self.Gs)
        file += "-Gg-" + str(self.Gg)
        file += "-Ge-" + str(self.Ge)
        file += "-nt-" + str(self.nt)
        file += "-N-"  + str(self.nmax)
        if Path(file + "I.npy").is_file():
            self.I = np.load(file + "I.npy")
            self.T = np.load(file + "T.npy")
        else:
            self.T = np.array(np.linspace(0, self.tf, self.nt + 1))

            seed_sequence = np.random.SeedSequence(42)
            seeds = [
                s.generate_state(1)[0]
                for s in seed_sequence.spawn(args["nMC"] if "nMC" in args else 10000)
            ]

            summe = np.zeros_like(self.T)
            if args.get("sequential", False):
                if args.get("status", False):
                    for seed in tqdm(seeds, total=len(seeds)):
                        summe += self.MC(seed)
                else:
                    for seed in seeds:
                        summe += self.MC(seed)
            else:
                summe = np.zeros_like(self.T)
                if self.Ge==0 and self.Gg==0:
                    with PPE() as executor:
                        for ergebnis in tqdm(executor.map(self.MC_Dicke, seeds, chunksize=1)):
                            summe += ergebnis
                else:
                    with PPE() as executor:
                        for ergebnis in tqdm(executor.map(self.MC, seeds, chunksize=100)):
                            summe += ergebnis

            self.I = np.array(summe/len(seeds))
            np.save(Path(file + "I.npy"), self.I)
            np.save(Path(file + "T.npy"), self.T)

    def getIMax(self, normalized=True):
        if normalized:
            return np.max(self.I)/(self.Gs*self.nmax**2/5)
        else:
            return np.max(self.I)

    def getI(self, normalized=True):
        if normalized:
            return self.I/(self.Gs*self.nmax**2/5)
        else:
            return self.I

    def getT(self, normalized=True):
        if normalized:
            return self.T * self.nmax*self.Gs/np.log(self.nmax)
        else:
            return self.T
