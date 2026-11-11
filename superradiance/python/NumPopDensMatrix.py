import scipy as sc
import scipy.integrate
import scipy.stats
from scipy.special import gammaln

import numpy as np
import mpmath as mp
from pathlib import Path as Path

from concurrent.futures import ProcessPoolExecutor as PPE
from tqdm import tqdm

T_global = None
args_global = None
Tmax_global = None

def init_worker(args, Tmax, T):
    global args_global, Tmax_global, T_global
    args_global = args
    Tmax_global = Tmax
    T_global = T

def calc_In(n):
    rho = np.zeros(n + 1)
    rho[n] = 1.0

    h0 = np.array([-args_global["Gs"] * m * (n - m + 1) for m in range(n + 1)])
    h1 = np.array([+args_global["Gs"] * (m + 1) * (n - m) for m in range(n + 1)])

    solver = sc.integrate.RK45(
        lambda t, x: h0 * x + h1 * np.roll(x, -1),
        0,
        rho,
        Tmax_global,
        Tmax_global / args_global["nt"]
    )
    intensity = np.array([m * (n - m + 1) for m in range(n + 1)])
    res = []
    for t in T_global:
        while solver.t < t and solver.status == "running":
            solver.step()
        res.append(intensity @ solver.y)
    return np.array(res)

class NumPopDensMatrix:
    def __init__(self, args):
        file = str(Path(__file__).resolve().parents[1] / "data/NumPopDensMatrix-nmax-")+str(args["nmax"])+"-Gs-"+str(args["Gs"])+"-nt-"+str(args["nt"])+"-n-"+str(args["n"])+"-"
        if "Tmax" in args:
            file = file + "Tmax-" + str(args["Tmax"]) + "-"
        if Path(file+"In.npy").is_file():
            self.In = np.load(file+"In.npy")
            self.T = np.load(file+"T.npy")
        else:
            if "Tmax" in args:
                Tmax = args["Tmax"]
            else:
                Tmax = float(3*mp.log(args["n"])/(args["n"]*args["Gs"]))

            self.In = [] # In[n, t]
            self.T = np.linspace(0, Tmax, args["nt"]+1)

            print("n-it: "+str(args["nmax"]+1))
            with PPE(
                initializer = init_worker,
                initargs = (args, Tmax, self.T)
            ) as executor:
                self.In = np.array(list(tqdm(executor.map(
                    calc_In,
                    range(args["nmax"]+1),
                    chunksize = 10
                ))))
            print("Done!")
            print()
            self.In = np.array(self.In)
            self.T = np.array(self.T)


            np.save(Path(file+"In.npy"), self.In)
            np.save(Path(file+"T.npy"), self.T)

    def intensity(self,args):
        if args["statistics"] == "Fock":
            return self.In[args["n"],:]
        elif args["statistics"] == "Coherent":
            return np.array([sc.stats.poisson.pmf(k, args["n"]) for k in range(len(self.In[:,0]))]) @ self.In
        elif args["statistics"] == "Cat":
            if args.get("normalized", True):
                n = float(mp.findroot(lambda x: x * mp.tanh(x) - args["n"], args["n"]))
            else:
                n = args["n"]
            return np.array([(1-k%2)*(1+np.tanh(n))*sc.stats.poisson.pmf(k, n) for k in range(len(self.In[:,0]))]) @ self.In
        elif args["statistics"] == "Gaussian":
            n = np.array([sc.stats.norm.pdf(k, loc=args["n"], scale=args["sigma"]) for k in range(len(self.In[:,0]))])
            n = n/np.sum(n)
            if args.get("normalized", True):
                n_exp = n @ np.arange(n.size)
                n = n * args["n"]/n_exp
                n[0]= 1 - (sum(n)-n[0])
                return n @ self.In
            else:
                return n @ self.In
        elif args["statistics"] == "Thermal":
            d = np.log((args["n"]+1)/args["n"])
            p = np.array([np.exp(-d*k) for k in range(self.In.shape[0])])
            p = p/np.sum(p)
            return p @ self.In
        elif args["statistics"] == "Squeezed":
            r = np.arcsinh(np.sqrt(args["n"]))
            p = np.zeros(self.In.shape[0])
            n = np.arange(self.In.shape[0])[::2] // 2
            log_p = (
                gammaln(2*n+1) # ln((2n)!)
                -2*gammaln(n+1) # -2ln(n!)
                -2*n*np.log(2) # -2n ln(2)
                +2*n*np.log(np.tanh(r)) # -2n ln(tanh(r))
                -(r+np.log(1+np.exp(-2*r))-np.log(2)) # -ln(cosh(r)) = -(r+ln(1+exp(-2r))-ln(2))
                )
            p[::2] = np.exp(log_p)
            return p @ self.In
        elif args["statistics"] == "Lorentz":
            gamma = args["gamma"]*args["n"]
            p = np.array([gamma/(np.pi*((n-args["n"])**2+gamma**2))for n in range(self.In.shape[0])])
            n_exp = p @ np.arange(p.size)
            p = p * args["n"] / n_exp
            p[0] = 1 - sum(p[1:])
            return p @ self.In
        else:
            return None
