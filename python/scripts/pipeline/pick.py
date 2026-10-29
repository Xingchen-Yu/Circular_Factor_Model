"""Best chain among results/full/c* and results/brute/c* for House c, by mean loglik of the saved draws.
Prints '<dir>' of the best chain (e.g. results/brute/c6), then one line per chain."""
import sys, glob
from circfactor.model import Posterior
c = int(sys.argv[1]); ll = {}
for r in sorted(glob.glob(f"results/full/c*/H{c}/posterior.npz") + glob.glob(f"results/brute/c*/H{c}/posterior.npz")):
    p = Posterior.load(r); ll[r.rsplit("/", 2)[0]] = p.loglik_trace[-p.beta.shape[0]:].mean()
best = max(ll, key=ll.get)
print(best)
for k, v in sorted(ll.items(), key=lambda kv: -kv[1]): print(f"  {k:18s} {v:10.0f}{'  <- best' if k == best else ''}")
