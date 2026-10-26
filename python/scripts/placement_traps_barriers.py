"""Escape barriers of the spurious peaks and slopes of the log density
(docs/cut_posterior_simulation.md, section 2.3).

  uv run python scripts/placement_traps_barriers.py
"""
import numpy as np, jax.numpy as jnp
from scipy.special import betainc
from circfactor.special import log_link_cdf
def wrap(a): return (np.asarray(a)+np.pi)%(2*np.pi)-np.pi
J=200; rng=np.random.default_rng(2)
yes,no,kappa=rng.uniform(-np.pi,np.pi,J),rng.uniform(-np.pi,np.pi,J),rng.gamma(2.0,20.0,J)
omega=3.0; b0=0.8
for nv in [5,30,200]:
    r=np.random.default_rng(3)
    e=wrap(no-b0)**2-wrap(yes-b0)**2
    p=np.exp(np.asarray(log_link_cdf(jnp.asarray(e),jnp.asarray(kappa))))
    y=(r.random(J)<p).astype(float); y[r.permutation(J)[nv:]]=np.nan
    G=20000; g=-np.pi+2*np.pi*(np.arange(G)+.5)/G; o=~np.isnan(y)
    lp=omega*np.cos(g)
    for j in np.flatnonzero(o):
        u=(2*y[j]-1)*(wrap(no[j]-g)**2-wrap(yes[j]-g)**2); lp=lp+np.log(betainc(kappa[j],kappa[j],(np.pi**2+u)/(2*np.pi**2)))
    left,right=np.roll(lp,1),np.roll(lp,-1)
    peaks=np.flatnonzero((lp>left)&(lp>=right)); mins=np.sort(np.flatnonzero((lp<left)&(lp<=right)))
    top=peaks[np.argmax(lp[peaks])]
    print(f"{nv} votes: {len(peaks)} peaks; spurious peaks: (location, height below top, escape barrier = drop to the lower adjacent minimum, basin width)")
    out=[]
    for pk in peaks:
        if pk==top: continue
        k=np.searchsorted(mins,pk)
        lm, rm = mins[(k-1)%len(mins)], mins[k%len(mins)]
        barrier=lp[pk]-max(lp[lm],lp[rm])
        width=((rm-lm)%G)*2*np.pi/G
        out.append((lp[pk]-lp[top], g[pk], barrier, width))
    for h,loc,b,wd in sorted(out,reverse=True)[:8]:
        print(f"   {loc:+.3f}  {h:8.1f}  barrier {b:7.2f} nats  width {wd:.3f} rad")
    if out: print(f"   median barrier {np.median([o[2] for o in out]):.2f} nats, max {max(o[2] for o in out):.1f}")
    # vote-count makeup: how many votes have their 'cut points' in total

print("\nslope of the log density in the spurious basins (|d log p / d beta|, nats/rad) and implied leapfrog jump eps^2*|grad|/2 at eps=0.05")
for nv in [5,30,200]:
    r=np.random.default_rng(3)
    e=wrap(no-b0)**2-wrap(yes-b0)**2
    p=np.exp(np.asarray(log_link_cdf(jnp.asarray(e),jnp.asarray(kappa))))
    y=(r.random(J)<p).astype(float); y[r.permutation(J)[nv:]]=np.nan
    G=20000; g=-np.pi+2*np.pi*(np.arange(G)+.5)/G; o=~np.isnan(y)
    lp=omega*np.cos(g)
    for j in np.flatnonzero(o):
        u=(2*y[j]-1)*(wrap(no[j]-g)**2-wrap(yes[j]-g)**2); lp=lp+np.log(betainc(kappa[j],kappa[j],(np.pi**2+u)/(2*np.pi**2)))
    grad=np.abs(np.gradient(lp,g))
    away=np.abs(wrap(g-g[np.argmax(lp)]))>0.5
    med=np.median(grad[away])
    print(f"  {nv:>3} votes: median slope away from the main peak {med:7.1f}  -> jump ~{0.05**2*med/2:.3f} rad per half-step kick (x up to 10 leapfrog steps)")
