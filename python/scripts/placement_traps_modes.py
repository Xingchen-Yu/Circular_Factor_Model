"""Local maxima of p(beta | theta, y) and where GHMC reference chains end
(docs/cut_posterior_simulation.md, section 2.3).

  uv run python scripts/placement_traps_modes.py
"""
import numpy as np, jax, jax.numpy as jnp
from scipy.special import betainc
from circfactor.special import log_link_cdf
from circfactor.ghmc import angle_to_circle, circle_to_angle, ghmc_step
from circfactor.model import Data, loglik_matrix
from circfactor.placement import place_legislators
def wrap(a): return (np.asarray(a)+np.pi)%(2*np.pi)-np.pi
J=200; rng=np.random.default_rng(2)
yes,no,kappa=rng.uniform(-np.pi,np.pi,J),rng.uniform(-np.pi,np.pi,J),rng.gamma(2.0,20.0,J)
N=5000; omega=3.0; b0=0.8
def ghmc_ref(y, n_iter=400, seed=5):
    mask=~np.isnan(y); data=Data(jnp.asarray(np.where(mask,2*np.nan_to_num(y)-1,0.))[None], jnp.asarray(mask)[None])
    Y,No,K=map(jnp.asarray,(yes,no,kappa))
    def lp(x):
        th=circle_to_angle(x)
        return jax.vmap(lambda b: loglik_matrix(b[None],Y,No,K,data)[0].sum())(th)+omega*jnp.cos(th)
    @jax.jit
    def body(x,key):
        k1,k2,k3=jax.random.split(key,3)
        x,_=ghmc_step(k3,x,lp,jax.random.uniform(k1,(N,),minval=.005,maxval=.05),jax.random.randint(k2,(N,),1,11),10); return x,None
    key=jax.random.key(seed); x0=angle_to_circle(jax.random.uniform(key,(N,),minval=-np.pi,maxval=np.pi))
    x,_=jax.lax.scan(body,x0,jax.random.split(key,n_iter)); return np.asarray(circle_to_angle(x)), np.asarray(circle_to_angle(x0))
for nv in [5,30,200]:
    r=np.random.default_rng(3)
    e=wrap(no-b0)**2-wrap(yes-b0)**2
    p=np.exp(np.asarray(log_link_cdf(jnp.asarray(e),jnp.asarray(kappa))))
    y=(r.random(J)<p).astype(float); y[r.permutation(J)[nv:]]=np.nan
    G=20000; g=-np.pi+2*np.pi*(np.arange(G)+.5)/G; o=~np.isnan(y)
    lp=omega*np.cos(g)
    for j in np.flatnonzero(o):
        u=(2*y[j]-1)*(wrap(no[j]-g)**2-wrap(yes[j]-g)**2); lp=lp+np.log(betainc(kappa[j],kappa[j],(np.pi**2+u)/(2*np.pi**2)))
    w=np.exp(lp-lp.max()); w/=w.sum()
    # local maxima on the circle and their basins (assign each grid point to the max reached by steepest ascent)
    left,right=np.roll(lp,1),np.roll(lp,-1)
    peaks=np.flatnonzero((lp>=left)&(lp>=right))
    # basin: between neighbouring minima
    mins=np.flatnonzero((lp<=left)&(lp<=right))
    def basin_of(i):
        # walk uphill
        while True:
            l,rr=(i-1)%G,(i+1)%G
            nxt=max([i,l,rr],key=lambda k:lp[k])
            if nxt==i: return i
            i=nxt
    # vectorized-ish: assign by nearest peak with ascent for sampled points
    pk_mass={pk:0.0 for pk in peaks}
    lab=np.empty(G,int)
    # sweep: assign using minima boundaries
    order=np.sort(mins)
    for pk in peaks:
        pass
    for i in range(0,G):
        lab[i]=-1
    # cheaper: label by segments between consecutive minima
    ms=np.sort(mins)
    seg=np.searchsorted(ms,np.arange(G))%len(ms)
    seg_peak={}
    for pk in peaks: seg_peak[seg[pk]]=pk
    ref,start=ghmc_ref(y)
    grid=place_legislators(np.broadcast_to(yes,(N,J)),np.broadcast_to(no,(N,J)),np.broadcast_to(kappa,(N,J)),np.full(N,omega),y[None],seed=4)[:,0]
    def frac_in(samples):
        idx=np.clip(((samples+np.pi)/(2*np.pi)*G).astype(int),0,G-1); s=seg[idx]
        return {sp: np.mean(s==sp) for sp in seg_peak}
    fr,fg,fs=frac_in(ref),frac_in(grid),frac_in(start)
    print(f"\n=== {nv} votes: {len(peaks)} local maxima ===")
    print(" peak(rad)  logdens-max  exact_mass  grid_draws  ghmc_chains  ghmc_starts(uniform)")
    for sp,pk in sorted(seg_peak.items(), key=lambda t:-lp[t[1]]):
        mass=w[seg==sp].sum()
        print(f"  {g[pk]:+.3f}    {lp[pk]-lp.max():9.1f}   {mass:9.2e}   {fg[sp]:8.4f}   {fr[sp]:8.4f}    {fs[sp]:8.4f}")
