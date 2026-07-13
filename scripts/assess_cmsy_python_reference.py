#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Data-limited CMSY/BSM (Schaefer) type assessment + CPUE.
Estimates r, k, MSY, Bmsy, Fmsy, B/Bmsy, F/Fmsy, depletion, Kobe status.
Monte-Carlo (SIR): samples (r,k), simulates biomass from catches,
filters viable trajectories, weights by the lognormal likelihood of CPUE.
"""
import pandas as pd, numpy as np, json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
rng=np.random.default_rng(42)
PROC="data/processed"; TAB="output/tables"; FIG="output/figures"
os.makedirs(TAB,exist_ok=True)

# resilience -> prior on r (Froese et al. 2017)
R_PRIOR={"Octopus vulgaris":(0.6,1.5),   # high (annual, r-strategist)
         "Penaeus notialis":(0.6,1.5),   # high (annual shrimp)
         "Cymbium spp":(0.2,0.8)}        # medium (gastropod)

def cmsy(catch, years, cpue=None, cpue_years=None, r_range=(0.2,1.0),
         start_dep=(0.4,0.9), end_dep=(0.05,0.6), n=300000):
    C=np.asarray(catch,float); T=len(C)
    kmax=200*C.max(); kmin=C.max()          # wide bounds on k
    r=np.exp(rng.uniform(np.log(r_range[0]),np.log(r_range[1]),n))
    k=np.exp(rng.uniform(np.log(kmin),np.log(kmax),n))
    b0=rng.uniform(*start_dep,n)             # B1/k
    B=np.empty((n,T)); B[:,0]=b0*k
    viable=np.ones(n,bool)
    for t in range(T-1):
        Bt=B[:,t]
        Bn=Bt+r*Bt*(1-Bt/k)-C[t]
        viable&=(Bn>0)&(Bn<k*1.05)
        Bn=np.clip(Bn,1e-6,k*1.05)
        B[:,t+1]=Bn
    dep=B[:,-1]/k
    viable&=(dep>=end_dep[0])&(dep<=end_dep[1])
    w=np.where(viable,1.0,0.0)
    # CPUE likelihood (biomass ~ index/q) on viable trajectories
    if cpue is not None:
        cy=np.asarray(cpue,float); idx=[years.index(y) for y in cpue_years if y in years]
        obs=np.log(cy[[cpue_years.index(years[i]) for i in idx]])
        Bsub=B[:,idx]
        vv=viable & (Bsub>0).all(1)
        ll=np.full(n,-1e9)
        Bl=np.log(np.clip(Bsub,1e-6,None))
        q=(obs[None,:]-Bl).mean(1)                       # log q per draw
        resid=obs[None,:]-(Bl+q[:,None])
        sig=resid[vv].std()+1e-3
        ll_v=-0.5*np.sum((resid/sig)**2,1)-len(idx)*np.log(sig)
        ll[vv]=ll_v[vv]
        w=np.where(vv,np.exp(ll-np.nanmax(ll[vv])),0.0)
    if w.sum()==0: raise RuntimeError("no viable trajectory — widen priors")
    w/=w.sum()
    # SIR resampling
    keep=rng.choice(n,size=min(20000,(w>0).sum()*5 or 20000),p=w)
    return dict(r=r[keep],k=k[keep],B=B[keep],dep=(B[keep,-1]/k[keep]))

def summ(x):
    return dict(median=float(np.median(x)),lo=float(np.percentile(x,2.5)),hi=float(np.percentile(x,97.5)))

catch=pd.read_csv(f"{PROC}/captures_totales_annuelles.csv")
rows=[]; traj={}
for sp in ["Octopus vulgaris","Penaeus notialis","Cymbium spp"]:
    g=catch[catch.espece==sp].sort_values("annee")
    years=g.annee.astype(int).tolist(); C=g.capture_t.values
    # standardized CPUE
    fn=f"{TAB}/cpue_std_{sp.split()[0].lower()}.csv"
    cp=pd.read_csv(fn) if os.path.exists(fn) else None
    cpv=cpy=None
    if cp is not None:
        cp=cp[cp.annee.isin(years)]; cpv=cp.index_scaled.values; cpy=cp.annee.astype(int).tolist()
    # depletion priors: start lightly exploited, end based on recent catch trend
    recent=C[-5:].mean()/C.max()
    end=(0.05,0.5) if recent<0.5 else (0.2,0.7)
    out=cmsy(C,years,cpv,cpy,r_range=R_PRIOR[sp],start_dep=(0.5,0.95),end_dep=end)
    r,k=out["r"],out["k"]; MSY=r*k/4; Bmsy=k/2; Fmsy=r/2
    Bcur=out["B"][:,-1]; BBmsy=Bcur/Bmsy
    Fcur=C[-1]/Bcur; FFmsy=Fcur/Fmsy
    traj[sp]=out
    rows.append(dict(stock=sp,
        r=summ(r)["median"], k=summ(k)["median"], MSY=summ(MSY)["median"],
        MSY_lo=summ(MSY)["lo"],MSY_hi=summ(MSY)["hi"],
        Bmsy=summ(Bmsy)["median"], Fmsy=summ(Fmsy)["median"],
        depletion=summ(out["dep"])["median"],
        B_Bmsy=summ(BBmsy)["median"], F_Fmsy=summ(FFmsy)["median"],
        C_last=float(C[-1]), year_last=years[-1]))
    print(f"{sp}: r={np.median(r):.2f} k={np.median(k):.0f} MSY={np.median(MSY):.0f}t "
          f"B/Bmsy={np.median(BBmsy):.2f} F/Fmsy={np.median(FFmsy):.2f} dep={np.median(out['dep']):.2f}")

res=pd.DataFrame(rows)
res.to_csv(f"{TAB}/stock_assessment_BRP.csv",index=False)

# ---- Kobe plot ----
fig,ax=plt.subplots(figsize=(7,6))
ax.axvspan(0,1,ymin=0.5,color="#FF7F7F",alpha=.35); ax.axvspan(0,1,ymax=0.5,color="#FFF59D",alpha=.5)
ax.axvspan(1,3,ymin=0.5,color="#FFF59D",alpha=.5); ax.axvspan(1,3,ymax=0.5,color="#90EE90",alpha=.4)
for _,r in res.iterrows():
    ax.scatter(r.B_Bmsy,r.F_Fmsy,s=90,zorder=5)
    ax.annotate(r.stock.split()[0],(r.B_Bmsy,r.F_Fmsy),xytext=(5,5),textcoords="offset points",fontsize=10)
ax.axhline(1,ls="--",c="grey"); ax.axvline(1,ls="--",c="grey")
ax.set_xlim(0,2.5); ax.set_ylim(0,max(3,res.F_Fmsy.max()*1.2))
ax.set_xlabel("B/Bmsy"); ax.set_ylabel("F/Fmsy"); ax.set_title("Kobe diagram — stock status (2024)")
plt.tight_layout(); plt.savefig(f"{FIG}/fig8_kobe.png",dpi=130); plt.close()

# ---- biomass trajectories ----
fig,axs=plt.subplots(1,3,figsize=(16,4.5))
for ax,sp in zip(axs,traj):
    g=catch[catch.espece==sp].sort_values("annee"); years=g.annee.astype(int).values
    B=traj[sp]["B"]; k=traj[sp]["k"]
    med=np.median(B,0); lo=np.percentile(B,2.5,0); hi=np.percentile(B,97.5,0)
    Bmsy=np.median(k)/2
    ax.fill_between(years,lo,hi,alpha=.25,color="#2c7fb8"); ax.plot(years,med,color="#12507b")
    ax.axhline(Bmsy,ls="--",color="#c0392b"); ax.text(years[0],Bmsy,"Bmsy",color="#c0392b",fontsize=8,va="bottom")
    ax.set_title(sp,style="italic",fontsize=10); ax.set_xlabel("Year"); ax.set_ylabel("Biomass (t)")
plt.suptitle("Biomass trajectories (CMSY/BSM Schaefer)",fontweight="bold"); plt.tight_layout()
plt.savefig(f"{FIG}/fig9_biomasse.png",dpi=130); plt.close()

print("\n=== BRP TABLE ===")
print(res.round(2).to_string(index=False))
print("\nfigures 8 (Kobe) & 9 (biomass) OK")
