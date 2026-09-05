"""Re-run only the draws admissible at gamma=0, at a given leak rate.
Reuses fused_sweep's cube/to_param/run_model/admissible, so the parameter
vectors are identical to the main evaluation. Resumable by line count of OUT."""
import os, sys, json, time
import numpy as np
sys.path.insert(0, os.path.join(os.environ['HOME'], 'mnt/delay/scripts'))
sys.argv = ['fused_sweep.py', '1200', '4000', '/dev/null', '0.0']
import fused_sweep as FS
GAMMA=float(os.environ['GAMMA']); OUT=os.environ['OUT']
IDX=json.load(open(os.path.join(os.environ['HOME'],'mnt/delay/scripts/results/leak/adm_idx.json')))
sh_n=int(os.environ.get('SHARD_N','1')); sh_i=int(os.environ.get('SHARD_I','0'))
IDX=[d for j,d in enumerate(IDX) if j % sh_n == sh_i]
rng=np.random.default_rng(FS.SEED0); cube=FS.lhs(FS.N_DRAWS, FS.DIMS, rng)
done=sum(1 for _ in open(OUT)) if os.path.exists(OUT) else 0
print(f"gamma={GAMMA} shard {sh_i}/{sh_n}: {len(IDX)} draws, resuming at {done}", flush=True)
with open(OUT,'a') as f:
    for i in IDX[done:]:
        t0=time.time()
        p,bf=FS.to_param(cube[i])
        beta2=np.exp(FS.LO['log_beta']+bf*(FS.HI['log_beta']-FS.LO['log_beta']))
        Delta=(FS.FIXED['delay']+FS.FIXED['LRS'])-FS.FIXED['Lmotion']
        r1=FS.run_model(p,0.0,1,FS.SEED0+1000*i,gamma=GAMMA)
        r2=FS.run_model(p,beta2/Delta,2,FS.SEED0+1000*i+500,gamma=GAMMA)
        a1,a2=FS.admissible(r1),FS.admissible(r2)
        f.write(json.dumps(dict(draw=i,params=p,beta=beta2,alpha=beta2/Delta,gamma=GAMMA,
            model1=dict(admissible=a1,rt=FS.rt_effects(r1),slope=FS.slope_criterion(r1),met=FS.met_invariance(r1)),
            model2=dict(admissible=a2,rt=FS.rt_effects(r2),slope=FS.slope_criterion(r2),met=FS.met_invariance(r2))))+"\n")
        f.flush()
        print(f"draw {i}: adm1={a1} adm2={a2} ({time.time()-t0:.1f}s)", flush=True)
print("SHARD DONE", flush=True)
