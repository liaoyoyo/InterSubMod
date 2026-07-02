import sys,os,glob,time
sys.path.insert(0,os.path.dirname(__file__))
import decisionflow_wg as D
WT=D.WT
paths=[(p,"chr22") for p in glob.glob(f"{WT}/output/_thresh_cal_2122/**/distance/BERNOULLI/matrix.csv",recursive=True)]
print(f"matrices={len(paths)} NPROC={D.NPROC}")
t0=time.time(); pool=D._mp.Pool(D.NPROC); rt=rm=err=0
for res in pool.imap_unordered(D._work,paths,chunksize=8):
    if res is None: continue
    if isinstance(res,dict) and "_err" in res: err+=1; print("ERR",res["_err"]); continue
    rt+=1; rm+=1
pool.close(); pool.join()
print(f"DONE rt={rt} rm={rm} err={err} elapsed={time.time()-t0:.1f}s ({len(paths)/(time.time()-t0):.0f} loci/s)")
