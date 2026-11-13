import pickle,glob,csv,collections,statistics as st,math,bisect
csv.field_size_limit(10**9)
SR=["Spot42","RyhB","GcvB","OmrA","OmrB","FnrS","CyaR","RprA"]
RR=list(csv.DictReader(open("regions.tsv"),delimiter="	"))
seen={};reg=[]
for r in RR:
    if r["sequence"] not in seen: seen[r["sequence"]]=len(reg); reg.append(r["region_id"])
uid={r["region_id"]:seen[r["sequence"]] for r in RR}
gene_of={r["region_id"]:r["gene"] for r in csv.DictReader(open("regions.tsv"),delimiter="\t")}
data={}
for f in sorted(glob.glob("work/windows_checkpoint/seg-*.pkl")):
    d=pickle.load(open(f,"rb"))
    for (c,i),v in d.items():
        E=v.get("intarna_E")
        data[(c,i)]=E
print("entries",len(data)); ks=sorted({k[0] for k in data}); print(ks[:5],len(ks))
# check mapping vs windows_scored
ws={}
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    ws[(r["candidate"],r["region_id"])]=r["intarna_E"]
bad=0;n=0
for i in range(0,len(reg),97):
    for c in SR:
        a=data.get((c,i)); b=ws.get((c,reg[i]))
        if a is None or b in(None,""): continue
        n+=1; bad+= abs(float(b)-a)>1e-6
print("mapping check",n,"mismatch",bad)
T1=set()
for l in open("../bact_ecoli/truth.tsv"):
    c=l.split("\t")
    if c[0] in SR: T1.add((c[0],c[1].lower()))
K=set()
for r in csv.DictReader(open("../litreview/ril/tier2_ecoli.tsv"),delimiter="\t"):
    if r["source"]=="known_S3": K.add((r["srna"],r["gene"].lower()))
def stats(c):
    out={}  # gene -> best stats over windows
    dec=[k for k in ks if k.startswith(c+"__")]
    for rid in gene_of:
        i=uid[rid]
        E=data.get((c,i))
        if E is None: continue
        D=[data[(d,i)] for d in dec if data.get((d,i)) is not None]
        if len(D)<10: continue
        mu=st.mean(D); sd=st.pstdev(D) or 1e-9
        z=(E-mu)/sd           # negative = real binds better than its decoys at the same window
        p=(1+sum(x<=E for x in D))/(1+len(D))
        g=gene_of[rid].lower()
        o=out.get(g)
        out[g]={"E":min(E,o["E"]) if o else E,"z":min(z,o["z"]) if o else z,"p":min(p,o["p"]) if o else p, "dmu":min(E-mu,o["dmu"]) if o else E-mu} if True else None
    return out
def auc(vals,pos):
    # lower = better; vals: gene->score
    s=sorted(vals.values()); P=[vals[g] for (c,g) in pos if g in vals]
    N=[v for g,v in vals.items() if (c0,g) not in pos]
    Ns=sorted(N)
    if len(P)<3: return None,len(P)
    a=sum(bisect.bisect_left(Ns,v)+0.5*(bisect.bisect_right(Ns,v)-bisect.bisect_left(Ns,v)) for v in P)
    return 1-a/(len(P)*len(Ns)),len(P)
print("sRNA\tn_pos\tAUC_E\tAUC_z(paired)\tAUC_p\tAUC_E-mu\t| curated n\tAUC_E\tAUC_z")
M=collections.defaultdict(list)
for c0 in SR:
    S=stats(c0)
    row=[c0]
    for name,pos in (("t1",T1),("k",K|T1)):
        pp={(c,g) for (c,g) in pos if c==c0}
        res=[]
        for key in ("E","z","p","dmu"):
            a,n=auc({g:v[key] for g,v in S.items()},pp); res.append(a)
        if name=="t1": row+= [str(n)]+["%.3f"%x for x in res]
        else: row+= ["|",str(n),"%.3f"%res[0],"%.3f"%res[1]]
        for k,x in zip(("E","z","p","dmu"),res): M[(name,k)].append(x)
    print("\t".join(row))
for k,v in sorted(M.items()): print(k,"macro %.3f"%st.mean(v))
# how many validated pairs beat all 20 decoys at their window?
S_all={c:stats(c) for c in SR}
b=sum(1 for (c,g) in T1 if g in S_all[c] and S_all[c][g]["p"]<=1/21+1e-9); tot=sum(1 for (c,g) in T1 if g in S_all[c])
allw=sum(1 for c in SR for g,v in S_all[c].items() if v["p"]<=1/21+1e-9); alln=sum(len(S_all[c]) for c in SR)
print("tier1 pairs whose best window beats all decoys at that window: %d of %d; all genes: %d of %d (%.1f%%)"%(b,tot,allw,alln,100*allw/alln))
# z-based p (normal) + BH
import math
def Phi(z): return 0.5*math.erfc(-z/math.sqrt(2))
L=[(Phi(v["z"]),c,g) for c in SR for g,v in S_all[c].items()]
L.sort(); N=len(L); prev=1; q=[0]*N
for i in range(N-1,-1,-1): prev=min(prev,L[i][0]*N/(i+1)); q[i]=prev
for thr in(0.05,0.1,0.25): 
    k=[(c,g) for (p,c,g),qq in zip(L,q) if qq<=thr]
    print("normal-approx paired z, BH q<=%.2f: %d calls, tier1 %d, curated %d"%(thr,len(k),sum(x in T1 for x in k),sum(x in (K|T1) for x in k)))
print("top 15:",[(c,g,"%.1f"%S_all[c][g]["z"],"T1" if (c,g) in T1 else "K" if (c,g) in K else "") for p,c,g in L[:15]])
