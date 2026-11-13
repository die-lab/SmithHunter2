import csv, math, collections, statistics as st
csv.field_size_limit(10**9)
SR=["Spot42","RyhB","GcvB","OmrA","OmrB","FnrS","CyaR","RprA"]
rows=collections.defaultdict(list)
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    rows[r["candidate"]].append(r)
truth=set()
for l in open("../bact_ecoli/truth.tsv"):
    c=l.split("\t")
    if c[0] in SR: truth.add((c[0],c[1]))
gcvb_t2="argT livJ livK cycA serA asnB cfa gatY gdhA ilvC lrp panD tcyJ thrL yifK ysgA argP asnA brnQ inaA metQ ndk nlpA tppB".split()
tier2={("GcvB",g) for g in gcvb_t2}
# self-locus windows (known coordinates, strand-agnostic overlap)
selfloc={"GcvB":(2942696,2942901),"RyhB":(3580922,3581016)}
def overlaps(r):
    c=r["candidate"]
    if c not in selfloc: return False
    g=r["genomic"].split(":")[1].split("(")[0].split("-")
    a,b=int(g[0]),int(g[1]); lo,hi=selfloc[c]
    # window span unknown; use hit span plus 300 nt window ~ use hit coordinates only
    return not(b<lo or a>hi)
def gene_best(rs):
    best={}
    for r in rs:
        if r["intarna_E"] in ("","NA"): continue
        e=float(r["intarna_E"]); g=r["gene"]
        if g not in best or e<best[g][0]: best[g]=(e,r)
    return best
def gumbel_p(E,x):  # lower-tail of E modelled as min-Gumbel via -E max-Gumbel, moment fit
    v=[-e for e in E]; m=st.mean(v); s=st.pstdev(v)
    beta=s*math.sqrt(6)/math.pi; mu=m-0.5772*beta
    z=(-x-mu)/beta
    return 1-math.exp(-math.exp(-z))
res=[]; allp=[]
for c in SR:
    b=gene_best(rows[c])
    b={g:v for g,v in b.items() if not overlaps(v[1])}
    E=[v[0] for v in b.values()]
    for g,(e,r) in b.items():
        allp.append((gumbel_p(E,e),c,g,e))
# BH
allp.sort(); N=len(allp)
qs=[0]*N; prev=1
for i in range(N-1,-1,-1):
    prev=min(prev,allp[i][0]*N/(i+1)); qs[i]=prev
print("tests",N)
for thr in (0.05,0.1,0.25):
    k=[(c,g) for (p,c,g,e),q in zip(allp,qs) if q<=thr]
    print("q<=%.2f: %d windows; validated tier1 %d; GcvB tier2 %d"%(thr,len(k),sum(x in truth for x in k),sum(x in tier2 for x in k)))
print("top calls:")
for (p,c,g,e),q in list(zip(allp,qs))[:25]:
    print(c,g,e,"%.2g"%p,"%.3f"%q,"T1" if (c,g) in truth else ("T2" if (c,g) in tier2 else ""))
# AUC of tier1 and tier1+tier2 within each sRNA (gene-level, self-locus removed)
def auc(c,pos):
    b=gene_best(rows[c]); b={g:v for g,v in b.items() if not overlaps(v[1])}
    s=sorted(b,key=lambda g:b[g][0]); n=len(s)
    p=[i for i,g in enumerate(s) if (c,g) in pos]
    if not p: return None,0
    neg=n-len(p)
    # mean fraction of negatives ranked below positive
    a=sum((n-1-i-(len(p)-1-j)) for j,i in enumerate(p))/(len(p)*neg)
    return a,len(p)
print("\nper-sRNA AUC (best window per gene, self-locus windows removed)")
for c in SR:
    a1,n1=auc(c,truth); a2,n2=auc(c,truth|tier2)
    print(c,n1,"%.3f"%a1 if a1 else "-",n2,"%.3f"%a2 if a2 else "-")
pool=[]
for c in SR:
    a,n=auc(c,truth); 
    if a: pool.append((a,n))
print("macro AUC tier1 %.3f ; pair-weighted %.3f"%(st.mean(a for a,n in pool),sum(a*n for a,n in pool)/sum(n for a,n in pool)))
