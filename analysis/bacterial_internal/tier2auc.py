import csv,collections,statistics as st
csv.field_size_limit(10**9)
SR=["Spot42","RyhB","GcvB","OmrA","OmrB","FnrS","CyaR","RprA"]
best=collections.defaultdict(dict)
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    if r["intarna_E"] in("","NA"): continue
    e=float(r["intarna_E"]); g=r["gene"].lower(); c=r["candidate"]
    if g not in best[c] or e<best[c][g]: best[c][g]=e
T1=set()
for l in open("../bact_ecoli/truth.tsv"):
    c=l.split("\t")
    if c[0] in SR: T1.add((c[0],c[1].lower()))
K=set();R2=set();R1=set()
for r in csv.DictReader(open("../litreview/ril/tier2_ecoli.tsv"),delimiter="\t"):
    k=(r["srna"],r["gene"].lower())
    if r["source"]=="known_S3": K.add(k)
    else:
        if int(r["libs"])>=2: R2.add(k)
        R1.add(k)
def auc(c,pos,neg_excl):
    b=best[c]; s=sorted(b,key=lambda g:b[g]); n=len(s)
    P=[i for i,g in enumerate(s) if (c,g) in pos]
    ex=sum(1 for g in s if (c,g) in neg_excl and (c,g) not in pos)
    if len(P)<3: return None,len(P)
    # negatives = genes neither in pos nor in neg_excl (other known/RIL pairs removed from negatives)
    N=[i for i,g in enumerate(s) if (c,g) not in pos and (c,g) not in neg_excl]
    import bisect
    Ns=sorted(N)
    a=sum(bisect.bisect_left(Ns,i) for i in P)  # negatives ranked better than positive
    return 1-a/(len(P)*len(Ns)),len(P)
sets={"tier1(36 truth)":(T1,set()),"known S3 (curated)":(K|T1,set()),"RIL-seq >=2 libs":(R2,set()),"RIL-seq any":(R1,set()),"known+RIL>=2, RIL-any removed from negatives":(K|T1|R2,R1)}
print("set\t"+"\t".join(SR)+"\tmacro")
for name,(pos,ex) in sets.items():
    v=[];row=[]
    for c in SR:
        a,n=auc(c,pos,ex); row.append("-" if a is None else "%.2f(%d)"%(a,n))
        if a is not None: v.append(a)
    print(name+"\t"+"\t".join(row)+"\t%.3f"%st.mean(v))
