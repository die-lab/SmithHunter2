import csv,re,collections,statistics as st,bisect
csv.field_size_limit(10**9)
SR=["Spot42","RyhB","GcvB","OmrA","OmrB","FnrS","CyaR","RprA"]
nm={"spf":"Spot42","ryhb":"RyhB","gcvb":"GcvB","omra":"OmrA","omrb":"OmrB","fnrs":"FnrS","cyar":"CyaR","rpra":"RprA"}
# window spans
span={}
for r in csv.DictReader(open("regions.tsv"),delimiter="\t"):
    xs=[]
    for b in r["blocks"].split(","):
        a,c=b.split("-"); xs.append((int(a),int(c)))
    span[r["gene"].lower()]=(r["strand"],xs)
def ov(xs,a,b): return any(not(y<a or x>b) for x,y in xs)
def genes(x): return [g.lower() for g in x.split(".") if g and g not in("IGT","EST3UTR","IGR","AS","5UTR","3UTR") and not g.isdigit()]
pairs=collections.defaultdict(lambda:{"libs":0,"frag":0,"inwin":False,"rows":0})
for r in csv.reader(open("../litreview/ril/S2_conditions_rows.tsv"),delimiter="\t"):
    if len(r)<18: continue
    for si,ti in ((0,1),(1,0)):
        s=r[si].lower()
        if s not in nm: continue
        a,b=int(r[12+3*ti]),int(r[13+3*ti]); lo,hi=min(a,b),max(a,b)
        try: nl=int(r[3])
        except: nl=1
        for g in genes(r[ti]):
            k=(nm[s],g); p=pairs[k]; p["rows"]+=1; p["frag"]+=int(r[2]); p["libs"]=max(p["libs"],nl)
            if g in span and ov(span[g][1],lo,hi): p["inwin"]=True
print("RIL pairs",len(pairs)," in window:",sum(v["inwin"] for v in pairs.values()))
best=collections.defaultdict(dict)
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    if r["intarna_E"] in("","NA"): continue
    e=float(r["intarna_E"]);g=r["gene"].lower();c=r["candidate"]
    if g not in best[c] or e<best[c][g]: best[c][g]=e
def auc(c,pos):
    b=best[c]; P=[b[g] for (cc,g) in pos if cc==c and g in b]
    Ns=sorted(v for g,v in b.items() if (c,g) not in pos)
    if len(P)<3: return None,len(P)
    a=sum(bisect.bisect_left(Ns,v) for v in P)
    return 1-a/(len(P)*len(Ns)),len(P)
def rep(name,pos):
    v=[];row=[]
    for c in SR:
        a,n=auc(c,pos); row.append("-" if a is None else "%.2f(%d)"%(a,n))
        if a is not None: v.append(a)
    print("%-42s %s  macro %.3f"%(name," ".join(row),st.mean(v)))
print(" "*42," ".join(SR))
A={k for k,v in pairs.items() if v["libs"]>=2}
rep("RIL>=2 libs, all",A)
rep("RIL>=2 libs, fragment inside window",{k for k in A if pairs[k]["inwin"]})
rep("RIL>=2 libs, fragment outside window",{k for k in A if not pairs[k]["inwin"]})
B={k for k,v in pairs.items() if v["libs"]>=2 and v["frag"]>=30}
rep("RIL>=2 libs & >=30 fragments, inside",{k for k in B if pairs[k]["inwin"]})
print("share of RIL>=2 pairs inside window: %.0f%%"%(100*sum(pairs[k]["inwin"] for k in A)/len(A)))
