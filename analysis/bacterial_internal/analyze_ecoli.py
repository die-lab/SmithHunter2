import csv, sys, collections, statistics as st
csv.field_size_limit(10**9)
want={"Spot42","RyhB","GcvB","OmrA","OmrB","FnrS","CyaR","RprA"}
V=[r for r in csv.DictReader(open("validation.tsv"),delimiter="\t") if r["candidate"] in want]
print("pairs for the 8 sRNAs:",len(V))
pct=[float(r["intarna_pct"]) for r in V if r["intarna_pct"] not in ("","NA")]
print("scored:",len(pct),"median pct %.3f"%st.median(pct),"top5%%: %d top10%%: %d top1%%: %d"%(sum(p<=0.05 for p in pct),sum(p<=0.10 for p in pct),sum(p<=0.01 for p in pct)))
# AUC = 1 - mean percentile (approx)
print("AUC approx %.3f"%(1-st.mean(pct)))
truth=collections.defaultdict(str)
for l in open("../bact_ecoli/truth.tsv"):
    c=l.rstrip("\n").split("\t")
    if c[0] in want: truth[(c[0],c[1])]=c[2].split("|")[0].strip()
W={}
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    if r["candidate"] in want: W[(r["candidate"],r.get("gene",""))]=r
print("\nsRNA\tgene\tclass\tpct\trank\tE\tq\thit_region\trel")
by=collections.defaultdict(list)
for r in sorted(V,key=lambda r:float(r["intarna_pct"] or 1)):
    w=W.get((r["candidate"],r["gene"]),{})
    cl=truth.get((r["candidate"],r["gene"]),"?")
    by[cl].append(float(r["intarna_pct"] or 1))
    print("\t".join([r["candidate"],r["gene"],cl,r["intarna_pct"],r["intarna_rank"],w.get("intarna_E",""),w.get("intarna_q",""),w.get("hit_region",""),w.get("rel_start","")+".."+w.get("rel_end","")]))
for k,v in by.items(): print(k,len(v),"median pct %.3f"%st.median(v),"top10%%: %d"%sum(p<=0.1 for p in v))
print("\ncalled windows (q<=0.05):")
for r in csv.DictReader(open("windows_scored.tsv"),delimiter="\t"):
    if r.get("pass_fdr") in ("1","True","true"): print(r["candidate"],r.get("gene"),r["intarna_E"],r["intarna_q"])
