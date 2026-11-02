#!/usr/bin/env python3
# Coverage of the pre-specified targets: genomes in the full library vs the two-pass selection
import csv, os, collections
D=os.environ.get("D","/data/refs/k2_gtdb207/")
SEL=os.environ.get("SEL","/data/refs/k2_gtdb207_full/selected_genomes.tsv")
par,rank,name={},{},{}
for l in open(D+"taxonomy/nodes.dmp"):
    p=l.split("\t|\t"); par[p[0]]=p[1]; rank[p[0]]=p[2]
for l in open(D+"taxonomy/names.dmp"):
    p=l.split("\t|\t")
    if "scientific name" in p[3]: name[p[0]]=p[1]
lib=set(l.rstrip("\n").split("\t")[1] for l in open(D+"seqid2taxid.map"))
sel=set(l.strip() for l in open(SEL))
R={"S":"species","G":"genus","F":"family","P":"phylum"}
cl,cs=collections.Counter(),collections.Counter()
for g in lib:
    t=g
    while par.get(t,t)!=t:
        t=par[t]
        if rank.get(t) in R.values():
            cl[(rank[t],name.get(t))]+=1
            if g in sel: cs[(rank[t],name.get(t))]+=1
rows=[]
for r in csv.DictReader(open("plan/target_keys.tsv"),delimiter="\t"):
    if r["status"]!="mapped": rows.append([r["id"],r["target"],r["status"],"","",""]); continue
    rk=R[r["match_rank"]]; ks=r["gtdb_keys"].split(";")
    rows.append([r["id"],r["target"],rk,len(ks),sum(cl[(rk,k)] for k in ks),sum(cs[(rk,k)] for k in ks)])
with open("plan/target_coverage.tsv","w") as f:
    f.write("id\ttarget\trank\tn_keys\tgenomes_library\tgenomes_selected\n")
    for x in rows: f.write("\t".join(map(str,x))+"\n")
m=[x for x in rows if x[3]!=""]
lost=[x for x in m if x[4]>0 and x[5]==0]; none=[x for x in m if x[4]==0]
thin=[x for x in m if 0<x[5]<=2]
print(f"mapped targets: {len(m)}; absent from GTDB: {len(none)}; in GTDB but lost by selection: {len(lost)}; only 1-2 genomes kept: {len(thin)}")
for lab,v in (("ABSENT",none),("LOST",lost),("FEW",thin)):
    for x in v: print(lab,x[0],x[1],f"library={x[4]} selected={x[5]}")
