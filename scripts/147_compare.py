#!/usr/bin/env python3
"""147_compare.py: comparison helper for the clean-environment reproduction check (scripts/147_reproduction_check.sh). Reporting only.

  147_compare.py pair  NEW ARCHIVED LABEL GROUP OUT.tsv   compare one regenerated file with the committed (archived) one and append a row to OUT.tsv
  147_compare.py perm  ROOT OUT.tsv                       recompute permutation-null SD, B and p of every rescaled W from the stored null fits
  147_compare.py report WORK                              write WORK/REPRODUCTION_CHECK.md from the rows and the run record

File comparison by type: .tsv numerically (strings exact, numbers to 1e-9; the 'commit' column of the manifests, which records the clone's HEAD,
is ignored); .png pixel by pixel; .pdf by extracted text (pdftotext; the files embed creation times); .xlsx cell by cell; .docx by plain text
(pandoc); anything else byte by byte. A byte-identical file is always reported as identical.
"""
import os, sys, glob, hashlib, json, math, subprocess, datetime

TOL = 1e-9


def md5(f):
    h = hashlib.md5()
    with open(f, "rb") as fh:
        for b in iter(lambda: fh.read(1 << 20), b""): h.update(b)
    return h.hexdigest()


def cmp_tsv(a, b):
    import pandas as pd
    x = pd.read_csv(a, sep="\t", dtype=str, keep_default_na=False); y = pd.read_csv(b, sep="\t", dtype=str, keep_default_na=False)
    ign = [c for c in ("commit",) if c in x.columns and c in y.columns]
    x = x.drop(columns=ign); y = y.drop(columns=ign)
    if list(x.columns) != list(y.columns): return "differs", f"columns differ: {sorted(set(x.columns) ^ set(y.columns))}"
    if x.shape != y.shape: return "differs", f"shape {x.shape} vs {y.shape}"
    worst, nstr = 0.0, 0; where = ""
    for c in x.columns:
        for i, (u, v) in enumerate(zip(x[c], y[c])):
            if u == v: continue
            try:
                fu, fv = float(u), float(v)
                if math.isnan(fu) and math.isnan(fv): continue
                d = abs(fu - fv)
                if d > worst: worst, where = d, f"{c} row {i + 1}: {u} vs {v}"
            except ValueError:
                nstr += 1; where = where or f"{c} row {i + 1}: '{u[:40]}' vs '{v[:40]}'"
    note = " (commit column ignored)" if ign else ""
    if nstr: return "differs", f"{nstr} text cells differ; e.g. {where}"
    if worst == 0: return "identical content", f"all cells equal{note}"
    if worst <= TOL: return "numerically identical", f"max abs difference {worst:.2e} <= {TOL:g}{note}"
    return "differs", f"max abs difference {worst:.3g} at {where}"


def cmp_png(a, b):
    from PIL import Image, ImageChops
    x = Image.open(a).convert("RGBA"); y = Image.open(b).convert("RGBA")
    if x.size != y.size: return "differs", f"size {x.size} vs {y.size}"
    d = ImageChops.difference(x, y); bbox = d.getbbox()
    if bbox is None: return "identical content", f"all {x.size[0]} x {x.size[1]} pixels equal"
    n = sum(1 for p in d.getdata() if any(p)); mx = max(max(p) for p in d.getdata())
    return "differs", f"{n} of {x.size[0] * x.size[1]} pixels differ (max channel difference {mx}), region {bbox}"


def text_of(f, kind):
    if kind == "pdf": return subprocess.run(["pdftotext", "-layout", f, "-"], capture_output=True, text=True, check=True).stdout
    return subprocess.run(["pandoc", f, "-t", "plain", "--wrap=none"], capture_output=True, text=True, check=True).stdout


def cmp_text(a, b, kind):
    ta, tb = text_of(a, kind), text_of(b, kind)
    if ta == tb: return "identical content", f"extracted text equal ({len(ta.split())} words)"
    la, lb = ta.splitlines(), tb.splitlines(); i = next((k for k, (u, v) in enumerate(zip(la, lb)) if u != v), min(len(la), len(lb)))
    return "differs", f"text differs from line {i + 1}: '{(la[i] if i < len(la) else '')[:60]}' vs '{(lb[i] if i < len(lb) else '')[:60]}'"


def cmp_xlsx(a, b):
    from openpyxl import load_workbook
    x = load_workbook(a, read_only=True, data_only=True); y = load_workbook(b, read_only=True, data_only=True)
    if x.sheetnames != y.sheetnames: return "differs", f"sheets differ: {sorted(set(x.sheetnames) ^ set(y.sheetnames))}"
    cells = 0
    for s in x.sheetnames:
        rx = list(x[s].iter_rows(values_only=True)); ry = list(y[s].iter_rows(values_only=True))
        if len(rx) != len(ry): return "differs", f"sheet '{s}': {len(rx)} vs {len(ry)} rows"
        for i, (u, v) in enumerate(zip(rx, ry)):
            if u != v:
                bad = [(j, p, q) for j, (p, q) in enumerate(zip(u, v)) if p != q and not (isinstance(p, float) and isinstance(q, float) and abs(p - q) <= TOL)]
                if bad or len(u) != len(v): return "differs", f"sheet '{s}' row {i + 1}: {bad[:2]}"
            cells += len(u)
    return "identical content", f"{len(x.sheetnames)} sheets, {cells} cells equal"


def pair(new, old, label, group, out):
    if not os.path.exists(new): st, det = "not produced", "the regenerated file is missing"
    elif not os.path.exists(old): st, det = "no archived file", ""
    elif md5(new) == md5(old): st, det = "identical", "byte-identical"
    else:
        ext = os.path.splitext(new)[1].lower()
        try:
            st, det = {".tsv": lambda: cmp_tsv(new, old), ".png": lambda: cmp_png(new, old), ".pdf": lambda: cmp_text(new, old, "pdf"),
                       ".xlsx": lambda: cmp_xlsx(new, old), ".docx": lambda: cmp_text(new, old, "docx")}.get(ext, lambda: ("differs", "bytes differ"))()
        except Exception as e:  # never stop the check on a comparison error; record it
            st, det = "comparison failed", f"{type(e).__name__}: {e}"
    new_hdr = not os.path.exists(out)
    with open(out, "a") as fh:
        if new_hdr: fh.write("group\tfile\tstatus\tdetail\n")
        fh.write(f"{group}\t{label}\t{st}\t{det}\n")
    print(f"  {st:24s} {label}  {det}")


def perm(root, out):
    """Every rescaled W in the two files the manuscript uses: recompute SD, B and b/B from the stored within-design null fits (as 89 and 83 do)."""
    import pandas as pd
    rows = []
    for f, nulldir in (("results/replication_rfi3_calibrated/estimates_calibrated_all.tsv", "results/rfi_null"), ("results/methane_pooled/estimates_calibrated.tsv", "results/methane_dryrun")):
        E = pd.read_csv(os.path.join(root, f), sep="\t", dtype={"id": str})
        for t, g in E.groupby("table"):
            fs = sorted(glob.glob(os.path.join(root, nulldir, f"{t}_seed1[0-9][0-9][0-9]_targets.tsv")))
            N = pd.concat([pd.read_csv(x, sep="\t", dtype={"id": str}) for x in fs]) if fs else pd.DataFrame(columns=["id", "W", "status"])
            N = N[N.status == "ok"]
            for r in g.itertuples():
                w = N.W[N.id == r.id].astype(float)
                k, sd = len(w), (w.std(ddof=1) if len(w) > 1 else float("nan"))
                p = float((w.abs() >= abs(r.W_raw)).mean()) if k else float("nan")
                ok = k == r.k_null and abs(sd - r.sd_null) <= TOL and abs(p - r.p_perm) <= TOL and abs(r.W_raw / max(1.0, sd) - r.W) <= TOL
                rows.append(dict(file=f, table=t, id=r.id, null_files=len(fs), B_file=r.k_null, B_recomputed=k, sd_file=r.sd_null, sd_recomputed=sd,
                                 p_file=r.p_perm, p_recomputed=p, match=ok))
    R = pd.DataFrame(rows); R.to_csv(out, sep="\t", index=False)
    print(f"  permutation nulls: {int(R.match.sum())} of {len(R)} rescaled W reproduced (SD, B, b/B and W) from {R.drop_duplicates('table').null_files.sum()} stored null fits")


def report(work):
    import pandas as pd
    run = json.load(open(os.path.join(work, "out", "run.json")))
    C = pd.read_csv(os.path.join(work, "out", "compare.tsv"), sep="\t", dtype=str, keep_default_na=False)
    Pm = pd.read_csv(os.path.join(work, "out", "perm.tsv"), sep="\t")
    S = pd.read_csv(os.path.join(work, "out", "steps.tsv"), sep="\t", dtype=str, keep_default_na=False)
    good = {"identical", "identical content", "numerically identical"}
    nbad = int((~C.status.isin(good)).sum()); sbad = int((S.rc != "0").sum()); pbad = int((~Pm.match).sum())
    verdict = "PASS: every regenerated output matches the committed file" if not (nbad or sbad or pbad) else f"DIFFERENCES: {nbad} file(s), {sbad} step(s) and {pbad} permutation row(s) listed below"
    L = [f"# Reproduction check (clean environment)", "",
         f"Draft 20, item 3, requested by Dr. Goswami. Run {run['date']} by scripts/147_reproduction_check.sh; reporting only (no new analysis, and nothing in the working repository was changed).", "",
         f"**Result: {verdict}.**", "",
         "## What was done", "",
         f"1. Fresh clone of {run['url']} (branch master) at commit `{run['commit']}` into a new folder ({run['work']}/repo).",
         f"2. The four data archives downloaded from the published Zenodo record {run['record']} (https://zenodo.org/records/{run['record']}), MD5 checked against the record ({run['archives']}), and unpacked at the clone root.",
         "3. Two new conda environments built only from the clone's environment files: `env/rumen-r.yml` (R analysis and figures) and `env/reporting.yml` (Python reporting, pandoc, rsvg-convert).",
         "4. From the committed results, regenerated and compared with the committed files: the number manifests behind Tables 1 to 5 (scripts 126, 143, 145, 146), Figures 1 to 5 (140 and rsvg-convert, 139, 103) and Additional files 1 to 4 (142).",
         "5. End to end from the processed inputs (genus and species count tables, harmonized per-animal phenotypes): the Amendment 46 model-input dumps (21, 65, 81), the per-cohort CLR partial correlations with 9,999 Freedman-Lane permutations (112), and the 48.1 seven-claim table with its Table 4 patterns (128), compared with the archived outputs.",
         "6. The permutation-null SD, number of permutations and b/B behind every rescaled W in the manuscript, recomputed from the stored null fits.", "",
         "## Machine and software", "",
         "| Item | Value |", "|---|---|"] + [f"| {k} | {v} |" for k, v in run["machine"].items()] + [f"| {k} | {v} |" for k, v in run["software"].items()] + ["",
         "## Commands", "", "```", *run["commands"], "```", "",
         "## Steps", "", "| Step | Exit code | Seconds |", "|---|---|---|"] + [f"| {r.step} | {r.rc} | {r.seconds} |" for r in S.itertuples()] + ["",
         "## Comparisons", "",
         "Status: identical = byte-identical; identical content = same pixels (PNG), same extracted text (PDF, DOCX; these files embed creation times), same cell values (XLSX) or same table values (TSV); numerically identical = numbers agree to 1e-9. The `commit` column of the manifests records the clone's HEAD and is not compared.", ""]
    for g, x in C.groupby("group", sort=False):
        L += [f"### {g}", "", "| File | Status | Detail |", "|---|---|---|"] + [f"| `{r.file}` | {r.status} | {r.detail.replace('|', '/')} |" for r in x.itertuples()] + [""]
    L += ["### Permutation nulls behind the rescaled W", "",
          f"{int(Pm.match.sum())} of {len(Pm)} target-cohort estimates reproduced (SD, B, b/B and rescaled W to 1e-9) from {int(Pm.drop_duplicates('table').null_files.sum())} stored null fits in `results/rfi_null/` and `results/methane_dryrun/`."]
    if pbad: L += ["", "| Table | id | B file / recomputed | SD file / recomputed | b/B file / recomputed |", "|---|---|---|---|---|"] + [
        f"| {r.table} | {r.id} | {r.B_file} / {r.B_recomputed} | {r.sd_file:.6g} / {r.sd_recomputed:.6g} | {r.p_file:.4g} / {r.p_recomputed:.4g} |" for r in Pm[~Pm.match].itertuples()]
    L += ["", "## Not covered", "",
          "Read processing and taxonomic profiling (quality control, host removal, Kraken2 and Bracken against the GTDB r207 database, DADA2 against SILVA 138.1) were not rerun: they start from about 540 GB of reads and the Kraken2 database, which are not deposited. Their per-sample reports and logs are in the deposit (`read_qc_and_classification.tar.gz`), and the count tables they produced are the processed inputs used in step 5. ANCOM-BC2 model fits and their 200-permutation nulls were not refitted; step 6 recomputes every rescaled W and permutation p from the stored null fits instead.", ""]
    open(os.path.join(work, "REPRODUCTION_CHECK.md"), "w").write("\n".join(L) + "\n")
    print(verdict)


if __name__ == "__main__":
    m = sys.argv[1]
    if m == "pair": pair(*sys.argv[2:7])
    elif m == "perm": perm(sys.argv[2], sys.argv[3])
    elif m == "report": report(sys.argv[2])
    else: sys.exit(__doc__)
