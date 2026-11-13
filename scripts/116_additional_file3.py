#!/usr/bin/env python3
"""116: Additional file 3 (deviation table and the full analysis plan) as a Word file, from docs/DEVIATIONS.md and plan/ANALYSIS_PLAN.md.

Part 1 is docs/DEVIATIONS.md from its "How to read this file" heading on; Part 2 is plan/ANALYSIS_PLAN.md, unedited, with every heading
moved down one level under "Part 2". Styles come from the previous Additional file 3 (--reference-doc), so the look is unchanged.
Run: python3 scripts/116_additional_file3.py   (needs pandoc >= 2.11; RUMEN_ROOT defaults to /data/rumen_multicohort)
Writes results/additional_files/Additional_file_3_analysis_plan.docx (git keeps earlier versions).
"""
import os, re, shutil, subprocess, tempfile, zipfile

ROOT = os.environ.get("RUMEN_ROOT", "/data/rumen_multicohort")
P = lambda *a: os.path.join(ROOT, *a)
OUT = P("results", "additional_files", "Additional_file_3_analysis_plan.docx")
TITLE = "Additional file 3. Analysis plan, amendments and deviation table"

dev = open(P("docs", "DEVIATIONS.md"), encoding="utf-8").read()
i = dev.index("# How to read this file"); part1 = dev[i:].rstrip() + "\n"
plan = open(P("plan", "ANALYSIS_PLAN.md"), encoding="utf-8").read()
part2 = "\n".join(("#" + l) if re.match(r"^#{1,5} ", l) else l for l in plan.splitlines()) + "\n"
md = part1 + "\n\n# Part 2. Analysis plan with all amendments (unedited)\n\n" + part2

with tempfile.TemporaryDirectory() as td:
    src = os.path.join(td, "af3.md"); open(src, "w", encoding="utf-8").write(md)
    ref = os.path.join(td, "reference.docx"); shutil.copy(OUT, ref)          # previous version supplies the styles
    tmp = os.path.join(td, "af3.docx")
    subprocess.run(["pandoc", src, "-f", "markdown-raw_html", "-t", "docx", "--reference-doc", ref, "-o", tmp], check=True)
    # set the document title in the file properties only (no title paragraph, as before)
    fixed = os.path.join(td, "af3_fixed.docx")
    with zipfile.ZipFile(tmp) as zi, zipfile.ZipFile(fixed, "w", zipfile.ZIP_DEFLATED) as zo:
        for it in zi.infolist():
            b = zi.read(it.filename)
            if it.filename == "docProps/core.xml":
                s = b.decode("utf-8")
                s = re.sub(r"<dc:title>.*?</dc:title>|<dc:title/>", f"<dc:title>{TITLE}</dc:title>", s) if "dc:title" in s else s.replace("</cp:coreProperties>", f"<dc:title>{TITLE}</dc:title></cp:coreProperties>")
                b = s.encode("utf-8")
            zo.writestr(it, b)
    shutil.move(fixed, OUT)
print("wrote", OUT, "| Part 1 from docs/DEVIATIONS.md, Part 2 plan/ANALYSIS_PLAN.md (", len(plan.splitlines()), "lines )")
