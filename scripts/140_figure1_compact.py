#!/usr/bin/env python3
"""140: Figure 1 for draft 18, the 111 schematic with the same content and a compact layout (Amendment 48.11, reconciliation item 13: figure and legend
to fit 170 x 225 mm). Only spacing changes: line height 1.10 x size (111: 1.22), panel gaps, table row height, box padding. Text sizes and the overflow check
are as in 111; one wording change for draft 18: panel d reads 'coverage rule met in simulation' (111: 'validated by simulation'). Writes results/figures_draft18/fig1_study_design.svg (refuses to overwrite). Usage: 140_figure1_compact.py [root]
111 docstring follows.
Figure 1 (manuscript): study design schematic. Writes fig1_study_design.svg in the current directory.
Render to PDF/PNG with: node scripts/111_render_svg.js fig1_study_design.svg out.png out.pdf <height>mm   (height printed below)
Amendment 46.5.7: redrawn so that every text is at least 8 pt at 170 mm. The canvas is 482 units wide = 170 mm, so 1 unit = 1 pt.
Text widths are measured with Arial/Liberation Sans advance widths (table below) and every text is checked against the space it is given;
the script stops if any text would overflow. Panel d describes the effect-size estimand behind the exclusion statements: the Amendment 46 CLR partial r
with the modified Knapp-Hartung interval validated in Amendment 47; panel e lists the Amendment 47 interval validation and null-scan calibration."""
import re, math, json, os, sys

W = 482
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#d6d4cc'
BLUE_D, BLUE, BLUE_P = '#104281', '#256abf', '#eaf1fb'
ORNG, ORNG_P = '#b8601d', '#fbf0e6'
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"
S_SM, S_BODY, S_HEAD, S_TITLE, S_LET = 8.0, 8.5, 9.0, 10.0, 12.0   # pt; smallest text 8 pt
assert min(S_SM, S_BODY, S_HEAD, S_TITLE, S_LET) >= 8.0
ADV = json.loads('{"r": {" ": 278, "!": 278, "\\"": 355, "#": 556, "$": 556, "%": 889, "&": 667, "\'": 191, "(": 333, ")": 333, "*": 389, "+": 584, ",": 278, "-": 333, ".": 278, "/": 278, "0": 556, "1": 556, "2": 556, "3": 556, "4": 556, "5": 556, "6": 556, "7": 556, "8": 556, "9": 556, ":": 278, ";": 278, "<": 584, "=": 584, ">": 584, "?": 556, "@": 1015, "A": 667, "B": 667, "C": 722, "D": 722, "E": 667, "F": 611, "G": 778, "H": 722, "I": 278, "J": 500, "K": 667, "L": 556, "M": 833, "N": 722, "O": 778, "P": 667, "Q": 778, "R": 722, "S": 667, "T": 611, "U": 722, "V": 667, "W": 944, "X": 667, "Y": 667, "Z": 611, "[": 278, "\\\\": 278, "]": 278, "^": 469, "_": 556, "`": 333, "a": 556, "b": 556, "c": 500, "d": 556, "e": 556, "f": 278, "g": 556, "h": 556, "i": 222, "j": 222, "k": 500, "l": 222, "m": 833, "n": 556, "o": 556, "p": 556, "q": 556, "r": 333, "s": 500, "t": 278, "u": 556, "v": 500, "w": 722, "x": 500, "y": 500, "z": 500, "{": 334, "|": 260, "}": 334, "~": 584, "→": 1000, "≥": 549, "√": 549, "–": 556, "—": 1000, "₄": 408, "†": 556, "‡": 556, "×": 584, "·": 333, "²": 333}, "b": {" ": 278, "!": 333, "\\"": 474, "#": 556, "$": 556, "%": 889, "&": 722, "\'": 238, "(": 333, ")": 333, "*": 389, "+": 584, ",": 278, "-": 333, ".": 278, "/": 278, "0": 556, "1": 556, "2": 556, "3": 556, "4": 556, "5": 556, "6": 556, "7": 556, "8": 556, "9": 556, ":": 333, ";": 333, "<": 584, "=": 584, ">": 584, "?": 611, "@": 975, "A": 722, "B": 722, "C": 722, "D": 722, "E": 667, "F": 611, "G": 778, "H": 722, "I": 278, "J": 556, "K": 722, "L": 611, "M": 833, "N": 722, "O": 778, "P": 667, "Q": 778, "R": 722, "S": 667, "T": 611, "U": 722, "V": 667, "W": 944, "X": 667, "Y": 667, "Z": 611, "[": 333, "\\\\": 278, "]": 333, "^": 584, "_": 556, "`": 333, "a": 556, "b": 611, "c": 556, "d": 611, "e": 556, "f": 333, "g": 611, "h": 611, "i": 278, "j": 278, "k": 556, "l": 278, "m": 889, "n": 611, "o": 611, "p": 611, "q": 611, "r": 389, "s": 556, "t": 333, "u": 611, "v": 556, "w": 778, "x": 556, "y": 556, "z": 500, "{": 389, "|": 280, "}": 389, "~": 584, "→": 1000, "≥": 549, "√": 549, "–": 556, "—": 1000, "₄": 477, "†": 556, "‡": 556, "×": 584, "·": 333, "²": 333}}')   # advance widths per 1000 em, Liberation Sans regular (r) and bold (b)

out = []; checks = []
def add(s): out.append(s)
def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
PARTS = r'(\*[^*]+\*|_\{[^}]+\}|\^\{[^}]+\})'
def width(s, size, bold=False):
    t = ADV['b' if bold else 'r']; w = 0.0
    for p in re.split(PARTS, s):
        if not p: continue
        if p.startswith('*'): q, f = p[1:-1], 1.0
        elif p[:2] in ('_{', '^{'): q, f = p[2:-1], 0.72
        else: q, f = p, 1.0
        w += f * sum(t.get(c, 600) for c in q)
    return w * size / 1000
def rich(s):
    o = ''
    for p in re.split(PARTS, s):
        if not p: continue
        if p.startswith('*'): o += f'<tspan font-style="italic">{esc(p[1:-1])}</tspan>'
        elif p.startswith('_{'): o += f'<tspan baseline-shift="sub" font-size="72%">{esc(p[2:-1])}</tspan>'
        elif p.startswith('^{'): o += f'<tspan baseline-shift="super" font-size="72%">{esc(p[2:-1])}</tspan>'
        else: o += esc(p)
    return o
def text(x, y, s, size=S_BODY, weight='normal', color=INK, anchor='start', lh=None, maxw=None, tag=''):
    lines = s if isinstance(s, list) else [s]
    lh = lh or size * 1.10
    for i, ln in enumerate(lines):
        if not ln: continue
        if maxw is not None:
            w = width(ln, size, weight == 'bold'); checks.append((w, maxw, ln))
        add(f'<text x="{x:.1f}" y="{y + i*lh:.1f}" font-size="{size}" font-weight="{weight}" fill="{color}" text-anchor="{anchor}">{rich(ln)}</text>')
    return y + len(lines) * lh
def rect(x, y, w, h, fill='white', stroke=GRID, sw=0.7, r=2.5):
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"/>')
def line(x1, y1, x2, y2, color=MUTED, sw=0.8, arrow=False):
    a = ' marker-end="url(#arr)"' if arrow else ''
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{sw}"{a}/>')
def panel(x, y, letter, title, maxw):
    text(x, y, letter, size=S_LET, weight='bold')
    text(x + 11, y, title, size=S_TITLE, weight='bold', maxw=maxw - 11)

M = 4; IW = W - 2 * M          # margins, inner width
LH = S_SM * 1.10               # small-text line height (111: 1.22)

# ================================================================ a. cohorts (table)
y = 13
panel(M, y, 'a', 'Five public cohorts, about 2,400 animals', IW)
cols = [('Cohort', 42), ('n', 27), ('Country', 63), ('Accession', 61), ('Data', 84), ('Trait', 87), ('Design terms', 0)]
cols[-1] = (cols[-1][0], IW - 6 - sum(w for _, w in cols[:-1]))
xs = [M + 6]
for _, w in cols[:-1]: xs.append(xs[-1] + w)
rows = [('Feed efficiency, beef cattle', BLUE_D, BLUE, BLUE_P, [
            ('Nelore', '52', 'Brazil', 'PRJNA987743', 'Shotgun metagenome', 'RFI, continuous', 'diet, slaughter batch'),
            ('Charolais', '30', 'France', 'PRJEB70895', 'Metatranscriptome', 'RFI, high vs low class', 'diet, year'),
            ('Alberta', '567', 'Canada', 'PRJNA393057', '16S rRNA gene', 'RFI, continuous', 'management group, breed')]),
        ('Methane, dairy cows', ORNG, ORNG, ORNG_P, [
            ('Difford', '749', 'Denmark', 'PRJEB28065', '16S rRNA gene', 'CH₄, g per day', 'herd, batch, parity, DIM'),
            ('Wallace', '1,009', 'Europe, 7 farms', 'PRJNA517480', '16S rRNA gene', 'CH₄, g per day', 'farm, parity, DIM')])]
y += 15
for i, (h, _) in enumerate(cols):
    text(xs[i], y, h, size=S_SM, weight='bold', color=INK2, maxw=cols[i][1] - 3)
line(M, y + 3.5, M + IW, y + 3.5, color=GRID, sw=0.7)
y += 5
RH = 10.0
for head, hc, col, pale, rr in rows:
    y += RH
    text(M + 6, y - 1, head, size=S_HEAD, weight='bold', color=hc, maxw=IW)
    top = y + 2.5
    rect(M, top, IW, RH * len(rr) + 3, fill=pale, stroke=col, sw=0.7)
    add(f'<rect x="{M:.1f}" y="{top:.1f}" width="2.8" height="{RH*len(rr)+3:.1f}" fill="{col}"/>')
    for r in rr:
        y += RH
        for i, v in enumerate(r):
            text(xs[i], y, v, size=S_SM, weight='bold' if i == 0 else 'normal', color=INK if i in (0, 4, 5) else INK2, maxw=cols[i][1] - 5)
    y += 6
y += 9
text(M, y, 'All sequence data are public; phenotypes come from the original publications. DIM, days in milk.', size=S_SM, color=INK2, maxw=IW)
a_bottom = y

# ================================================================ b. targets
y = a_bottom + 10
panel(M, y, 'b', 'Published targets', IW)
bw = (IW - 8) / 2
def tset(x, y0, col, pale, head, rows, result):
    h = 14 + LH * 1.08 * len(rows) + 15
    rect(x, y0, bw, h, fill='white', stroke=col, sw=0.8)
    text(x + 6, y0 + 12, head, size=S_HEAD, weight='bold', color=col, maxw=bw - 12)
    yy = y0 + 12 + LH * 1.25
    for code, desc in rows:
        if code:
            add(f'<rect x="{x+6:.1f}" y="{yy-7.4:.1f}" width="10.5" height="9.4" rx="1.5" fill="{pale}" stroke="{col}" stroke-width="0.6"/>')
            text(x + 11.25, yy, code, size=S_SM, weight='bold', color=col, anchor='middle')
        text(x + 21, yy, desc, size=S_SM, color=INK, maxw=bw - 27)
        yy += LH * 1.08
    text(x + 6, yy + 3, result, size=S_BODY, weight='bold', color=INK, maxw=bw - 12)
    return y0 + h
y0 = y + 4
h1 = tset(M, y0, BLUE, BLUE_P, 'Feed efficiency: 76 targets',
          [('T', 'Li et al. 2019: 6 taxa'), ('N', 'Nelore discovery study: 64 genera'), ('C', 'Canonical rumen taxa: 6'), ('', ''), ('', '')],
          '27 testable')
h2 = tset(M + bw + 8, y0, ORNG, ORNG_P, 'Methane: 6 targets + Nelore genera',
          [('M', 'Reported in more than one country: 3'), ('', '*M. gottschalkii* : *M. ruminantium* clade ratio;'), ('', 'Succinivibrionaceae; total archaea'),
           ('S', 'Secondary: 3 genera'), ('N', 'Nelore genera reported for residual methane')],
          '11 testable')
y = max(h1, h2) + 9
text(M, y, ["Testable: mapped feature in at least 10% of a cohort's animals. A cohort that reported a target, or helped choose it,",
            'counts as re-analysis, not replication.'], size=S_SM, color=INK2, lh=LH, maxw=IW)
b_bottom = y + LH

# ================================================================ c. estimation
y = b_bottom + 10
panel(M, y, 'c', 'Estimation in each cohort, for each target', IW)
sy = y + 4; ns = 4; sg = 9
sw_ = (IW - (ns - 1) * sg) / ns
steps = [
    ('1  Profiling', [('b', 'Shotgun and'), ('b', 'metatranscriptome'), ('', 'Host reads removed'), ('', '(minimap2, ARS-UCD2.0)'), ('', 'Kraken 2 + Bracken,'), ('', 'GTDB r207'),
                      ('b', '16S rRNA gene'), ('', 'cutadapt, DADA2,'), ('', 'SILVA 138.1'), ('', 'Genera in ≥10% of animals')]),
    ('2  Detection check†', [('', 'Same animals, table and'), ('', 'model, with the design'), ('', 'factor as the tested term'), ('', ''),
                              ('', 'Null results count only if'), ('', '≥1 genus reaches q < 0.05'), ('', ''), ('', 'Targets mapped by name;'), ('', 'family targets summed')]),
    ('3  ANCOM-BC2', [('', 'genus ~ design terms + trait'), ('', ''), ('', '*W* = log fold change / SE'), ('', ''),
                      ('', 'Continuous traits z-scored'), ('', 'within cohort, and within'), ('', 'farm for Wallace')]),
    ('4  Calibration', [('b', 'Permutation null'), ('', 'Trait shuffled 200 times'), ('', 'within design groups'), ('', ''), ('', ''), ('', ''), ('', ''),
                                    ('', '*W*_{cal} = *W* / max(1, SD_{null})'), ('', '*r* = *W*_{cal} / √(*W*_{cal}^{2} + df)'), ('', 'plus permutation *p*')]),
]
SL = LH * 1.0
sh = 15 + 5 + SL * max(len(r) for _, r in steps) + 3
for i, (head, rows_) in enumerate(steps):
    x = M + i * (sw_ + sg)
    rect(x, sy, sw_, sh, fill='white', stroke=INK2, sw=0.7)
    add(f'<path d="M{x:.1f},{sy+15:.1f} v-12.5 a2.5,2.5 0 0 1 2.5,-2.5 h{sw_-5:.1f} a2.5,2.5 0 0 1 2.5,2.5 v12.5 z" fill="{INK2}"/>')
    text(x + 5, sy + 11, head, size=S_HEAD, weight='bold', color='white', maxw=sw_ - 8)
    yy = sy + 15 + 5 + SL * 0.75
    for kind, s in rows_:
        text(x + 5, yy, s, size=S_SM, weight='bold' if kind == 'b' else 'normal', color=INK if kind == 'b' else INK2, maxw=sw_ - 8)
        yy += SL
    if i < ns - 1:
        line(x + sw_ + 0.8, sy + sh / 2, x + sw_ + sg - 0.8, sy + sh / 2, color=INK2, sw=1.0, arrow=True)
# null glyph in step 4
gx = M + 3 * (sw_ + sg) + 7; gw = sw_ - 14; gy = sy + 15 + 5 + SL * 5.6
def curve(sd, col, dash=None):
    pts = [(gx + gw * k / 60, gy - 17 * math.exp(-0.5 * ((-4 + 8 * k / 60) / sd) ** 2) / sd) for k in range(61)]
    d = 'M' + ' L'.join(f'{a:.1f},{b:.1f}' for a, b in pts)
    da = f' stroke-dasharray="{dash}"' if dash else ''
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="1.0"{da}/>')
line(gx, gy, gx + gw, gy, color=GRID, sw=0.7)
curve(1.0, MUTED, '2,1.5'); curve(2.0, BLUE)
text(gx, gy - 17.5, 'nominal', size=S_SM, color=MUTED)
text(gx, gy - 8.5, 'observed', size=S_SM, color=BLUE)
c_bottom = sy + sh

# ================================================================ d. pooling, e. robustness
y = c_bottom + 10
dw = (IW - 8) / 2; ex = M + dw + 8; ew = IW - dw - 8
panel(M, y, 'd', 'Pooling, replication label and effect size', dw)
panel(ex, y, 'e', 'Robustness checks', ew)
py = y + 4
yy = py + 13
text(M + 6, yy, 'Fisher *z*  →  REML random effects  →  BH within set', size=S_SM, color=INK, maxw=dw - 12)
yy += LH + 4
crit = ['Pooled *q* < 0.05', 'Same sign as published', 'Holds with each cohort left out']
cbw = dw - 12 - 70
for k, c in enumerate(crit):
    rect(M + 6, yy - 7.6 + k * (LH + 3), cbw, LH + 1.5, fill=BLUE_P, stroke=BLUE, sw=0.6, r=1.5)
    text(M + 10, yy + k * (LH + 3), c, size=S_SM, color=INK, maxw=cbw - 8)
cmid = yy - 7.6 + (3 * (LH + 3) - 3) / 2
line(M + 6 + cbw + 2, cmid, M + 6 + cbw + 10, cmid, color=INK2, sw=1.0, arrow=True)
rect(M + 6 + cbw + 12, cmid - 8, dw - 12 - cbw - 12, 16, fill=BLUE_D, stroke=BLUE_D)
text(M + 6 + cbw + 12 + (dw - 12 - cbw - 12) / 2, cmid + 3, 'Replicated', size=S_BODY, weight='bold', color='white', anchor='middle', maxw=dw - 12 - cbw - 14)
yy += 3 * (LH + 3) + 5
yy = text(M + 6, yy, ['Otherwise: single-cohort support, opposite sign,', 'not replicated, or untestable'], size=S_SM, color=INK2, lh=LH, maxw=dw - 12)
yy += 6
text(M + 6, yy, 'Effect size and benchmarks', size=S_SM, weight='bold', maxw=dw - 12)
yy = text(M + 6, yy + LH, ['CLR partial *r* per cohort; modified Knapp-Hartung', 'pooled intervals, coverage rule met in simulation;',
                           'does the 95% interval exclude *r* = 0.15 or 0.30', 'in the published direction?'], size=S_SM, color=INK2, lh=LH, maxw=dw - 12)
d_end = yy - LH + 6
def bullets(y0, head, items, filled):
    text(ex + 7, y0, head, size=S_SM, weight='bold', maxw=ew - 12)
    y0 += LH + 2
    for s in items:
        if filled: add(f'<circle cx="{ex+10:.1f}" cy="{y0-2.6:.1f}" r="1.9" fill="{INK2}"/>')
        else: add(f'<circle cx="{ex+10:.1f}" cy="{y0-2.6:.1f}" r="1.9" fill="white" stroke="{INK2}" stroke-width="0.8"/>')
        text(ex + 16, y0, s, size=S_SM, color=INK2, maxw=ew - 20); y0 += LH
    return y0
yy = bullets(py + 13, 'Specified before the affected results',
             ['Wallace without its Swedish farm', 'Methane per kg dry matter intake', 'Wallace as a fourth RFI cohort', 'Alberta without possible discovery animals'], True)
yy = bullets(yy + 5, 'Run after the primary results',
             ['Sequencing depth as a covariate', 'Second estimator: CLR linear model‡', 'Leave-one-farm-out; untargeted genus scan', 'LinDA and ALDEx2',
              'Freedman-Lane permutation; spike-in', 'Interval validation; null-scan calibration'], False)
ph = max(d_end, yy - LH + 6) - py
rect(M, py, dw, ph, fill='none', stroke=GRID, sw=0.7)
rect(ex, py, ew, ph, fill='none', stroke=GRID, sw=0.7)
H = py + ph + 4

bad = [(round(w, 1), round(m, 1), s) for w, m, s in checks if w > m + 0.01]
assert not bad, "text overflow (width, space, text): " + repr(bad)
svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="170mm" height="{H/W*170:.1f}mm" viewBox="0 0 {W} {H:.0f}" font-family="{FONT}">'
       '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
       f'<path d="M0,0 L10,5 L0,10 z" fill="{INK2}"/></marker></defs><rect width="{W}" height="{H:.0f}" fill="white"/>'
       + '\n'.join(out) + '</svg>')
root = sys.argv[1] if len(sys.argv) > 1 else '/data/rumen_multicohort'
od = os.path.join(root, 'results', 'figures_draft18'); os.makedirs(od, exist_ok=True); fo = os.path.join(od, 'fig1_study_design.svg')
if os.path.exists(fo): sys.exit(f'output exists, never overwritten: {fo}')
open(fo, 'w').write(svg)
print(f'height {H:.0f} pt = {H/W*170:.1f} mm; texts checked {len(checks)}; smallest text {S_SM} pt')
