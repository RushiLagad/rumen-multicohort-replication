#!/usr/bin/env python3
"""Figure 1 (manuscript): study design schematic. Writes fig1_study_design.svg.
Render to PDF/PNG with: node scripts/111_render_svg.js fig1_study_design.svg out.png out.pdf 160.2mm
Canvas 642 px wide = 170 mm, so 1 px = 0.75 pt. Smallest text 8 px = 6 pt."""
import re, math

W = 642
INK, INK2, MUTED, GRID = '#0b0b0b', '#52514e', '#898781', '#d6d4cc'
BLUE_D, BLUE, BLUE_P = '#104281', '#256abf', '#eaf1fb'
ORNG, ORNG_P = '#b8601d', '#fbf0e6'
FONT = "Arial, 'Liberation Sans', Helvetica, sans-serif"
S_SM, S_BODY, S_HEAD, S_TITLE, S_LET = 8.0, 8.6, 9.2, 10.6, 13.3   # px; x0.75 = pt

out = []
def add(s): out.append(s)
def esc(s): return s.replace('&', '&amp;').replace('<', '&lt;').replace('>', '&gt;')
def rich(s):
    parts = re.split(r'(\*[^*]+\*|_\{[^}]+\}|\^\{[^}]+\})', s)
    o = ''
    for p in parts:
        if not p: continue
        if p.startswith('*'): o += f'<tspan font-style="italic">{esc(p[1:-1])}</tspan>'
        elif p.startswith('_{'): o += f'<tspan baseline-shift="sub" font-size="72%">{esc(p[2:-1])}</tspan>'
        elif p.startswith('^{'): o += f'<tspan baseline-shift="super" font-size="72%">{esc(p[2:-1])}</tspan>'
        else: o += esc(p)
    return o
def text(x, y, s, size=S_BODY, weight='normal', color=INK, anchor='start', lh=None):
    lines = s if isinstance(s, list) else [s]
    lh = lh or size * 1.27
    for i, ln in enumerate(lines):
        if ln:
            add(f'<text x="{x:.1f}" y="{y + i*lh:.1f}" font-size="{size}" font-weight="{weight}" fill="{color}" text-anchor="{anchor}">{rich(ln)}</text>')
    return y + len(lines) * lh
def rect(x, y, w, h, fill='white', stroke=GRID, sw=0.8, r=3, dash=None):
    d = f' stroke-dasharray="{dash}"' if dash else ''
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="{w:.1f}" height="{h:.1f}" rx="{r}" fill="{fill}" stroke="{stroke}" stroke-width="{sw}"{d}/>')
def line(x1, y1, x2, y2, color=MUTED, sw=0.9, arrow=False):
    a = ' marker-end="url(#arr)"' if arrow else ''
    add(f'<line x1="{x1:.1f}" y1="{y1:.1f}" x2="{x2:.1f}" y2="{y2:.1f}" stroke="{color}" stroke-width="{sw}"{a}/>')
def panel(x, y, letter, title):
    text(x, y, letter, size=S_LET, weight='bold')
    text(x + 13, y, title, size=S_TITLE, weight='bold')

# ================================================================ a. cohorts
ax, ay = 6, 16
panel(ax, ay, 'a', 'Five public cohorts')
cw, ch, gap = 124, 96, 6
def card(x, y, name, n, place, acc, data, trait, design, col, pale):
    rect(x, y, cw, ch, fill=pale, stroke=col, sw=0.9)
    add(f'<rect x="{x:.1f}" y="{y:.1f}" width="3.2" height="{ch}" fill="{col}"/>')
    text(x + 9, y + 13, name, size=S_HEAD, weight='bold')
    text(x + cw - 6, y + 13, n, size=S_HEAD, weight='bold', color=col, anchor='end')
    text(x + 9, y + 25, [place, acc], size=S_SM, color=INK2, lh=10)
    text(x + 9, y + 50, [data, trait], size=S_BODY, color=INK, lh=11)
    text(x + 9, y + 78, ['Design terms:', design], size=S_SM, color=INK2, lh=10)
ly = ay + 17
text(ax, ly, 'Feed efficiency, beef cattle', size=S_HEAD, weight='bold', color=BLUE_D)
cy = ly + 6
card(ax, cy, 'Nelore', 'n = 52', 'Brazil', 'PRJNA987743', 'Shotgun metagenome', 'RFI, continuous', 'diet, slaughter batch', BLUE, BLUE_P)
card(ax + cw + gap, cy, 'Charolais', 'n = 30', 'France', 'PRJEB70895', 'Metatranscriptome', 'RFI, high vs low class', 'diet, year', BLUE, BLUE_P)
card(ax + 2*(cw + gap), cy, 'Alberta', 'n = 567', 'Canada', 'PRJNA393057', '16S rRNA gene', 'RFI, continuous', 'management group, breed', BLUE, BLUE_P)
ly2 = cy + ch + 16
text(ax, ly2, 'Methane, dairy cows', size=S_HEAD, weight='bold', color=ORNG)
cy2 = ly2 + 6
card(ax, cy2, 'Difford', 'n = 749', 'Denmark', 'PRJEB28065', '16S rRNA gene', 'CH₄, g per day', 'herd, batch, parity, DIM', ORNG, ORNG_P)
card(ax + cw + gap, cy2, 'Wallace', 'n = 1,009', 'Europe, 7 farms', 'PRJNA517480', '16S rRNA gene', 'CH₄, g per day', 'farm, parity, DIM', ORNG, ORNG_P)
nx = ax + 2*(cw + gap) + 6
text(nx, cy2 + 14, ['About 2,400 animals.', '', 'All sequence data are', 'public; phenotypes come', 'from the original', 'publications.', '', 'DIM, days in milk'], size=S_BODY, color=INK2, lh=11)
a_bottom = cy2 + ch

# ================================================================ b. targets
bx = ax + 3*cw + 2*gap + 12
bw = W - 6 - bx
panel(bx, ay, 'b', 'Published targets')
def tset(y, col, pale, head, rows, result):
    h = 18 + 11.5*len(rows) + 16
    rect(bx, y, bw, h, fill='white', stroke=col, sw=0.9)
    text(bx + 7, y + 13, head, size=S_HEAD, weight='bold', color=col)
    yy = y + 27
    for code, desc in rows:
        if code:
            add(f'<rect x="{bx+7:.1f}" y="{yy-8:.1f}" width="11" height="10" rx="1.5" fill="{pale}" stroke="{col}" stroke-width="0.6"/>')
            text(bx + 12.5, yy, code, size=S_SM, weight='bold', color=col, anchor='middle')
        text(bx + 22, yy, desc, size=S_SM, color=INK)
        yy += 11.5
    text(bx + 7, yy + 4, result, size=S_BODY, weight='bold', color=INK)
    return y + h
y1 = tset(cy, BLUE, BLUE_P, 'Feed efficiency: 76 targets',
          [('T', 'Li et al. 2019: 6 taxa'), ('N', 'Nelore discovery study: 64 genera'), ('C', 'Canonical rumen taxa: 6')],
          '27 testable')
y2 = tset(y1 + 7, ORNG, ORNG_P, 'Methane: 6 targets + Nelore genera',
          [('M', 'Reported in >1 country: 3'), ('', '*M. gottschalkii* : *M. ruminantium*'), ('', 'clade ratio; Succinivibrionaceae;'), ('', 'total archaea'),
           ('S', 'Secondary: 3 genera')],
          '11 testable')
text(bx, y2 + 14, ['Testable: mapped feature in at', "least 10% of a cohort's animals.", 'A cohort that reported a target,', 'or helped choose it, counts as', 're-analysis, not replication.'], size=S_SM, color=INK2, lh=10)

# ================================================================ c. estimation
cy0 = a_bottom + 26
panel(ax, cy0, 'c', 'Estimation in each cohort, for each target')
sy = cy0 + 9
ns = 4
sw_ = (W - 12 - (ns - 1)*13) / ns
sh = 132
steps = [
    ('1  Profiling', [('b', 'Shotgun, metatranscriptome'), ('', 'Host reads removed'), ('', '(minimap2, ARS-UCD2.0)'), ('', 'Kraken 2 + Bracken, GTDB r207'),
                      ('b', '16S rRNA gene'), ('', 'cutadapt, DADA2, SILVA 138.1'), ('', ''), ('', 'Genus table: genera in'), ('', '≥10% of animals')]),
    ('2  Detection check†', [('', 'Same animals, table and'), ('', 'model, with the design'), ('', 'factor as the tested term'), ('', ''),
                                ('', 'Null results count only if'), ('', '≥1 genus reaches q < 0.05'), ('', ''), ('', 'Targets mapped by name;'), ('', 'family targets summed')]),
    ('3  ANCOM-BC2', [('', 'genus ~ design terms + trait'), ('', ''), ('', '*W* = log fold change / SE'), ('', ''),
                      ('', 'Continuous traits z-scored'), ('', 'within cohort, and within'), ('', 'farm for Wallace')]),
    ('4  Permutation calibration', [('', 'Trait shuffled 200 times'), ('', 'within design groups'), ('', ''), ('', ''), ('', ''), ('', ''), ('', ''),
                                    ('', '*W*_{cal} = *W* / max(1, SD_{null})'), ('', '*r* = *W*_{cal} / √(*W*_{cal}^{2} + df)'), ('', 'plus permutation *p*')]),
]
for i, (head, rows) in enumerate(steps):
    x = 6 + i*(sw_ + 13)
    rect(x, sy, sw_, sh, fill='white', stroke=INK2, sw=0.8)
    add(f'<path d="M{x:.1f},{sy+16:.1f} v-13 a3,3 0 0 1 3,-3 h{sw_-6:.1f} a3,3 0 0 1 3,3 v13 z" fill="{INK2}"/>')
    text(x + 7, sy + 11.5, head, size=S_HEAD, weight='bold', color='white')
    yy = sy + 29
    for kind, s in rows:
        text(x + 7, yy, s, size=S_SM, weight='bold' if kind == 'b' else 'normal', color=INK if kind == 'b' else INK2)
        yy += 10.2
    if i < ns - 1:
        line(x + sw_ + 1.5, sy + sh/2, x + sw_ + 12, sy + sh/2, color=INK2, sw=1.1, arrow=True)
# null glyph in step 4
gx = 6 + 3*(sw_ + 13) + 10; gy = sy + 86; gw = sw_ - 20
def curve(sd, col, dash=None):
    pts = [(gx + gw*k/60, gy - 24*math.exp(-0.5*((-4 + 8*k/60)/sd)**2)/sd) for k in range(61)]
    d = 'M' + ' L'.join(f'{a:.1f},{b:.1f}' for a, b in pts)
    da = f' stroke-dasharray="{dash}"' if dash else ''
    add(f'<path d="{d}" fill="none" stroke="{col}" stroke-width="1.2"{da}/>')
line(gx, gy, gx + gw, gy, color=GRID, sw=0.8)
curve(1.0, MUTED, '2.2,1.6'); curve(2.0, BLUE)
text(gx + gw, gy - 25, 'nominal', size=S_SM, color=MUTED, anchor='end')
text(gx + gw, gy - 14, 'observed', size=S_SM, color=BLUE, anchor='end')
c_bottom = sy + sh

# ================================================================ d. pooling
dy0 = c_bottom + 26
panel(ax, dy0, 'd', 'Pooling and replication label')
dw = 330; py = dy0 + 9; ph = 142
rect(6, py, dw, ph, fill='white', stroke=GRID, sw=0.8)
text(14, py + 15, 'Fisher *z*  →  REML random effects  →  BH within target set', size=S_BODY, color=INK)
crit = [['Pooled', '*q* < 0.05'], ['Same sign as', 'published'], ['Holds with each', 'cohort left out']]
bxw, bxh, cx0, cyy = 74, 32, 14, py + 26
for k, c in enumerate(crit):
    x = cx0 + k*(bxw + 13)
    rect(x, cyy, bxw, bxh, fill=BLUE_P, stroke=BLUE, sw=0.8)
    text(x + bxw/2, cyy + 13.5, c, size=S_SM, color=INK, anchor='middle', lh=10)
    line(x + bxw + 1.5, cyy + bxh/2, x + bxw + 12, cyy + bxh/2, color=INK2, sw=1, arrow=True)
rx = cx0 + 3*(bxw + 13)
rect(rx, cyy + 5, 6 + dw - rx - 8, bxh - 10, fill=BLUE_D, stroke=BLUE_D)
text(rx + (6 + dw - rx - 8)/2, cyy + bxh/2 + 3, 'Replicated', size=S_BODY, weight='bold', color='white', anchor='middle')
text(14, cyy + bxh + 15, ['Otherwise: single-cohort support, opposite sign,', 'not replicated, or untestable'], size=S_SM, color=INK2, lh=10)
text(14, cyy + bxh + 44, 'Power and equivalence', size=S_BODY, weight='bold')
text(14, cyy + bxh + 56, ['Effective *n*, smallest *r* detectable at 80% power, and', 'whether the 95% CI excludes *r* = 0.15 or 0.30 in the', 'published direction'], size=S_SM, color=INK2, lh=10)

# ================================================================ e. robustness
ex = 6 + dw + 12; ew = W - 6 - ex
panel(ex, dy0, 'e', 'Robustness checks')
rect(ex, py, ew, ph, fill='white', stroke=GRID, sw=0.8)
def bullets(y, head, items, filled):
    text(ex + 9, y, head, size=S_BODY, weight='bold')
    y += 12.5
    for s in items:
        if filled: add(f'<circle cx="{ex+12:.1f}" cy="{y-2.8:.1f}" r="2.1" fill="{INK2}"/>')
        else: add(f'<circle cx="{ex+12:.1f}" cy="{y-2.8:.1f}" r="2.1" fill="white" stroke="{INK2}" stroke-width="0.9"/>')
        text(ex + 19, y, s, size=S_SM, color=INK2); y += 10.4
    return y
yy = bullets(py + 15, 'Specified before the affected results',
             ['Wallace without its Swedish farm', 'Methane per kg dry matter intake', 'Wallace as a fourth RFI cohort', 'Alberta without possible discovery animals'], True)
bullets(yy + 6, 'Run after the primary results',
        ['Sequencing depth as a covariate', 'Second estimator: CLR linear model‡', 'Leave-one-farm-out; untargeted genus scan', 'LinDA and ALDEx2'], False)
d_bottom = py + ph

H = d_bottom + 8

svg = (f'<svg xmlns="http://www.w3.org/2000/svg" width="170mm" height="{H/W*170:.1f}mm" viewBox="0 0 {W} {H:.0f}" font-family="{FONT}">'
       '<defs><marker id="arr" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="5" markerHeight="5" orient="auto-start-reverse">'
       f'<path d="M0,0 L10,5 L0,10 z" fill="{INK2}"/></marker></defs><rect width="{W}" height="{H:.0f}" fill="white"/>'
       + '\n'.join(out) + '</svg>')
open('fig1_study_design.svg', 'w').write(svg)
print(f'height {H:.0f} px = {H/W*170:.0f} mm')
