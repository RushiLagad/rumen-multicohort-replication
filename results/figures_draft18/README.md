# Draft 18 figures, Table 4 data and Additional file 4

Built for draft 18 (Amendment 48.11; reporting only, no new fit). The draft 17 figures are in `results/figures_manuscript/`, unchanged.

| File | Manuscript | Script | Notes |
|---|---|---|---|
| fig1_study_design.{svg,png,pdf} | Fig. 1 | scripts/140_figure1_compact.py | Same content as 111, compact layout, 170 x 178.6 mm. PNG and PDF rendered from the SVG with rsvg-convert (`-d 300 -p 300 -b white`), the renderer that reproduces the draft 17 PNG exactly. |
| fig2_forest_pooled.{pdf,png} | Fig. 2 | scripts/139_figures_draft18.R | As 102 (Amendment 47 CLR, modified Knapp-Hartung), 170 x 150 mm (FIG2_H_MM). |
| fig4_clade_ratio.{pdf,png} | Fig. 4 | scripts/139_figures_draft18.R | Pooled row = final CLR modified Knapp-Hartung interval; the mixed-outcome pooled row is in Additional file 2. |
| fig5_dairy_effects.{pdf,png} | Fig. 5 | scripts/139_figures_draft18.R | From 48.3 (results/late_a48/dairy/); b, days-in-milk positive control. |
| af4_seven_claims.{pdf,png} | Additional file 4 | scripts/139_figures_draft18.R | 48.1 forest plot. |
| table4_seven_claims.tsv | Table 4 | scripts/139_figures_draft18.R | 48.1; permutation p reported as (b + 1)/(B + 1), B = 200. |

Additional files 1 to 4 for draft 18 are in `results/additional_files_draft18/` (scripts 141 and 142).
