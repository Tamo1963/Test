# Napier-Grass CSTR Design Suite

Streamlit application + calculation engine that sizes an industrial cascade of
CSTR anaerobic digesters for crushed **Napier grass** (Pakchong-1 lineage),
scaled from the lab-validated mixing study:

> P. Satjaritanun, Y. Khunatorn, N. Vorayos, S. Shimpalee, E. Bringley,
> "Numerical analysis of the mixing characteristic for napier grass in the
> continuous stirring tank reactor for biogas production",
> *Biomass and Bioenergy* 86 (2016) 53–64,
> doi:10.1016/j.biombioe.2016.01.018

## Run

```bash
pip install -r requirements.txt
streamlit run app.py
```

## What it does

- **Left panel (sidebar inputs):** feed rate (t/d), TS (%), VS (% of TS),
  particle density (kg/m³) + advanced set (HRT, module count, target slurry TS,
  CH₄ fraction, yield band, mixing power density).
- **Main workspace:** computed headline metrics and the complete **feed
  pre-treatment specification** (silage bunkers, chopping, hammer milling,
  dilution, pumping, floatation countermeasures).
- **Tabs:**
  1. *Reaction kinetics & assumptions* — first-order model, live conversion
     curve, full assumption list.
  2. *Outputs (in sequence)* — the entire calculation trace, step by step.
  3. *Recommendations* — engineering actions generated from the current inputs.
  4. *Conclusions* — summary metrics and closing inferences.
  5. *Design report & download* — one-click **markdown design document**
     (preamble, inputs table, equations, assumptions, step-by-step
     calculations, result tables, inferences & conclusions, references).

## Module layout

| File | Role |
|---|---|
| `app.py` | Streamlit front end |
| `design_core.py` | Parametric design engine + markdown report generator |
| upstream repo docs | `../napier-cstr-design-120tpd.md` (fixed 120 t/d recast), `../napier-grass-bioreactors.md` (installation register) |

Key defaults reproduce the fixed 120 t/day recast: 120 t/d feed, TS 15 %,
VS 85 % of TS, particle density 843 kg/m³, HRT 40 d, 4 modules, 12 % slurry TS,
Np = 2.00 dual-propeller + baffles mixing concept.
