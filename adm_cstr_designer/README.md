# ADM1 Anaerobic CSTR Designer

Standalone Streamlit tool: from a substrate (flow, temperature, TS, VS, feedstock blend) it solves an
ADM1-based anaerobic CSTR to steady state, sizes the agitators and heating coil, shows the kinetics,
builds a 3D model (FreeCAD script + in-app preview) and exports a Word design report.

```bash
pip install -r requirements.txt
streamlit run app.py          # run from this folder
pytest -q                     # tests (run from this folder with PYTHONPATH=.)
```

## Layout
| Path | Purpose |
|---|---|
| `core/adm1.py` | ADM1 (BSM2 form): Petersen matrix with C/N closure, rates, algebraic pH, gas transfer |
| `core/cstr.py` | Reactor + headspace ODE, steady-state solve, COD balance |
| `core/feed.py` | TS/VS/blend -> ADM1 influent (BMP-consistent biodegradable fraction) |
| `core/mixing.py`, `core/thermal.py`, `core/geometry.py` | Agitator, heat load / coil, tank geometry |
| `core/scenarios.py` | HRT / temperature / TS / blend sweeps |
| `data/kinetics_db.json`, `data/feedstocks.json` | Editable libraries |
| `cad/` | FreeCAD build script (`freecad_body.py`), generator, preview mesh |
| `report/docx_report.py` | Word report |
| `tests/` | BSM2 validation, sizing identities, UI smoke test |

## Validation
`tests/test_bsm2_validation.py` reproduces the BSM2 benchmark digester (Q = 178.5 m3/d, V = 3400 m3, 35 C)
to within 1 % on all states (pH 7.26) with COD balance closure within 1 %. Reference values are the MATLAB
steady state quoted in the QSDsan test-suite.

## FreeCAD
Tab 6 offers a self-contained script with this design's dimensions embedded. Run it with
`freecadcmd build_cstr_freecad.py`, or paste it into FreeCAD / a FreeCAD MCP execute-code tool. It writes
`.FCStd`, `.step` and `.stl`. The in-app viewer is an independent preview built from the same dimensions.
**The FreeCAD script has only been checked against stub modules, not run inside FreeCAD.**

## Limitations
- Feedstock data (TS, VS, COD/VS, BMP, composition, hydrolysis constants, ammonia) are indicative
  typical-range defaults, not cited measurements. Replace them with lab data and a BMP test.
- ADM1 constants are the BSM2 mesophilic set; temperature correction outside ~30-42 C is an extrapolation.
- ADM1 does not wash out methanogens until HRT ~ 3 d, so short-HRT results are optimistic.
- Viscosity, wall/coil U-values and impeller power/flow numbers are indicative; confirm with vendors.
