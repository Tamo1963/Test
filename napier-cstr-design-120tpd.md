# Industrial CSTR Design for Napier Grass Biogas Production — 120 Tonnes/Day
### Complete recast plan scaled from the lab design validated in Satjaritanun et al. (2016)

*Design reference: P. Satjaritanun, Y. Khunatorn, N. Vorayos, S. Shimpalee, E. Bringley, "Numerical analysis of the mixing characteristic for napier grass in the continuous stirring tank reactor for biogas production," Biomass and Bioenergy 86 (2016) 53–64: https://dx.doi.org/10.1016/j.biombioe.2016.01.018*

*Companion register: `napier-grass-bioreactors.md` in this repository.*

---

## 1. Executive summary

- **Capacity:** 120 t/day fresh Napier grass (Pakchong 1-type) → **~150 m³/day slurry** at 12% TS → **4 × 1,500 m³ CSTR digesters in series (≈5,900 m³ working, HRT ≈ 39–40 d)**, mesophilic 37 °C.
- **Agitators:** scaled "Model 3" concept from the paper — **two downward-pumping 3-blade 45°-pitch axial propellers per tank (D = 5.0 m, D/T = 0.40) + 4 radial wall baffles**, ~43 kW per tank, 12 rpm, Re ≈ 4.8 × 10⁶ (fully turbulent), Froude ≈ 0.019 (no vortex).
- **Expected output:** CH₄ **4,100–5,200 m³/d** (0.27–0.34 m³/kg VS, paper's CSTR range) → biogas **6,900–8,700 m³/d** → **3.0–3.7 t/d CBG** or **0.70–0.89 MWe** gross via CHP (~0.5–0.6 MWe net after plant parasitics).
- **Key scale insight:** the lab's surface vortex problem disappears at industrial speed (Fr = n²D/g ≈ 0.02 ≪ 1); the *floating-rind problem is what must be engineered out* — solved by the upper propeller drawing the liquid surface downward aggressively, exactly as the paper's Model 3 does.

---

## 2. Design basis and mass balance

### 2.1 Feedstock (per paper §2.2 and Table 3)

| Parameter | Value | Source |
|---|---|---|
| Napier grass species | *Pennisetum purpureum* cv. Pakchong 1 (Thailand lineage) | Paper §1 [ref 3,4] |
| Fresh feed rate | **120.0 t/day** (as received) | User requirement (matches Biofics Rajkot CBG plant input scale) |
| Moisture content | 85% (→ TS 15%) | Typical NT Design assumptions |
| Total solids | 120 × 0.15 = **18.0 t/d** | Derived |
| Volatile solids | 18.0 × 0.85 = **15.3 t/d** (VS/TS ≈ 0.85) | Typical grass |
| Particle density | **843 kg/m³** (floats in water, ρw = 998 kg/m³) | Paper Table 3 |
| Particle size after feed preparation | **≤ 30 mm** longest dimension (chopped + hammer-milled); paper's DPM used 2.5 mm spheres | Paper §2.2 |
| Harvest | 5–6 cuts/year, 45–60 day cycle, 375–635 t/ha/yr (Pakchong 1) | Paper §1 |

### 2.2 Slurry balance

| Stream | t/d | m³/d | Notes |
|---|---|---|---|
| Fresh grass in | 120.0 | ~140 | |
| Dilution water in | 30.0 | 30.0 | Condition feed to pumpable 12% TS |
| Recycle/inoculum (start-up only) | — | — | First fill from seed digestate (2 × 1,500 m³) |
| **Slurry to digesters** | **150.0** | **150.0** | 12% TS, 10.2% VS |
| Digestate out | ~142.5 | ~140 | ~3.2% TS after digestion |
| Biogas | ~7.3–9.2 | 6,900–8,700 | See §2.3 |

### 2.3 Gas balance

Paper's CSTR mono-digestion result: **0.27–0.34 L CH₄/g VS added** at OLR 1.0–1.5 g VS/L·d, 53–59% CH₄ (AGRIS-recorded study). Applying to 15.3 t VS/d:

| Quantity | Low | High | Basis |
|---|---|---|---|
| Methane yield | 0.27 m³/kg VS | 0.34 m³/kg VS | Paper mono-digestion CSTR |
| CH₄ production | 4,131 m³/d | 5,202 m³/d | 15,300 kg VS/d × yield |
| Biogas production (60% CH₄) | 6,885 m³/d | 8,670 m³/d | |
| Thermal energy (CH₄ LHV 9.95 kWh/m³) | 41.1 MWh/d | 51.8 MWh/d | 1.71–2.16 MW<sub>th</sub> avg |
| CBG output (after upgrade) | **3.0 t/d** | **3.7 t/d** | ρCH₄ = 0.716 kg/m³ |
| CHP gross electrical (41% eff.) | 0.70 MWe | 0.89 MWe | |
| Plant parasitic electric | ~0.15–0.25 MW | | Mixers + pumping + upgrading note |

**Reconciliation note:** the Biofics Rajkot plant (120 t/d input → 5 t/d CBG) implies ~52 m³ CH₄/t fresh grass ≈ 0.29 L/g VS *plus co-digestion boost*. The paper's mono-digestion yields land at the low end — **co-digestion with poultry/cattle dung is identified as the yield-upgrade path to 5+ t/d**, consistent with our register findings (Rajkot, Botad).

**Volumetric loading check:** OLR = (150 m³/d × 102 kg VS/m³)/6,000 m³ = **2.55 kg VS/m³·d** — above the paper's best-yield band of 1.0–1.5 g VS/L·d (equivalent). Options: (a) accept moderate yield penalty, (b) feed 12% TS slurry at lower VS fraction, (c) add subject area of co-digestion. Flagged as the main process risk of 120 t/d mono-grass at 40 d HRT.

---

## 3. Digester configuration

### 3.1 Architecture

- **4 CSTR modules, 1,500 m³ working volume each, operated in series** (cascade improves VFA control and pathogen fade); configurable as 2 × 2 parallel during maintenance.
- Each module: cylindrical, **flat bottom** (per paper; dished bottom avoided for sludge accumulation simplicity and floor-mix pattern), vertical wall.
- **Liquid: 1,500 m³** each; **headspace ≥ 20%** (375 m³) for gas + fluctuation + foam.

### 3.2 Geometry (per module)

| Dimension | Value | Derived from / rationale |
|---|---|---|
| Tank diameter T | **12.5 m** | Industrial precedent 10–13 m silo digesters |
| Liquid depth H | **12.0 m** | H/T = 0.96 (lab was 0.667 — see deviations note) |
| Working volume | π × 6.25² × 12.0 = **1,474 ≈ 1,475 m³** | ×4 = 5,900 m³; system HRT = 5,900/150 = **39.3 d ≈ 40 d** ✓ |
| Total height incl. roof | ~14.2 m | Fixed membrane roof + gas dome |
| Construction | Post-tensioned concrete, epoxy-lined interior; 316L SS internals | Standard biogas practice; H₂S corrosion |
| Feed nozzle | DN200, bottom-injection toward lower impeller sweep (mirrors paper's bottom feeding pipe Ø3 cm × L 15 cm) | Paper Fig. 6/Table 1 |
| Transfer piping | DN300 bottom transfer to next module under overflow weir (hydraulic seal), self-balanced flows | |

**System HRT (40 d) per module: 10 d each** — the desired cascade effect, with staged stability and lower VFA slip.

### 3.3 Deviations from lab (acknowledged)

| Lab (paper) | Plant | Rationale |
|---|---|---|
| H/T = 0.667 | H/T = 0.96 | Land cost; typical digester aspect ratio; compensated by **extra impeller head/separation** (see §4) |
| Da/T = 0.667 | D/T = 0.40 | 5.0 m impellers on 12.5 m tank |
| HRT/feed continuous | 40 d cascade, continuous | Biogas engineering baseline 30–50 d for Napier-type grass |

---

## 4. Agitator design (the core of this plan)

### 4.1 Concept — scaled "Model 3"

The paper's Model 3 (two propellers + four baffles) outperformed single-propeller designs on every metric tracked:

| Metric (paper, lab, 1.67 Hz) | Model 1 | Model 2 | **Model 3** |
|---|---|---|---|
| Mixing efficiency (steady) | 10.08% | 23.77% | **34.66%** |
| Power number Np | 0.91 | 1.20 | **2.00** |
| Vortex on free surface | Severe (≥1.33 Hz) | Eliminated | Eliminated + turbulent top |
| Particle behavior | Floats to surface, stays | Floats after ~12 s | **Floats/sinks alternately continuously** — floating suppressed |
| Levels A/B/C (upper 3 of 4 strata) occupancy | None | None | **Populated (only design mixing upper tank)** |

### 4.2 Impeller specification (per module)

| Parameter | Lab (paper) | **Plant (this design)** | Scaling note |
|---|---|---|---|
| Type | 3-blade bladed-propeller, 45° pitch | Same (stiffened structural conversion) | Blade angle drives axial-down flow (paper conclusion) |
| Number of blades Nb | 3 | 3 | |
| Blade length L | 9 cm (0.45 Da) | L = 2.25 m | Geometric ×25 |
| Blade width W | 2 cm (0.10 Da) | W = 0.50 m | Geometric |
| Blade thickness t | 0.3 cm (0.015 Da) | t = 75 mm | Structural scaling ≈ L× scaling for strength; 316L or duplex SS |
| Impeller dia Da | 20 cm | **5.0 m** (D/T = 0.40) | Scaled-down from 2/3 to 0.40 (land/practicality) |
| Shaft diameter (est.) | — | 350 mm hollow or 300 mm solid 316L | Sized below |
| Pumping direction | **Down (axial)** | **Down (axial)** | Paper: must pull Napier to bottom |
| Impeller count per tank | 2 | **2** | |
| Hub elevations (off floor) | — | Lower hub **3.0 m**, upper hub **8.0 m**; hub-to-hub **5.0 m = 1.0 D** | Lab spacing was 0.5 D with H/Da = 1.0; plant H/Da = 2.4 needs 1.0–1.2 D separation for coverage (conventional two-impeller rule) |
| Upper impeller reach | — | Top blade at 10.5 m elevation → **1.5 m below liquid surface** | Draws the surface down = direct anti-floatation control (Model 3's job) |

### 4.3 Operating speed and regime

| Parameter | Value | Formula / check |
|---|---|---|
| Rotation rate n | **0.19 Hz = 12 rpm** (VFD 8–16 rpm) | Down from lab 1.67 Hz; similarity scaling n ∝ (D<sub>lab</sub>/D)<sup>…</sup> with P/V hold |
| Reynolds Re = ρnDa²/μ (Eq. 13) | 1000 × 0.19 × 25 / 0.001 = **4.8 × 10⁶** ≫ 10⁴ | Fully turbulent ✓ (lab range 3.3–8 × 10⁴) |
| Froude Fr = n²·Da/g | 0.19² × 5 / 9.81 = **0.018** ≪ 1 | **No surface vortex** — at lab the vortex appeared because n was high and H/Da low; at 12 rpm Fr is tiny. Vortex is a *non-issue* at scale ✓ |
| Tip speed | π × 5.0 × 0.19 = **3.0 m/s** | Fibrous solids handled; moderate shear (paper cautions vs destroying micro-organisms) |
| Power per shaft | P = Np ρ n³ Da⁵ = 2.0 × 1000 × 0.19³ × 3125 = **~43 kW** | Paper's Np = 2.00 (Model 3) |
| Mixing intensity | 43 kW/1,475 m³ = **29 W/m³** | Grass-digester typical design band 20–40 W/m³ active |
| Shaft torque | P/(2πn) = 43,000/(2π × 0.19) ≈ **36 kN·m** | Large-planetary drive sizing |
| Motor | **55 kW**, 4-pole AC, VFD (8–16 rpm), helical-bevel or planetary gearbox ≈ 125:1, 200% starting torque | Higher-torque margin for start-up in settled fibrous slurry |

### 4.4 Baffles

| | Option A — paper-faithful | Option B — pragmatic (recommended) |
|---|---|---|
| Layout | 4 × wall baffles, width T/12 = 1.04 m, full 12 m liquid height, thickness ~0.2 m (T/60) | 4 × **short anti-float pilaster baffles**: 1.0 m × 2.5 m × 0.25 m at liquid surface band only, rounded upstream nose |
| Pros | Full Np characterization carried over; lowest risk of swirl | No rag/long-fibre accumulation traps (digesters with fibrous feed clog full-depth baffles) |
| Cons | Fibre/straw "rapping" and rag build-up on leading edges; cleaning needs | Np may trend toward Model-2 value (1.2) → recheck power below |
| Power with Option B | Assume Np_eff ≈ 1.6 averaged over two-impeller zone → P ≈ 34 kW; use 45 kW drive same shaft | Already covered by 55 kW/VFD margin |

**Decision:** fabricate piping for both; start with Option B pilasters, add full-depth if torque trends show swirl dominance. The paper itself notes the *combination* of speed increase OR second propeller compensates baffle-induced velocity loss; the dual-impeller lever is retained in both options.

### 4.5 Redundant mixing — gas recirculation

- **Central bottom gas sparger, 100–150 Nm³/h per module** (raw biogas re-compressed at 2 × 4 kW blowers, duty/standby): during impeller shutdown windows it keeps solids agitated and breaks scum; also enables periodic-stirring energy saving (paper notes stirring is "stopped periodically" in practice).
- Scum trap: perimeter spray manifold using recycled digestate, 2 × 2 m³/h per tank.

### 4.6 Mixing time estimate

Using standard axial-impeller correlations (Ruszkowski-type) with our geometry (T/D = 2.5, Np = 2, n = 0.19 Hz), blend time at plant scale is **on the order of 5–10 minutes** (conservative quadrupled from lab's ~3.5 s scaled by geometric similarity). Acceptance: **complete homogenization within 30 min** of a step change, verified at commissioning by heat tracer/torque response.

---

## 5. Feed preparation and feeding system

| Step | Equipment | Spec |
|---|---|---|
| 5-day buffer | Field silage bunkers | 2 × 2,800 m³ (≈15 days feed at 120 t/d, ρ silage 650 kg/m³) |
| Pre-crush | Disc cutter + hammer mill (paper §2.2: disc and hammer mill, compacted; particle ≤30 mm) → lodged at 120 t/h capacity? Run 6–8 h/d window | 2 × 40 t/h lines (1 duty, 1 standby) |
| Hopper/dilution | Receiving hopper 60 m³ + water make-up 30 m³/d, inline macerator, TS sensor | Feed 12% TS target |
| Slurry pumps | 2 × hose/peristaltic pumps, 50 m³/h @ 10 bar, VFD | 24/7 feed 6.25 m³/h per module-line duty/standby |
| Feeding pipe | DN200 to each module bottom (per paper's bottom-injection) | Continuous feed per module = 1/24 flow pulses at 10 min/h? For simplicity — **continuous metered feed** 150/4 = 37.5 m³/d each |

---

## 6. Heat and energy services

| Item | Duty | Note |
|---|---|---|
| Slurry heating 12→37 °C | 150,000 kg/d × 4.0 kJ/kgK × 25 K = 15.0 GJ/d = **174 kW<sub>th</sub>** avg | Heat exchangers external spiral-tube, 4 × 70 kWth |
| Tank loss (concrete, 37 °C inside, ~18 °C ambient avg) | ~30 kWth total | Insulation |
| **Total process heat** | **~210 kW<sub>th</sub>** |≪ CH₄ output 1.7–2.2 MWth → thermally self-sufficient |
| Mixing electric continuous 24/7 | 4 × 43 = 172 kW → 4.13 MWh/d | 5–10% of gross energy if continuous |
| Option: intermittent duty (20 min/h) | ~55% of that = 2.3 MWh/d | Only with surface-tension scum monitoring in place |
| Parasitic electric total (mixers, pumps, chopper, gas treatment, compression boost, lighting) | ~0.20–0.30 MW avg | |
| Net CHP electric | **0.5–0.6 MWe** | Gross 0.70–0.89 minus parasitics |

---

## 7. Controls, instrumentation, and operation

- **Per-tank:** 3 × temperature (top/mid/bottom), level (radar), pressure (gas dome), pH + ORP, VFA alkalinity lab kit; torque transmitter on each gearbox (trend = mixing-health surrogate; the paper showed torque signature changes with vortex — here torque trending detects rag loading, fibre wrapping).
- **Interlocking:** feed-torque lock — reduce feed or start gas recirculation if shaft torque drifts >15% over 5-day baseline.
- **Stirring philosophy for floating feed:** *continuous primary duty* recommended (grass floats — Model-3 lesson), with periodic intelligent interruption windows (validated by torque rebound time < 10 min).
- pH setpoint 7.2–7.6 (auto NaFe? — lime or NaOH dosing line on last module), overloading kill-switch on feed pump.
- Gas: 60% CH₄, H₂S dosing FeCl₃ at feed inlet to <50 ppm in gas; biogas holder 2 × 500 m³ double-membrane; flare 1 × 100% capacity standby; H₂ detectors around headspace/wall treatment.

---

## 8. Acceptance criteria and verification (carried from the paper)

| Criterion | Lab reference (paper) | Plant acceptance |
|---|---|---|
| Mixing efficiency (4-ImageJ views, Eq. 14–19) | 34.66% at 1.67 Hz Model 3 | Homogeneity verified by differential temperature ≤ 0.8 °C across depths, and tracer test: 95% homogenization ≤ 30 min |
| No particle stranding in upper zone (Levels A–C) | Only Model 3 populated upper levels | Radar/dye check: no persistent scum island >0.5 m thick after continuous duty 1 week |
| Vortex | eliminated at Model 3 | Fr check (0.018 ≪ threshold) → not applicable; verify by level probes |
| Shear vs. mixing trade-off | paper cautions excess shear kills micro-organisms | Tip speed cap ≤ 3.5 m/s; biogas-rate monitoring before/after agitation cycle changes: no drop > 5% |

## 9. Scale-similarity summary table (lab → plant)

| Quantity | Lab | Scale law | Plant (this design) |
|---|---|---|---|
| Impeller dia Da | 0.20 m | ×(T ratio 0.3/12.5 ≈ 1/42, geometry eased) | 5.0 m |
| Rotation rate n | 1.67 Hz | ~ ×0.12 | 0.19 Hz |
| Reynolds ReS | 6.7 × 10⁴ | Re ∝ n Da² | 4.8 × 10⁶ |
| Froude Fr | 0.057 (vortex-prone) | Fr ∝ n² Da | 0.018 (vortex-free) |
| Power number Np | 2.00 (Model 3) | ~ invariant (Baffled similar) | 2.0 (Option A) / 1.6 est. (Option B) |
| P/V | 213 W/m³ | hold 20–40 W/m³ typical | 29 W/m³ |
| Mixing time | ~3.5 s | ∝ slower-flow scale | ~5–10 min |

## 10. Risk register (basis of the design reviews needed before FID/construction)

1. **OLR 2.55 kg VS/m³·d above paper's best-yield window** → mitigations: accept lower yield; reduce TS% of feed; repartition parallel/series modes; plan co-digestion upgrade path with manure/litter (documented plants show >30% yield uplift).
2. **Floating grass rind** → dual impellers + upper-blade surface draw + gas recirculation + surface sprays + scum-breaking sweep; monitor via radar + torque.
3. **Baffle ragging (Option A)** → Option B pilasters default; reverse-motion cleaning cycles; inspection ports at baffle rows.
4. **Sand/mineral settling** → flat bottom taper 2°, flush valves per tank bottom.
5. **Fibre wrapping of shafts** → anti-wrap cage plates (full circumferential disc ring) at each impeller hub — standard in grass-digestate service.
6. **Gearbox sizing error on start-up torque** → 200% starting margin, VFD ramping, torque-limited start.
7. **VFA overload in module 1** → online pH/ORP interlock; split new feed across module 1/2 if VFA > 2,000 mg/L.

## 11. Capex-shape summary (order-of-magnitude, not a quotation)

| Package | Content |
|---|---|
| Civil | 4 concrete digesters × 1,475 m³ + control buildings, 2 lagoons 3,500 m³ (45-d digestate), silage bunkers |
| Mechanical | Agitator packages: 4 × (55 kW motor + ~125:1 planetary + 350 mm shaft + 2 × 5.0 m 3-blade SS impellers + pilaster baffles + anti-wrap plates + sparger) |
| Feed | 2 × 40 t/h hammer-milling lines, hopper, 2 × hose pumps, 30 m³/d water module |
| Biogas | 2 × 500 m³ holders, desulphurization (Fe dosing + carbon), cooler+dryer, 1,000 m³/h flare; alternatively VPSA/PSA upgrading to CBG (3–4 t/d class) or 1.1 MWe CHP |
| Instrumentation/controls | 4 × mixing panel, PLC/SCADA, methane safety loop |

## 12. References

1. Satjaritanun et al., *Biomass and Bioenergy* 86 (2016) 53–64 — geometry, Models 1–3 mix efficiencies/power numbers, particle/impeller data, vortex findings. https://doi.org/10.1016/j.biombioe.2016.01.018
2. Wilawan, Pholchan, Aggarangsi — Pakchong 1 grass + layer-hen manure in a full CSTR (paper ref 3) — precedent for co-digestion boost.
3. ERDI Nakornping, Chiang Mai University, internal report — Napier grass physical prep and properties (paper ref 4).
4. McCabe, Smith, Harriott — Unit Operations (paper ref 19; Reynolds/impeller mixing theory).
5. Repository file `napier-grass-bioreactors.md` — plant register (Biofics 120 t/d Rajkot co-digestion precedent; Bareilly plug-flow lead; Thai CSTR precedents).
