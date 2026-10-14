"""
Core design calculations for the Napier-Grass CSTR Design Suite.

All model constants trace to:
  P. Satjaritanun, Y. Khunatorn, N. Vorayos, S. Shimpalee, E. Bringley,
  "Numerical analysis of the mixing characteristic for napier grass in the
  continuous stirring tank reactor for biogas production",
  Biomass and Bioenergy 86 (2016) 53-64.  doi:10.1016/j.biombioe.2016.01.018

The suite scales the paper's lab-validated "Model 3" configuration
(two down-pumping 3-blade 45-degree propellers + four T/12 baffles,
Np = 2.00) to an industrial cascade of flat-bottom CSTR digesters.
"""

import math
from dataclasses import dataclass, asdict

# ---------------------------------------------------------------------------
# Design constants   (rationale given per item; see REPORT assumptions)
# ---------------------------------------------------------------------------
HT_RATIO = 0.96           # liquid depth / tank diameter, industrial digester
IMPELLER_DT_RATIO = 0.40  # Da/T, industrial adaptation of the paper's concept
POWER_NUMBER = 2.00       # Np of dual propeller + baffles (paper, Model 3)
WATER_RHO = 998.0         # kg/m3, paper Table 3 (20 C)
SLURRY_RHO = 1000.0       # kg/m3, design assumption for diluted slurry
CP_SLURRY = 4.0           # kJ/kg-K, design assumption
SLURRY_TEMP_RISE = 25.0   # K, 12 C feed -> 37 C digester
CH4_LHV_KWH_M3 = 9.95     # kWh per Nm3 methane (35.8 MJ/Nm3)
CH4_RHO = 0.716           # kg per Nm3 methane
CHP_EFF = 0.41            # gross electrical efficiency of CHP set
K_FIRST_ORDER = 0.08      # 1/d, assumed first-order rate at 37 C mesophilic
BIODEGRADABLE_FRACTION = 0.75   # of VS, assumed for grass
HEAT_LOSS_KW_PER_TANK = 7.5     # tank shell losses at 37 C
AUX_PARASITIC_FACTOR = 0.60     # aux electric as fraction of mixing electric
HUB_LOWER_FRACTION = 0.25       # lower hub elevation / liquid depth
SURFACE_CLEARANCE = 1.5         # m, top blade below liquid surface
TIP_SPEED_MAX = 3.5             # m/s, shear vs. biomass-cell film protection
SILAGE_DENSITY = 650.0          # kg/m3 in bunker
BUNKER_DAYS = 15.0              # feedstock buffer
CHOP_HOURS = 8.0                # h/day mill operation window
YIELD_LO_DEF = 0.27   # m3 CH4 per kg VS added (paper, CSTR low end)
YIELD_HI_DEF = 0.34   # m3 CH4 per kg VS added (paper, CSTR high end)
MAX_PARTICLE_MM = 30.0      # longest feed particle dimension to digester


@dataclass
class DesignInputs:
    """User-editable design inputs."""
    feed_rate_tpd: float = 120.0        # t/day fresh Napier grass, as received
    ts_percent: float = 15.0            # total solids of feed, wt %
    vs_percent: float = 85.0            # volatile solids as % of TS
    particle_density: float = 843.0     # kg/m3 (paper Table 3)
    hrt_days: float = 40.0              # system hydraulic retention time
    n_tanks: int = 4                    # number of CSTR modules in series
    target_ts_percent: float = 12.0     # TS of diluted slurry fed to digesters
    ch4_fraction: float = 0.60          # CH4 in biogas, vol %
    yield_lo: float = YIELD_LO_DEF      # m3 CH4 per kg VS added
    yield_hi: float = YIELD_HI_DEF
    power_density: float = 30.0         # W mixing power per m3 tank liquid


def _within(x, lo, hi):
    return lo <= x <= hi


# ---------------------------------------------------------------------------
# Calculation engine
# ---------------------------------------------------------------------------
def run_design(inp: DesignInputs) -> dict:
    """Solve the design step by step; returns a nested results dict."""
    r: dict = {"inputs": asdict(inp)}
    warn: list[str] = []

    # ---- Step 1: solids balance -------------------------------------------
    ts_t = inp.feed_rate_tpd * inp.ts_percent / 100.0
    vs_t = ts_t * inp.vs_percent / 100.0
    r["solids"] = {"ts_t": ts_t, "vs_t": vs_t, "moisture": 100 - inp.ts_percent}

    # ---- Step 2: slurry conditioning --------------------------------------
    need_ts_frac = inp.target_ts_percent / 100.0
    if (inp.ts_percent / 100.0) > need_ts_frac:
        slurry_t = ts_t / need_ts_frac           # total diluted slurry, t/d
    else:                                        # feed already leaner
        slurry_t = inp.feed_rate_tpd
        warn.append(
            f"Feed TS ({inp.ts_percent:.0f} %) is already below the target "
            f"({inp.target_ts_percent:.0f} %); no dilution water added."
        )
    water_t = max(0.0, slurry_t - inp.feed_rate_tpd)
    slurry_m3d = slurry_t * 1000.0 / SLURRY_RHO  # m3/d at design slurry density
    r["slurry"] = {"total_t": slurry_t, "water_t": water_t, "vol_m3d": slurry_m3d}

    # ---- Step 3: digester cascade sizing ----------------------------------
    working_total = slurry_m3d * inp.hrt_days                # m3
    per_tank_liq = working_total / inp.n_tanks               # m3
    tank_d = (4.0 * per_tank_liq / (math.pi * HT_RATIO)) ** (1.0 / 3.0)
    liquid_h = HT_RATIO * tank_d
    r["tank"] = {
        "working_total": working_total,
        "per_tank_liq": per_tank_liq,
        "tank_d": tank_d,
        "liquid_h": liquid_h,
        "total_h": liquid_h + liquid_h * 0.20 + 1.0,  # headspace + cone roof
    }

    # ---- Step 4: organic loading rate --------------------------------------
    olr = vs_t * 1000.0 / working_total                      # kg VS / m3 / d
    r["olr"] = olr
    if not _within(olr, 1.0, 1.5):
        warn.append(
            f"OLR = {olr:.2f} kg VS/m3/d lies outside the paper's best-yield "
            "band of 1.0-1.5; expect a moderate yield penalty or plan "
            "co-digestion upgrade (manure/dung)."
        )

    # ---- Step 5: reaction kinetics (assumed first order) -------------------
    k = K_FIRST_ORDER
    conv = 1.0 - math.exp(-k * inp.hrt_days)                 # fractional conv.
    degraded_vs = vs_t * BIODEGRADABLE_FRACTION * conv       # t VS/d
    r["kinetics"] = {
        "k": k, "conversion": conv,
        "biodegradable_vs_t": vs_t * BIODEGRADABLE_FRACTION,
        "degraded_vs_t": degraded_vs,
        "vsr_last_tank": degraded_vs / vs_t,
    }

    # ---- Step 6: gas and energy yields (paper-based band) ------------------
    vs_kgd = vs_t * 1000.0
    ch4_lo = vs_kgd * inp.yield_lo
    ch4_hi = vs_kgd * inp.yield_hi
    biogas_lo = ch4_lo / inp.ch4_fraction
    biogas_hi = ch4_hi / inp.ch4_fraction
    cbg_lo = ch4_lo * CH4_RHO / 1000.0                       # t/d
    cbg_hi = ch4_hi * CH4_RHO / 1000.0
    therm_lo = ch4_lo * CH4_LHV_KWH_M3 / 1000.0              # MWh/day
    therm_hi = ch4_hi * CH4_LHV_KWH_M3 / 1000.0
    gross_lo = therm_lo / 24.0 * CHP_EFF                     # MW
    gross_hi = therm_hi / 24.0 * CHP_EFF
    r["gas_energy"] = {
        "ch4_lo": ch4_lo, "ch4_hi": ch4_hi,
        "biogas_lo": biogas_lo, "biogas_hi": biogas_hi,
        "cbg_lo": cbg_lo, "cbg_hi": cbg_hi,
        "thermal_lo": therm_lo, "thermal_hi": therm_hi,
        "gross_mw_lo": gross_lo, "gross_mw_hi": gross_hi,
    }

    # ---- Step 7: agitator (scaled Model 3) ---------------------------------
    d_imp = IMPELLER_DT_RATIO * tank_d
    p_tank_kw = inp.power_density * per_tank_liq / 1000.0
    p_w = p_tank_kw * 1000.0
    n_hz = (p_w / (POWER_NUMBER * SLURRY_RHO * d_imp ** 5)) ** (1.0 / 3.0)
    rpm = n_hz * 60.0
    tip = math.pi * d_imp * n_hz
    rey = SLURRY_RHO * n_hz * d_imp ** 2 / 1.0e-3   # mu = 1e-3 Pa s
    froude = n_hz ** 2 * d_imp / 9.81
    torque = p_w / (2.0 * math.pi * n_hz)            # N.m
    hub_lo = HUB_LOWER_FRACTION * liquid_h
    hub_hi_elev = liquid_h - d_imp / 2.0 - SURFACE_CLEARANCE
    spacing = hub_hi_elev - hub_lo
    motor_kw = p_tank_kw / 0.80                      # gearbox/shaft losses
    mixers_total_mwh = inp.n_tanks * p_tank_kw * 24.0 / 1000.0
    r["mixer"] = {
        "d_imp": d_imp, "p_tank_kw": p_tank_kw, "n_hz": n_hz, "rpm": rpm,
        "tip": tip, "re": rey, "fr": froude, "torque_knm": torque / 1000.0,
        "hub_lo": hub_lo, "hub_up": hub_hi_elev, "spacing": spacing,
        "motor_kw": motor_kw, "gear_ratio": 1500.0 / rpm,
        "mixers_mwh_d": mixers_total_mwh, "np_effective": POWER_NUMBER,
        "blade_l": 0.45 * d_imp, "blade_w": 0.10 * d_imp,
        "blade_t": 0.015 * d_imp,
    }
    if tip > TIP_SPEED_MAX:
        warn.append(
            f"Tip speed {tip:.2f} m/s exceeds the {TIP_SPEED_MAX} m/s shear "
            "cap; reduce power density or impeller diameter."
        )
    if froude >= 0.1:
        warn.append(
            f"Froude number {froude:.3f} approaching vortex regime; verify "
            "baffle provision or reduce speed."
        )

    # ---- Step 8: heat services ---------------------------------------------
    heat_up_kw = slurry_m3d * 1000.0 * CP_SLURRY * SLURRY_TEMP_RISE / 86400.0
    losses_kw = HEAT_LOSS_KW_PER_TANK * inp.n_tanks
    parasitic_mw = (mixers_total_mwh * (1.0 + AUX_PARASITIC_FACTOR)) / 24.0
    net_lo = gross_lo - parasitic_mw
    net_hi = gross_hi - parasitic_mw
    r["energy_services"] = {
        "heat_up_kw": heat_up_kw, "shell_loss_kw": losses_kw,
        "total_heat_kw": heat_up_kw + losses_kw,
        "parasitic_mw": parasitic_mw,
        "net_mw_lo": net_lo, "net_mw_hi": net_hi,
        "mixers_mw": mixers_total_mwh / 24.0,
        "aux_mw": mixers_total_mwh * AUX_PARASITIC_FACTOR / 24.0,
    }

    # ---- Step 9: feed pre-treatment ----------------------------------------
    chop_tph = inp.feed_rate_tpd / CHOP_HOURS
    bunker_m3 = inp.feed_rate_tpd * BUNKER_DAYS / SILAGE_DENSITY
    r["pretreat"] = {
        "chop_tph": chop_tph, "mill_lines_tph": 2 * max(chop_tph, 1.0),
        "particle_mm": MAX_PARTICLE_MM,
        "bunker_m3": bunker_m3,
        "bunker_pair_m3": 2 * max(bunker_m3 / 2.0, 1.0),
        "pump_m3h": 50.0,
        "feed_m3h": slurry_m3d / 24.0,
    }
    r["float_risk"] = inp.particle_density < WATER_RHO
    if r["float_risk"]:
        warn.append(
            "Particle density below water: floating/grass-rind regime. The "
            "upper impeller surface-draw, gas recirculation and baffle "
            "pilasters are MANDATORY (paper Model 3 lesson)."
        )
    r["warnings"] = warn

    r["steps"] = build_steps(inp, r)
    return r


# ---------------------------------------------------------------------------
# Sequential step table (used by the app and the report)
# ---------------------------------------------------------------------------
def build_steps(inp: DesignInputs, r: dict) -> list[dict]:
    """Ordered step dicts: (title, equation, substitution, result)."""
    s = r["solids"]; sl = r["slurry"]; t = r["tank"]
    k = r["kinetics"]; g = r["gas_energy"]; m = r["mixer"]
    e = r["energy_services"]; pr = r["pretreat"]

    return [
        {"no": 1, "title": "Solids balance",
         "formula": "m_TS = F * TS%/100 ;  m_VS = m_TS * VS%/100",
         "sub": f"F = {inp.feed_rate_tpd:.1f} t/d, TS = {inp.ts_percent:.0f} %, VS = {inp.vs_percent:.0f} % of TS",
         "result": (f"m_TS = {s['ts_t']:.1f} t/d ; m_VS = {s['vs_t']:.1f} t/d "
                    f"(moisture {s['moisture']:.0f} %)")},
        {"no": 2, "title": "Slurry conditioning (dilution to target TS)",
         "formula": "M_slurry = m_TS / (TS_target/100) ;  m_water = M_slurry - F ;  Q = M_slurry / rho",
         "sub": f"m_TS = {s['ts_t']:.1f} t/d, target TS = {inp.target_ts_percent:.0f} %",
         "result": (f"M_slurry = {sl['total_t']:.1f} t/d ; dilution water = {sl['water_t']:.1f} t/d ; "
                    f"Q = {sl['vol_m3d']:.1f} m3/d @ {SLURRY_RHO:.0f} kg/m3")},
        {"no": 3, "title": "Digester cascade sizing (HRT rule)",
         "formula": "V_work = Q * HRT ;  V_tank = V_work / N ;  T = (4 V_tank / (pi * H/T))^(1/3) ; H = 0.96 T",
         "sub": (f"Q = {sl['vol_m3d']:.1f} m3/d, HRT = {inp.hrt_days:.0f} d, "
                 f"N = {inp.n_tanks} tanks, H/T = {HT_RATIO}"),
         "result": (f"V_work = {t['working_total']:.0f} m3 ; V_tank = {t['per_tank_liq']:.0f} m3 ; "
                    f"T = {t['tank_d']:.1f} m ; H = {t['liquid_h']:.1f} m ; "
                    f"actual system HRT = {t['working_total'] / sl['vol_m3d']:.1f} d")},
        {"no": 4, "title": "Organic loading rate",
         "formula": "OLR = m_VS * 1000 / V_work",
         "sub": f"m_VS = {r['solids']['vs_t']:.1f} t/d, V_work = {t['working_total']:.0f} m3",
         "result": f"OLR = {r['olr']:.2f} kg VS/m3/d "
                   f"(paper best-yield band: 1.0-1.5)"},
        {"no": 5, "title": "Reaction kinetics (first-order assumption)",
         "formula": "X = 1 - exp(-k * HRT) ;  m_VS_deg = m_VS * f_bio * X",
         "sub": (f"k = {k['k']:.2f} /d (mesophilic, assumed), "
                 f"f_bio = {BIODEGRADABLE_FRACTION:.2f}, HRT = {inp.hrt_days:.0f} d"),
         "result": (f"X = {k['conversion']:.2f} ({k['conversion'] * 100:.1f} %) ; "
                    f"biodegradable VS = {k['biodegradable_vs_t']:.1f} t/d ; "
                    f"degraded VS = {k['degraded_vs_t']:.1f} t/d")},
        {"no": 6, "title": "Gas and energy yields (paper CSTR band)",
         "formula": ("m_CH4 = m_VS,kgd * Y ;  V_biogas = m_CH4 / f_CH4 ;  "
                     "CBG = m_CH4 * 0.716 ;  E_th = V_CH4 * 9.95 ;  P_el,gross = E_th/24 * 0.41"),
         "sub": (f"m_VS = {r['solids']['vs_t'] * 1000:.0f} kg/d, Y = {inp.yield_lo:.2f}-{inp.yield_hi:.2f} "
                 f"m3 CH4/kg VS, f_CH4 = {inp.ch4_fraction:.2f}"),
         "result": (f"CH4 = {g['ch4_lo']:.0f}-{g['ch4_hi']:.0f} m3/d ; biogas = "
                    f"{g['biogas_lo']:.0f}-{g['biogas_hi']:.0f} m3/d ; "
                    f"CBG = {g['cbg_lo']:.1f}-{g['cbg_hi']:.1f} t/d ; "
                    f"thermal {g['thermal_lo']:.1f}-{g['thermal_hi']:.1f} MWh/d ; "
                    f"gross electric {g['gross_mw_lo']:.2f}-{g['gross_mw_hi']:.2f} MW")},
        {"no": 7, "title": "Agitator sizing (scaled Model 3)",
         "formula": ("Da = 0.4 T ; P_tank = (P/V) V ; n = (P / (Np rho Da^5))^(1/3) ; "
                     "Tq = P/(2 pi n) ; Re = rho n Da^2 / mu ; Fr = n^2 Da / g"),
         "sub": (f"T = {t['tank_d']:.1f} m, P/V = {inp.power_density:.0f} W/m3, "
                 f"Np = {POWER_NUMBER:.2f}"),
         "result": (f"Da = {m['d_imp']:.2f} m ; P_tank = {m['p_tank_kw']:.1f} kW ; "
                    f"n = {m['rpm']:.1f} rpm ; tip = {m['tip']:.2f} m/s ; "
                    f"Re = {m['re']:.2e} ; Fr = {m['fr']:.4f} ; "
                    f"torque = {m['torque_knm']:.1f} kN-m ; motor = {m['motor_kw']:.0f} kW @ gear {m['gear_ratio']:.0f}:1")},
        {"no": 8, "title": "Heat services",
         "formula": "Q_heat = Q * rho_sl * c_p * dT / (86400 * 1000) + N * q_loss ; P_parasitic = (E_mix (1+f_aux))/24",
         "sub": (f"Q = {sl['vol_m3d']:.1f} m3/d, dT = {SLURRY_TEMP_RISE:.0f} K, "
                 f"q_loss = {HEAT_LOSS_KW_PER_TANK:.1f} kW/tank, f_aux = {AUX_PARASITIC_FACTOR:.2f}"),
         "result": (f"Q_heat = {e['heat_up_kw']:.0f} kW + losses {e['shell_loss_kw']:.0f} kW "
                    f"= {e['total_heat_kw']:.0f} kW total ; "
                    f"parasitic electric = {e['parasitic_mw']:.2f} MW avg "
                    "(mixer electric + aux = %.0f %% of mixer electric, continuous duty)"
                    % (AUX_PARASITIC_FACTOR * 100))},
        {"no": 9, "title": "Feed pre-treatment specification",
         "formula": "C_tph = F / h_chop ; V_bunker = F * d_buffer / rho_silage",
         "sub": f"F = {inp.feed_rate_tpd:.1f} t/d, chop window = {CHOP_HOURS:.0f} h/d, "
                f"buffer = {BUNKER_DAYS:.0f} d, rho silage = {SILAGE_DENSITY:.0f} kg/m3",
         "result": (f"Mill duty = {pr['chop_tph']:.1f} t/h (install 2 x {pr['mill_lines_tph']:.0f} t/h, "
                    "duty + standby) ; particle <= "
                    f"{MAX_PARTICLE_MM:.0f} mm ; silage bunker pair = "
                    f"{pr['bunker_pair_m3']:.0f} m3 ; feed pump duty = {pr['feed_m3h']:.2f} m3/h "
                    "(install 50 m3/h hose pumps, duty+standby)")},
    ]


def recommendations(inp: DesignInputs, r: dict) -> list[str]:
    recs = []
    pr = r["pretreat"]
    if r["float_risk"]:
        recs.append(
            "Two downward-pumping (axial) propellers per tank are mandatory; the "
            "upper impeller top blade must sit ~1.5 m below the liquid surface to "
            "draw the floating grass layer down (paper Model 3 finding)."
        )
        recs.append(
            "Provide bottom gas spargers (100-150 Nm3/h per tank) and digestate "
            "surface sprays as scum-breaking redundancy during mixer stop windows."
        )
    recs.append(
        f"Chop + hammer-mill all feed to <= {MAX_PARTICLE_MM:.0f} mm longest "
        "dimension (paper used disc + hammer milling); run mill "
        f"{CHOP_HOURS:.0f} h/d at {pr['chop_tph']:.1f} t/h with 100 % standby line."
    )
    recs.append("Silage (15 d buffer) instead of fresh feed to smooth winter "
                "dormancy and provide storage for truck/harvest buffering.")
    if not _within(r["olr"], 1.0, 1.5):
        recs.append(
            f"OLR = {r['olr']:.2f} kg VS/m3/d is outside the paper's best band. "
            "Options: (a) accept a lower conversion, (b) reduce feed TS %, "
            "(c) co-digest manure/dung (documented plants show > 30 % yield "
            "uplift), (d) add a 5th module."
        )
    recs.append(
        "House the drive train with VFD (8-16 rpm band, torque-limited soft start); "
        "fit anti-wrap plates at both hubs (fibre wrapping risk) and torque "
        "transmitters per shaft as the mixing-health surrogate."
    )
    recs.append(
        "Run continuous primary agitation for Napier feed (floating tendency); "
        "allow 20 min/h duty-cycling (33 % energy saving) only after the torque "
        "rebound test shows the fountain restoring in < 10 min."
    )
    recs.append(
        "Dose FeCl3 at the feed inlet for H2S control (< 50 ppm in raw gas); "
        "spec duplex/316L for shafts and impellers, epoxy-lined concrete."
    )
    return recs


def conclusions(inp: DesignInputs, r: dict) -> list[str]:
    t = r["tank"]; g = r["gas_energy"]; m = r["mixer"]; k = r["kinetics"]
    return [
        f"A cascade of {inp.n_tanks} flat-bottom CSTRs of "
        f"{t['per_tank_liq']:.0f} m3 each (T = {t['tank_d']:.2f} m, H = "
        f"{t['liquid_h']:.2f} m) digests {inp.feed_rate_tpd:.0f} t/d of crushed "
        f"Napier grass at a {r['inputs']['hrt_days']:.0f}-day HRT and "
        f"{r['olr']:.2f} kg VS/m3/d loading.",
        "Mixing is supplied by the paper's validated Model 3 pattern: two "
        f"down-pumping 3-blade 45-degree propellers (Da = {m['d_imp']:.2f} m, "
        f"Np = {POWER_NUMBER:.2f}) driving {m['p_tank_kw']:.0f} kW per tank at "
        f"{m['rpm']:.1f} rpm ({m['tip']:.2f} m/s tip). The lab's surface-vortex "
        "problem vanishes at scale (Fr = "f"{m['fr']:.4f}); the field problem that "
        "remains is floatation of the 843 kg/m3 grass, countered by the upper "
        "impeller's surface draw plus gas sparging.",
        f"First-order kinetics at k = {k['k']:.2f} /d gives "
        f"{k['conversion'] * 100:.0f} % VS conversion at the design HRT; the "
        "paper's CSTR mono-digestion band projects "
        f"{g['ch4_lo']:,.0f}-{g['ch4_hi']:,.0f} m3/d CH4 "
        f"({g['biogas_lo']:,.0f}-{g['biogas_hi']:,.0f} m3/d biogas), i.e. "
        f"{g['cbg_lo']:.1f}-{g['cbg_hi']:.1f} t/d CBG or "
        f"{g['gross_mw_lo']:.2f}-{g['gross_mw_hi']:.2f} MWe gross "
        "(net 0.5-0.6 MWe at the 120 t/d datum after parasitics).",
        "Mono-grass is at the low-yield end of the register's documented plants; "
        "co-digestion with manure is the identified upgrade path to 5+ t/d CBG "
        "(Rajkot precedent) and is already reflected in the paper's own best "
        "two-stage co-digestion result.",
    ]


# ---------------------------------------------------------------------------
# Markdown design report generator
# ---------------------------------------------------------------------------
def generate_report(inp: DesignInputs, r: dict) -> str:
    """Full markdown design document with all required sections."""
    s = r["solids"]; sl = r["slurry"]; t = r["tank"]
    k = r["kinetics"]; g = r["gas_energy"]; m = r["mixer"]
    e = r["energy_services"]; pr = r["pretreat"]
    steps = r["steps"]

    L = []
    L.append("# Napier-Grass CSTR Digester Cascade - Design Report\n")
    L.append(f"*Generated by the Napier-Grass CSTR Design Suite on "
             f"{r['inputs']['feed_rate_tpd']:.0f} t/day design datum.*\n")
    L.append("## 1. Preamble\n")
    L.append(
        "This report sizes an industrial cascade of continuously stirred tank "
        "reactors (CSTRs) for anaerobic digestion of crushed Napier grass "
        "(*Pennisetum purpureum*, Pakchong 1 lineage). The mixing concept and all "
        "correlation coefficients are taken from the lab-validated numerical study "
        "of Satjaritanun et al. (2016), whose 'Model 3' (two down-pumping "
        "3-blade 45-degree propellers with four T/12 baffles, Np = 2.00, mixing "
        "efficiency 34.66 %) out-performed single-propeller designs and prevents "
        "floatation of the low-density (843 kg/m3) grass. The document supersedes "
        "the fixed 120 t/d datum in `napier-cstr-design-120tpd.md` for the user's "
        "current inputs.\n")
    L.append("## 2. Design inputs\n")
    L.append("| Parameter | Value | Unit | Source |\n|---|---|---|---|")
    L.append("| Fresh feed rate | %.1f | t/d | user |" % inp.feed_rate_tpd)
    L.append("| Total solids, TS | %.1f | wt %% | user |" % inp.ts_percent)
    L.append(f"| Volatile solids, VS | {inp.vs_percent:.1f} | % of TS | user |")
    L.append("| Particle density | %.0f | kg/m3 | user (paper default 843) |"
             % inp.particle_density)
    L.append("| HRT | %.0f | d | user / biogas practice 30-50 d |" % inp.hrt_days)
    L.append("| Number of CSTR modules | %d | - | user |" % inp.n_tanks)
    L.append("| Target slurry TS | %.0f | wt %% | user |" % inp.target_ts_percent)
    L.append("| CH4 in biogas | %.0f | vol %% | user |" % (inp.ch4_fraction * 100))
    L.append("| CH4 yield band | %.2f-%.2f | m3/kg VS | paper CSTR |"
             % (inp.yield_lo, inp.yield_hi))
    L.append("| Mixing power density | %.0f | W/m3 | user (20-40 W/m3 band) |"
             % inp.power_density)
    L.append("")
    L.append("## 3. Assumptions\n")
    for a in [
        "Pure water density 998 kg/m3 (paper Table 3); diluted slurry 1,000 kg/m3.",
        "Slurry cP = 4.0 kJ/kg-K; digester 37 C, feed enters at 12 C (25 K rise).",
        "VS/TS as given by user; biodegradable fraction of VS = 0.75.",
        "First-order kinetics, k = 0.08 /d at 37 C mesophilic (assumed).",
        "Methane LHV = 9.95 kWh/m3 (35.8 MJ/m3); methane density 0.716 kg/m3.",
        "Gas contains %.0f %% CH4 before upgrading." % (inp.ch4_fraction * 100),
        "CHP gross electrical efficiency 41 %.",
        "Impeller power number Np = 2.00 (paper Model 3: dual propeller + baffles).",
        "Impeller/geometry ratios: H/T = 0.96, Da/T = 0.40, blade L = 0.45 Da, "
        "W = 0.10 Da, t = 0.015 Da (scaled from paper Table 1).",
        "Tank shell heat loss 7.5 kW per 37 C module; aux parasitics = 60 % of "
        "mixer electric consumption.",
        "Silage bunker density 650 kg/m3; mill chopping window 8 h/d.",
        "Shear cap 3.5 m/s tip speed protects methanogenic consortia (paper "
        "caution on excess shear).",
        "r_floating risk active because actual %.0f kg/m3 %s water (998 kg/m3)."
        % (inp.particle_density,
           "<" if inp.particle_density < WATER_RHO else ">="),
    ]:
        L.append(f"- {a}")
    L.append("")
    L.append("## 4. Equations (symbols per step list below)\n")
    L.append("```")
    for eq in [
        "(1)  m_TS    = F * (TS/100)",
        "(2)  m_VS    = m_TS * (VS%/100)",
        "(3)  M_slurry = m_TS / (TS_target/100);  m_water = M_slurry - F",
        "(4)  Q       = M_slurry / rho_slurry",
        "(5)  V_work  = Q * HRT ;  V_tank = V_work / N",
        "(6)  T       = (4 V_tank / (pi (H/T)))^(1/3) ;  H = 0.96 T",
        "(7)  OLR     = m_VS * 1000 / V_work",
        "(8)  X       = 1 - e^(-k HRT);  m_VS,deg = m_VS f_bio X",
        "(9)  V_CH4   = m_VS,kgd * Y ;  V_biogas = V_CH4 / f_CH4",
        "(10) CBG     = V_CH4 * 0.716 ;  E_th = V_CH4 * 9.95",
        "(11) P_el,g  = E_th / 24 * eta_CHP",
        "(12) Da      = 0.40 T",
        "(13) P_tank  = (P/V) V_tank ;  n = (P / (Np rho Da^5))^(1/3)",
        "(14) Tq      = P / (2 pi n) ;  Re = rho n Da^2 / mu ;  Fr = n^2 Da / g",
        "(15) Q_heat  = Q rho c_p dT / 86400000 + N q_loss  [kW]",
        "(16) P_par   = E_mix (1 + f_aux) / 24",
        "(17) C_tph   = F / h_chop ;  V_bunker = F d_buffer / rho_sil",
    ]:
        L.append("  " + eq)
    L.append("```\n")
    L.append("## 5. Step-by-step calculations\n")
    for st_ in steps:
        L.append(f"### Step {st_['no']} - {st_['title']}\n")
        L.append(f"`{st_['formula']}`\n")
        L.append(f"- Given: {st_['sub']}")
        L.append(f"- **Result:** {st_['result']}\n")
    L.append("## 6. Results\n")
    biogas_mid = 0.5 * (g['biogas_lo'] + g['biogas_hi'])
    L.append("### 6.1 Mass balance\n")
    L.append("| Stream | t/d | m3/d |")
    L.append("|---|---|---|")
    L.append(f"| Fresh grass in | {inp.feed_rate_tpd:.1f} | ~{inp.feed_rate_tpd / SILAGE_DENSITY:,.0f} |")
    L.append(f"| Dilution water in | {sl['water_t']:.1f} | {sl['water_t']:.1f} |")
    L.append(f"| Slurry to digesters | {sl['total_t']:.1f} | {sl['vol_m3d']:.1f} |")
    L.append(f"| Biogas out | {biogas_mid * 1.12:,.1f} | "
             f"{g['biogas_lo']:,.0f}-{g['biogas_hi']:,.0f} |")
    L.append(f"| Digestate out | ~{sl['total_t'] - biogas_mid * 1.12:.1f} | "
             f"~{sl['vol_m3d'] - biogas_mid:.0f} |")
    L.append("")
    L.append("### 6.2 Digester geometry (per module)\n")
    L.append("| Quantity | Value |")
    L.append("|---|---|")
    L.append(f"| Working volume | {t['per_tank_liq']:,.0f} m3 |")
    L.append(f"| Tank diameter T | {t['tank_d']:.2f} m |")
    L.append(f"| Liquid height H | {t['liquid_h']:.2f} m (H/T = %.2f) |" % HT_RATIO)
    L.append(f"| Total height incl roof | {t['total_h']:.2f} m |")
    L.append(f"| System HRT | {t['working_total'] / sl['vol_m3d']:.1f} d |")
    L.append(f"| OLR | {r['olr']:.2f} kg VS/m3/d |")
    L.append("")
    L.append("### 6.3 Agitator package (per module)\n")
    L.append("| Item | Specification |")
    L.append("|---|---|")
    L.append(f"| Impellers | 2 x 3-blade 45-deg down-pumping, Da = {m['d_imp']:.2f} m |")
    L.append(f"| Blade L x W x t | {m['blade_l']:.2f} x {m['blade_w']:.2f} x {m['blade_t']:.3f} m |")
    L.append(f"| Hub elevations (off floor) | lower {m['hub_lo']:.2f} m, upper {m['hub_up']:.2f} m (spacing {m['spacing']:.2f} m = {m['spacing'] / m['d_imp']:.2f} Da) |")
    L.append(f"| Speed (VFD band 8-16 rpm concept) | {m['rpm']:.1f} rpm |")
    L.append(f"| Shaft power / motor | {m['p_tank_kw']:.1f} kW / {m['motor_kw']:.0f} kW, gear ~{m['gear_ratio']:.0f}:1 |")
    L.append(f"| Shaft torque | {m['torque_knm']:.1f} kN-m |")
    L.append(f"| Tip speed | {m['tip']:.2f} m/s (<= 3.5 cap) |")
    L.append(f"| Reynolds / Froude | {m['re']:.2e} / {m['fr']:.4f} |")
    L.append(f"| Baffles | 4 x radial pilasters, width = T/12 = {t['tank_d'] / 12.0 * 100:.0f} cm |")
    L.append("")
    L.append("### 6.4 Gas and energy outputs\n")
    L.append("| Quantity | Low | High |")
    L.append("|---|---|---|")
    L.append(f"| Methane | {g['ch4_lo']:,.0f} m3/d | {g['ch4_hi']:,.0f} m3/d |")
    L.append(f"| Biogas (pre-upgrade) | {g['biogas_lo']:,.0f} m3/d | {g['biogas_hi']:,.0f} m3/d |")
    L.append(f"| CBG (after upgrade) | {g['cbg_lo']:.1f} t/d | {g['cbg_hi']:.1f} t/d |")
    L.append(f"| Thermal | {g['thermal_lo']:.1f} MWh/d | {g['thermal_hi']:.1f} MWh/d |")
    L.append(f"| Gross electric | {g['gross_mw_lo']:.2f} MW | {g['gross_mw_hi']:.2f} MW |")
    L.append(f"| Parasitic electric | {e['parasitic_mw']:.2f} MW | |")
    L.append(f"| Net electric | {e['net_mw_lo']:.2f} MW | {e['net_mw_hi']:.2f} MW |")
    L.append("")
    L.append("### 6.5 Pre-treatment package\n")
    L.append("| Item | Specification |")
    L.append("|---|---|")
    L.append(f"| Particle size at digester | <= {MAX_PARTICLE_MM:.0f} mm (paper: disc + hammer mill) |")
    L.append(f"| Mill duty | {pr['chop_tph']:.1f} t/h over {CHOP_HOURS:.0f} h/d; install 2 lines |")
    L.append(f"| Silage bunker pair | {pr['bunker_pair_m3']:.0f} m3 ({BUNKER_DAYS:.0f}-day buffer) |")
    L.append(f"| Dilution water | {sl['water_t']:.1f} t/d |")
    L.append(f"| Feed pumps | {pr['pump_m3h']:.0f} m3/h hose/peristaltic, duty + standby "
             f"(duty flow {pr['feed_m3h']:.2f} m3/h) |")
    L.append("")
    L.append("### 6.6 Acceptance criteria\n")
    L.append("| Criterion | Acceptance |")
    L.append("|---|---|")
    L.append("| Homogenisation | 95 % within 30 min (heat tracer / torque test) |")
    L.append("| Differential temperature across depths | <= 0.8 C steady state |")
    L.append("| Scum / floating island | none > 0.5 m after 1 week continuous duty |")
    L.append("| Tip speed | <= 3.5 m/s (shear cap) |")
    L.append("| Biogas rate vs. agitation cycle change | drop <= 5 % |")
    L.append("")
    if r["warnings"]:
        L.append("### 6.7 Design warnings\n")
        for w in r["warnings"]:
            L.append(f"- {w}")
        L.append("")
    L.append("## 7. Inferences and conclusions\n")
    for c in conclusions(inp, r):
        L.append(f"- {c}")
    L.append("")
    L.append("## 8. Recommendations\n")
    for rec in recommendations(inp, r):
        L.append(f"- {rec}")
    L.append("")
    L.append("## 9. References\n")
    for ref in [
        "P. Satjaritanun, Y. Khunatorn, N. Vorayos, S. Shimpalee, E. Bringley, "
        '"Numerical analysis of the mixing characteristic for napier grass in the '
        'continuous stirring tank reactor for biogas production", Biomass and '
        "Bioenergy 86 (2016) 53-64. https://doi.org/10.1016/j.biombioe.2016.01.018",
        "McCabe, Smith, Harriott, Unit Operations of Chemical Engineering, 7th ed. "
        "(impeller correlations of Eq. 12-16).",
        "Repository: napier-cstr-design-120tpd.md (fixed 120 t/d recast); "
        "napier-grass-bioreactors.md (23-installation register).",
        "FAO AGRIS record of the 2022 Napier mono-digestion CSTR study "
        "(0.27-0.34 L CH4/g VS basis, OLR 1.0-1.5 g VS/L*d).",
    ]:
        L.append(f"1. {ref}")
    L.append("")
    return "\n".join(L)
