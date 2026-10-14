"""
Napier-Grass CSTR Design Suite - Streamlit front end.

Run with:
    streamlit run app.py
"""
import math

import pandas as pd
import streamlit as st

from design_core import (
    DesignInputs, run_design, generate_report,
    recommendations, conclusions, WATER_RHO,
)

st.set_page_config(page_title="Napier Grass CSTR Design Suite",
                   page_icon="🌿", layout="wide", initial_sidebar_state="expanded")

# ---------------------------------------------------------------------------
# Sidebar: design inputs
# ---------------------------------------------------------------------------
st.sidebar.title("🌿 Design inputs")
st.sidebar.caption("Napier-grass CSTR digester cascade (Satjaritanun 2016 basis)")

feed_rate = st.sidebar.number_input(
    "Feed rate (t/day, as received)", min_value=10.0, max_value=1_000.0,
    value=120.0, step=5.0, help="Fresh crushed Napier grass fed per day.")
ts_pct = st.sidebar.slider("Total Solids, TS (wt %)", 5.0, 30.0, 15.0, 0.5)
vs_pct = st.sidebar.slider("Volatile Solids, VS (% of TS)", 50.0, 100.0, 85.0, 1.0)
rho_p = st.sidebar.number_input(
    "Particle density (kg/m³)", min_value=600.0, max_value=1_200.0,
    value=843.0, step=5.0, help="Paper default 843 kg/m³ — below water → "
    "floating-rind regime.")

with st.sidebar.expander("Advanced parameters"):
    hrt_days = st.slider("HRT (days)", 20.0, 60.0, 40.0, 1.0)
    n_tanks = st.selectbox("Number of CSTR modules (series)", [2, 3, 4, 5, 6], 3)
    target_ts = st.slider("Target slurry TS (wt %)", 8.0, 20.0, 12.0, 0.5)
    ch4_frac = st.slider("CH₄ in biogas (vol %)", 50.0, 70.0, 60.0, 1.0) / 100.0
    yield_lo = st.number_input("CH₄ yield low (m³/kg VS)", 0.15, 0.40, 0.27, 0.01)
    yield_hi = st.number_input("CH₄ yield high (m³/kg VS)", 0.15, 0.50, 0.34, 0.01)
    power_density = st.slider("Mixing power density (W/m³)", 15.0, 60.0, 30.0, 1.0)

inputs = DesignInputs(
    feed_rate_tpd=float(feed_rate), ts_percent=float(ts_pct),
    vs_percent=float(vs_pct), particle_density=float(rho_p),
    hrt_days=float(hrt_days), n_tanks=int(n_tanks),
    target_ts_percent=float(target_ts), ch4_fraction=float(ch4_frac),
    yield_lo=float(yield_lo), yield_hi=float(yield_hi),
    power_density=float(power_density),
)
r = run_design(inputs)
report_md = generate_report(inputs, r)

# ---------------------------------------------------------------------------
# Main workspace
# ---------------------------------------------------------------------------
st.title("🌿 Napier-Grass CSTR Digester Cascade — Design Suite")
st.caption(
    "Scales the lab-validated “Model 3” mixing concept of Satjaritanun et al. "
    "*Biomass and Bioenergy* 86 (2016) 53–64 to plant scale.")

if r["warnings"]:
    for w in r["warnings"]:
        st.warning(w)

pr = r["pretreat"]
sg = r["slurry"]; tk = r["tank"]
c1, c2, c3, c4 = st.columns(4)
c1.metric("Slurry to digesters", f"{sg['vol_m3d']:,.0f} m³/d",
          f"{sg['total_t']:,.0f} t/d")
c2.metric("Total working volume", f"{tk['working_total']:,.0f} m³",
          f"HRT {tk['working_total']/sg['vol_m3d']:.0f} d")
c3.metric("OLR", f"{r['olr']:.2f} kg VS/m³·d",
          "band 1.0–1.5" if 1.0 <= r["olr"] <= 1.5 else "⚠ outside band")
c4.metric("Floating risk", "ACTIVE" if r["float_risk"] else "low",
          f"ρ particle {inputs.particle_density:.0f} kg/m³")

# ---- Pre-treatment specifications block ------------------------------------
st.subheader("⚙️ Feed pre-treatment specifications required")
st.markdown(
    f"""
| Stage | Specification required |
|---|---|
| **Storage / silage** | {pr['bunker_pair_m3']:,.0f} m³ bunker pair ({15:.0f}-day buffer at {pr['chop_tph']:.1f} t/d → silaged at 650 kg/m³); covers winter dormancy & harvest gaps |
| **Chopping** | Disc cutter: reduce stalk to ≤ 50 mm; duty {pr['chop_tph']:.1f} t/h over {8:.0f} h/d |
| **Hammer milling** | Install 2 × {pr['mill_lines_tph']:.0f} t/h lines (duty + standby); final particle **≤ {pr['particle_mm']:.0f} mm longest dimension** (paper used disc + hammer milling; 2.5 mm model particles) |
| **Compaction/handling** | Airtight silage handling to limit aerobic losses (paper §2.2) |
| **Dilution** | {sg['water_t']:,.1f} t/d make-up water inline to **{inputs.target_ts_percent:.0f} % TS** slurry before pumps |
| **Maceration in line** | Inline macerator upstream of hose pumps to stop strand re-growth |
| **Pumping** | {pr['pump_m3h']:.0f} m³/h hose/peristaltic pumps, duty + standby; duty flow {pr['feed_m3h']:.2f} m³/h continuous |
| **Feeding point** | DN200 bottom-injection toward the lower-impeller sweep (mirrors paper's bottom feed pipe) |
| **Floatation countermeasures** | {'**Mandatory** — particle density ' + f'{inputs.particle_density:.0f}' + ' < water 998 kg/m³' if r['float_risk'] else 'Conventional setup sufficient'} |
""")
st.info("Downstream digester results update live with the left-panel inputs — "
        "see the tabs below.", icon="👇")

tabs = st.tabs([
    "🧪 Reaction kinetics & assumptions",
    "➡️ Outputs (in sequence)",
    "🧭 Recommendations",
    "✅ Conclusions",
    "📄 Design report & download",
])

# ---- Tab 1: kinetics --------------------------------------------------------
with tabs[0]:
    k = r["kinetics"]
    st.subheader("Reaction kinetics model")
    st.markdown(
        "**First-order (substrate-limited) model** — standard for high-fibre "
        "mesophilic digestion:\n\n"
        "$X = 1 - e^{-\\,k \\, HRT}$  ·  "
        "$\\dot m_{VS,deg} = \\, \\dot m_{VS} f_{bio} X$"
    )
    ka1, ka2, ka3, ka4 = st.columns(4)
    ka1.metric("Rate constant k (assumed)", f"{k['k']:.2f} /d", "37 °C mesophilic")
    ka2.metric("VS conversion X", f"{k['conversion']*100:.1f} %")
    ka3.metric("Biodegradable VS", f"{k['biodegradable_vs_t']:.1f} t/d", "f_bio = 0.75")
    ka4.metric("Degraded VS", f"{k['degraded_vs_t']:.1f} t/d", f"VSr = {k['vsr_last_tank']*100:.0f} % of feed VS")

    st.markdown("**Model assumptions**")
    st.markdown(
        f"""
- HRT-based first-order degradation; hydrolysis of lignocellulose is
  rate-limiting (paper: 41 % cellulose / 14 % hemicellulose / 13 % lignin).
- k = **{k['k']:.2f} /d** (assumed, mesophilic 37 °C, grass class 30–50 d HRT).
- Biodegradable fraction of VS **f_bio = 0.75** (grass assumption).
- Empirical CH₄ yield band **{inputs.yield_lo:.2f}–{inputs.yield_hi:.2f} m³/kg VS**
  (paper's CSTR mono-digestion at OLR 1.0–1.5 g VS/L·d) is used for output
  sizing — it already embeds real digester inhibition, so it is authoritative
  over the theoretical conversion result.
- Loading-rate check: paper yields fall at OLR ≥ 2.0 kg VS/m³·d
  (CSTR study); current design OLR is **{r['olr']:.2f}** —
  see warning area if flagged.
- Gypsum/trace nutrients, alkalinity control (pH 7.2–7.6) assumed available
  via feed conditioning; VFA slip monitored per-stage.
""")
    hrt_curve = pd.DataFrame({"HRT (d)": list(range(15, 61, 5))})
    hrt_curve["Conversion X (%)"] = hrt_curve["HRT (d)"].apply(
        lambda h: round(100.0 * (1.0 - math.exp(-k["k"] * h)), 1))
    hrt_curve = hrt_curve.set_index("HRT (d)")
    st.line_chart(hrt_curve, height=300)

# ---- Tab 2: outputs in sequence ---------------------------------------------
with tabs[1]:
    st.subheader("Design outputs — sequential calculation trace")
    for stp in r["steps"]:
        with st.container(border=True):
            st.markdown(f"**Step {stp['no']} — {stp['title']}**")
            st.code(stp["formula"], language="text")
            st.markdown(f"- Given: {stp['sub']}")
            st.markdown(f"- **Result:** {stp['result']}")

# ---- Tab 3: recommendations ---------------------------------------------------
with tabs[2]:
    st.subheader("Recommendations")
    for rec in recommendations(inputs, r):
        st.markdown(f"- {rec}")

# ---- Tab 4: conclusions -------------------------------------------------------
with tabs[3]:
    st.subheader("Conclusions")
    g = r["gas_energy"]; m = r["mixer"]; e = r["energy_services"]
    gc1, gc2, gc3, gc4 = st.columns(4)
    gc1.metric("Biogas (pre-upgrade)",
               f"{g['biogas_lo']:,.0f}–{g['biogas_hi']:,.0f} m³/d")
    gc2.metric("CBG potential", f"{g['cbg_lo']:.1f}–{g['cbg_hi']:.1f} t/d")
    gc3.metric("Gross electric", f"{g['gross_mw_lo']:.2f}–{g['gross_mw_hi']:.2f} MWe")
    gc4.metric("Process heat duty", f"{e['total_heat_kw']:,.0f} kW<sub>th</sub>",
               "self-sufficient", help=None)
    for c in conclusions(inputs, r):
        st.markdown(f"- {c}")
    st.download_button(
        "⬇️ Download design document (.md)",
        data=report_md, file_name="napier_cstr_design_report.md",
        mime="text/markdown", type="primary")

# ---- Tab 5: design report ------------------------------------------------------
with tabs[4]:
    st.subheader("Generated design report (markdown)")
    st.download_button(
        "⬇️ Download design document (.md)",
        data=report_md, file_name="napier_cstr_design_report.md",
        mime="text/markdown", type="primary",
        help="Complete report: preamble, inputs table, equations, assumptions, "
             "step-by-step calculations, result tables, inferences, "
             "recommendations, references.")
    st.markdown(report_md)
