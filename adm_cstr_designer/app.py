"""Anaerobic CSTR designer - Streamlit front end.

Run:  streamlit run app.py   (from the adm_cstr_designer folder)
"""
from __future__ import annotations

import dataclasses
import json
import sys
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st

sys.path.insert(0, str(Path(__file__).resolve().parent))

from core.data import load_feedstocks, load_kinetics, feedstock_meta
from core.design import DesignInputs, run_design
from core.feed import FeedInputs, blend_properties
from core.kin_curves import inhibition_curves, kinetics_table, monod_curves
from core.scenarios import sweep
from cad.geometry_model import cad_params
from cad.preview import build_parts, plotly_figure, render_png, to_stl_bytes
from cad.freecad_build import freecad_executable, make_script, run_freecad
from report.docx_report import build_report

st.set_page_config(page_title='Anaerobic CSTR Designer', layout='wide', page_icon=':material/science:')

FEED = load_feedstocks()
KIN = load_kinetics()

# ---------------------------------------------------------------- state
ss = st.session_state
ss.setdefault('design', None)
ss.setdefault('scenarios', None)
ss.setdefault('stamp', None)


def pct_inputs() -> dict:
    return {k: ss.get(f'blend_{k}', 0.0) for k in FEED}


# ---------------------------------------------------------------- sidebar
with st.sidebar:
    st.title('CSTR Designer')
    project = st.text_input('Project name', 'Anaerobic CSTR design')
    run_clicked = st.button('Run design', type='primary', use_container_width=True)
    st.divider()
    st.caption('Preliminary design tool. Feedstock data and some constants are indicative defaults; '
               'verify with lab analysis and BMP tests.')
    if ss.design:
        st.caption(f"Last run: {ss.stamp}")

tabs = st.tabs(['1 Inputs', '2 Kinetics', '3 Performance', '4 Agitator', '5 Heating coil',
                '6 CAD', '7 Scenarios', '8 Report'])

# ---------------------------------------------------------------- 1 Inputs
with tabs[0]:
    st.subheader('Substrate (five required inputs)')
    c1, c2, c3, c4 = st.columns(4)
    Q_tpd = c1.number_input('Substrate flow rate (t/d, wet)', 1.0, 5000.0, 100.0, 5.0, key='Q_tpd')
    T_feed = c2.number_input('Feed temperature (°C)', -5.0, 80.0, 25.0, 1.0, key='T_feed')
    TS = c3.number_input('Total solids, TS (% wet)', 0.5, 40.0, 8.0, 0.5, key='TS')
    VS = c4.number_input('Volatile solids, VS (% of TS)', 10.0, 99.0, 85.0, 1.0, key='VS')

    st.markdown('**Feedstock blend** (% of VS; normalised to 100 %)')
    defaults = {'Napier grass': 40.0, 'Cow dung': 40.0, 'Chicken litter': 20.0}
    cols = st.columns(len(FEED))
    for col, k in zip(cols, FEED):
        col.number_input(k, 0.0, 100.0, defaults.get(k, 0.0), 5.0, key=f'blend_{k}')
    tot = sum(pct_inputs().values())
    if tot <= 0:
        st.error('Enter at least one feedstock above 0 %.')
    elif abs(tot - 100) > 0.01:
        st.info(f'Blend entries sum to {tot:.0f} %; they will be normalised to 100 %.')

    st.subheader('Design basis')
    d1, d2, d3, d4 = st.columns(4)
    HRT = d1.number_input('HRT (d)', 5.0, 90.0, 20.0, 1.0, key='HRT')
    T_op = d2.number_input('Operating temperature (°C)', 25.0, 60.0, 37.0, 1.0, key='T_op')
    T_amb = d3.number_input('Ambient temperature (°C)', -20.0, 45.0, 15.0, 1.0, key='T_amb')
    ins_t = d4.number_input('Insulation thickness (mm)', 0.0, 300.0, 100.0, 10.0, key='ins') / 1000
    e1, e2, e3, e4 = st.columns(4)
    n_tanks = e1.number_input('Number of digesters', 1, 8, 1, key='ntanks')
    HD = e2.number_input('Liquid depth / diameter', 0.5, 2.0, 1.0, 0.1, key='HD')
    alk = e3.number_input('Feed net cations (kmol/m³)', 0.0, 0.2, 0.045, 0.005, format='%.3f', key='alk',
                          help='Net alkalinity from dilution water / buffer. Higher = more pH buffering.')
    kset_name = e4.selectbox('Kinetics set', list(KIN), format_func=lambda k: KIN[k]['label'], key='kset')
    with st.expander('Advanced (mixing, heating, overrides)'):
        a1, a2, a3, a4 = st.columns(4)
        eps = a1.number_input('Specific mixing power (W/m³)', 2.0, 20.0, 6.0, 0.5, key='eps')
        mu_ov = a2.number_input('Viscosity override (Pa·s, 0 = auto)', 0.0, 10.0, 0.0, 0.01, key='mu')
        Tw_in = a3.number_input('Hot-water supply (°C)', 40.0, 95.0, 70.0, 1.0, key='twi')
        Tw_out = a4.number_input('Hot-water return (°C)', 35.0, 90.0, 55.0, 1.0, key='two')
        U_ov = st.number_input('Coil U override (W/m²K, 0 = auto)', 0.0, 800.0, 0.0, 10.0, key='uov')

    # live summary (no ODE solve)
    blend_now = {k: v for k, v in pct_inputs().items() if v > 0}
    if blend_now:
        fi = FeedInputs(Q_tpd, T_feed, TS, VS, blend_now, alk)
        fp = blend_properties(fi, FEED)
        st.subheader('Live summary')
        m1, m2, m3, m4, m5 = st.columns(5)
        m1.metric('Feed volume (m³/d)', f"{fp['Q_m3d']:.0f}")
        m2.metric('VS load (t VS/d)', f"{fp['VS_load_kg_d'] / 1000:.2f}")
        m3.metric('Feed COD (kg/m³)', f"{fp['COD_kg_m3']:.0f}")
        m4.metric('Digester volume (m³)', f"{fp['Q_m3d'] * HRT:,.0f}")
        m5.metric('OLR (kg VS/m³/d)', f"{fp['VS_load_kg_d'] / (fp['Q_m3d'] * HRT):.2f}")
        st.caption(f"Library-predicted blend: TS {fp['TS_lib']:.0f} % as received, VS {fp['VS_lib']:.0f} % of TS "
                   f"(informational). Your entered TS/VS ({TS:g} % / {VS:g} %) drive the design; the blend sets "
                   f"composition, COD/VS, biodegradability and hydrolysis rates.")
    st.caption(feedstock_meta())

# ---------------------------------------------------------------- 2 Kinetics
with tabs[1]:
    ks = KIN[kset_name]
    st.subheader('Monod kinetic constants')
    st.caption(f"Source: {ks['reference']}")
    edit_df = pd.DataFrame([{'Group': v['label'], 'k_m (1/d)': v['k_m'], 'K_S (kg COD/m³)': v['K_S'],
                             'Y': v['Y'], 'b (1/d)': v['b']} for v in ks['uptake'].values()])
    edited = st.data_editor(edit_df, hide_index=True, disabled=['Group'], key='kin_editor',
                            use_container_width=True)
    theta = st.number_input('Temperature coefficient θ (rate × θ^(T−T_ref))', 1.0, 1.15, float(ks['theta']), 0.01,
                            key='theta')
    kt = kinetics_table(ks, T_op, theta)
    st.markdown(f'**Effective values at {T_op:g} °C** (edit the table above to change the base constants)')
    st.dataframe(kt, hide_index=True, use_container_width=True)
    overrides = {}
    for key, (_, row) in zip(ks['uptake'], edited.iterrows()):
        for col, field in (('k_m (1/d)', 'k_m'), ('K_S (kg COD/m³)', 'K_S'), ('Y', 'Y'), ('b (1/d)', 'b')):
            if abs(row[col] - ks['uptake'][key][field]) > 1e-12:
                overrides[f'uptake.{key}.{field}'] = float(row[col])
    if overrides:
        st.warning(f'{len(overrides)} constant(s) edited from the library values; edits apply on the next run.')

    left, right = st.columns(2)
    fig = go.Figure()
    for lab, (S, r) in monod_curves(ks, T_op, theta).items():
        fig.add_trace(go.Scatter(x=S / S.max() * 5, y=r, name=lab, mode='lines'))
    fig.update_layout(title='Monod uptake rate vs substrate (x-axis in multiples of K_S)', height=380,
                      xaxis_title='S / K_S', yaxis_title='specific uptake rate (1/d)', margin=dict(t=40))
    left.plotly_chart(fig, use_container_width=True)
    fig = go.Figure()
    for lab, (x, y) in inhibition_curves(ks).items():
        fig.add_trace(go.Scatter(x=x, y=y, name=lab, mode='lines'))
    fig.update_layout(title='Inhibition factors', height=380, yaxis_title='factor (1 = no inhibition)',
                      margin=dict(t=40))
    right.plotly_chart(fig, use_container_width=True)
    if ss.design:
        h = ss.design['feed']['hyd']
        st.markdown(f"**Blend hydrolysis constants in the current design** (feedstock-weighted): "
                    f"carbohydrate {h['k_hyd_ch']:.2f}, protein {h['k_hyd_pr']:.2f}, lipid {h['k_hyd_li']:.2f} 1/d "
                    f"(the library default of 10 1/d is for sewage sludge, not fibrous feedstocks).")
    with st.expander('Feedstock library'):
        st.dataframe(pd.DataFrame(FEED).T, use_container_width=True)
    st.download_button('Download kinetics library (JSON)',
                       json.dumps({'sets': KIN, 'feedstocks': FEED}, indent=2), 'kinetics_feedstock_library.json')

# ---------------------------------------------------------------- run
if run_clicked:
    if sum(pct_inputs().values()) <= 0:
        st.sidebar.error('Enter a blend first.')
    else:
        di = DesignInputs(
            Q_tpd=Q_tpd, T_feed_C=T_feed, VS_pct_TS=VS, TS_pct=TS, blend={k: v for k, v in pct_inputs().items() if v > 0},
            HRT_d=HRT, T_op_C=T_op, T_amb_C=T_amb, insul_t_m=ins_t, n_tanks=int(n_tanks), H_over_D=HD,
            alk_cat_kmol_m3=alk, kinetics_set=kset_name, theta=theta, eps_W_m3=eps,
            mu_override_Pa_s=mu_ov or None, T_w_in=Tw_in, T_w_out=Tw_out, U_coil_override=U_ov or None,
            k_overrides=overrides)
        with st.spinner('Solving ADM1 CSTR to steady state...'):
            try:
                ss.design = run_design(di, KIN, FEED)
                ss.scenarios = None
                ss.stamp = pd.Timestamp.now().strftime('%Y-%m-%d %H:%M')
                ss.pop('report_bytes', None)
                st.sidebar.success('Design complete')
            except Exception as e:
                st.sidebar.error(f'Design failed: {e}')

D = ss.design


def need_run():
    st.info('Enter the inputs on tab 1 and press **Run design** in the sidebar.')


# ---------------------------------------------------------------- 3 Performance
with tabs[2]:
    if not D:
        need_run()
    else:
        p = D['performance']
        if D['flags']:
            for fl in D['flags']:
                st.warning(fl)
        else:
            st.success('No process flags raised at these conditions.')
        k = st.columns(5)
        k[0].metric('Biogas (Nm³/d)', f"{p['biogas_Nm3_d']:,.0f}")
        k[1].metric('Methane (Nm³/d)', f"{p['CH4_Nm3_d']:,.0f}", f"{p['CH4_pct']:.1f} % CH₄")
        k[2].metric('Methane yield (Nm³/kg VS)', f"{p['SMY_Nm3_kgVS']:.3f}", f"BMP {p['BMP_blend']:.3f}", delta_color='off')
        k[3].metric('pH', f"{p['pH']:.2f}")
        k[4].metric('VFA (kg COD/m³)', f"{p['VFA_kgCOD_m3']:.2f}")
        k = st.columns(5)
        k[0].metric('Volume (m³)', f"{D['spec'].V_liq:,.0f}")
        k[1].metric('OLR (kg VS/m³/d)', f"{p['OLR_kgVS_m3_d']:.2f}")
        k[2].metric('Free NH₃-N (mg/L)', f"{p['NH3_mgN_L']:.0f}")
        k[3].metric('COD removal (%)', f"{p['COD_removal_pct']:.0f}")
        k[4].metric('Net energy (kW)', f"{D['energy']['net_kW']:.0f}", f"{D['energy']['parasitic_pct']:.0f} % parasitic", delta_color='off')
        h = D['history']
        c1, c2, c3 = st.columns(3)
        for col, key, title in ((c1, 'ch4', 'Methane (Nm³/d)'), (c2, 'pH', 'pH'), (c3, 'vfa', 'VFA (kg COD/m³)')):
            f = go.Figure(go.Scatter(x=h['t'], y=h[key], mode='lines'))
            f.update_layout(title=title, height=280, xaxis_title='days from start-up', margin=dict(t=40, l=10, r=10))
            col.plotly_chart(f, use_container_width=True)
        a, b = st.columns(2)
        with a:
            st.markdown('**Effluent state (kg COD/m³ unless noted)**')
            s = pd.Series(D['state']).rename('value').to_frame()
            st.dataframe(s.style.format('{:.4g}'), height=360, use_container_width=True)
        with b:
            st.markdown('**Feed to ADM1**')
            fr = D['feed']['fractions']
            st.dataframe(pd.DataFrame({'COD fraction': {'Carbohydrate (biodegradable)': fr['ch'], 'Protein': fr['pr'],
                                                       'Lipid': fr['li'], 'Acetate': fr['vfa'], 'Inert / slowly degradable': fr['inert']}}).style.format('{:.1%}'),
                         use_container_width=True)
            st.markdown(f"**COD balance closure:** {100 * p['cod_closure']:.1f} %  |  "
                        f"**Conversion of BMP:** {p['conversion_pct']:.0f} %")
            st.caption('ADM1 does not wash methanogens out until HRT ≈ 3 d, so short-HRT results are optimistic. '
                       'The model reproduces the BSM2 benchmark within 1 % (see tests).')

# ---------------------------------------------------------------- 4 Agitator
with tabs[3]:
    if not D:
        need_run()
    else:
        m = D['mixing']
        for fl in m['flags']:
            st.warning(fl)
        k = st.columns(5)
        k[0].metric('Agitators per tank', m['n_agitators'])
        k[1].metric('Impeller Ø (m)', f"{m['D_impeller_m']:.2f}", f"D/T {m['D_over_T']:.2f}", delta_color='off')
        k[2].metric('Speed (rpm)', f"{m['N_rpm']:.1f}")
        k[3].metric('Shaft power each (kW)', f"{m['P_shaft_kW']:.1f}")
        k[4].metric('Motor, IEC (kW)', f"{m['motor_kW']}")
        st.dataframe(pd.DataFrame({
            'Parameter': ['Specific power', 'Apparent viscosity', 'Reynolds number', 'Tip speed', 'Pumping capacity (each)',
                          'Turnover time', 'Velocity gradient G', 'Torque', 'Shaft diameter (min.)'],
            'Value': [f"{m['eps_W_m3']:.1f} W/m³", f"{m['mu_Pa_s']:.3f} Pa·s", f"{m['Re']:,.0f}", f"{m['tip_speed_m_s']:.2f} m/s",
                      f"{m['pumping_m3_h']:,.0f} m³/h", f"{m['turnover_min']:.1f} min", f"{m['G_1_s']:.1f} 1/s",
                      f"{m['torque_Nm']:,.0f} N·m", f"{m['shaft_diameter_mm']:.0f} mm"]}), hide_index=True)
        st.caption('Method: top-entry axial hydrofoil, N = (P/(Np·ρ·D⁵))^(1/3) with Np = 0.30, Nq = 0.55. '
                   'Viscosity is an indicative function of TS. Confirm with the agitator vendor.')

# ---------------------------------------------------------------- 5 Heating coil
with tabs[4]:
    if not D:
        need_run()
    else:
        ht, cl = D['heat'], D['coil']
        for fl in cl['flags']:
            st.warning(fl)
        k = st.columns(4)
        k[0].metric('Design duty per tank (kW)', f"{ht['Q_design_kW']:.0f}")
        k[1].metric('Coil area (m²)', f"{cl['area_m2']:.1f}")
        k[2].metric('Pipe length (m)', f"{cl['length_m']:.0f}", f"{cl['n_circuits']} circuit(s)", delta_color='off')
        k[3].metric('Water flow (m³/h)', f"{cl['water_flow_m3_h']:.1f}", f"{cl['water_velocity_m_s']:.2f} m/s", delta_color='off')
        a, b = st.columns(2)
        f = go.Figure(go.Bar(x=[ht['Q_feed_kW'], ht['Q_wall_kW'], ht['Q_roof_kW'], ht['Q_floor_kW']],
                             y=['Feed heating', 'Wall loss', 'Roof loss', 'Floor loss'], orientation='h'))
        f.update_layout(title='Heat load (kW, before margin)', height=300, margin=dict(t=40))
        a.plotly_chart(f, use_container_width=True)
        b.dataframe(pd.DataFrame({'Parameter': ['Water supply / return', 'U (indicative)', 'LMTD', 'Pipe OD / ID',
                                                'Coil diameter', 'Total turns', 'Pitch', 'Margin'],
                                  'Value': [f"{cl['T_w_in']:.0f} / {cl['T_w_out']:.0f} °C", f"{cl['U_W_m2K']:.0f} W/m²K",
                                            f"{cl['LMTD_K']:.1f} K", f"{cl['OD_mm']:.1f} / {cl['ID_mm']:.1f} mm",
                                            f"{cl['coil_diameter_m']:.2f} m", f"{cl['turns_total']:.1f}",
                                            f"{cl['pitch_m']:.2f} m", f"{100 * ht['margin']:.0f} %"]}), hide_index=True)

# ---------------------------------------------------------------- 6 CAD
with tabs[5]:
    if not D:
        need_run()
    else:
        cp = cad_params(D)
        parts = build_parts(cp)
        left, right = st.columns([3, 1])
        left.plotly_chart(plotly_figure(parts), use_container_width=True, height=620)
        with right:
            g = D['geometry']
            st.dataframe(pd.DataFrame({'Dimension': ['Diameter', 'Liquid depth', 'Free wall', 'Dome rise', 'Total height',
                                                     'Wall thickness', 'Impeller Ø', 'Coil Ø'],
                                       'm': [g['D'], g['H_liq'], g['H_free'], g['h_dome'], g['H_total'], g['wall_t'],
                                             cp['D_imp'], cp['coil_D']]}).style.format({'m': '{:.2f}'}), hide_index=True)
            st.download_button('FreeCAD script (.py)', make_script(cp), 'build_cstr_freecad.py',
                               help='Run with freecadcmd, or paste into FreeCAD / your FreeCAD MCP execute-code tool.')
            st.download_button('Preview mesh (.stl)', to_stl_bytes(parts), 'cstr_preview.stl')
            if freecad_executable():
                if st.button('Build FCStd / STEP with FreeCAD'):
                    out = run_freecad(cp, tempfile.mkdtemp(prefix='cstr_cad_'))
                    st.success('FreeCAD model built') if out['ok'] else st.error(out['message'])
                    for name, path in (out.get('files') or {}).items():
                        st.download_button(name, path.read_bytes(), name, key=f'dl_{name}')
            else:
                st.caption('FreeCAD not found on this host. Use the script above in your FreeCAD session / MCP; '
                           'the viewer shows an independent preview built from the same dimensions.')

# ---------------------------------------------------------------- 7 Scenarios
with tabs[6]:
    if not D:
        need_run()
    else:
        st.markdown('How can the system be modified? Pick a lever and compare outcomes against the current design.')
        kind = st.selectbox('Lever', ['HRT', 'Temperature', 'TS', 'Blend'])
        base = DesignInputs(**{k: v for k, v in D['inputs'].items()})
        if kind == 'HRT':
            vals = st.text_input('HRT values (d)', '10, 15, 20, 30, 40')
            vals = [float(x) for x in vals.split(',') if x.strip()]
        elif kind == 'Temperature':
            vals = st.text_input('Temperatures (°C)', '30, 35, 37, 45, 55')
            vals = [float(x) for x in vals.split(',') if x.strip()]
        elif kind == 'TS':
            vals = st.text_input('TS values (%)', '4, 6, 8, 10, 12')
            vals = [float(x) for x in vals.split(',') if x.strip()]
        else:
            vals = [(f'{k} only', {k: 100.0}) for k in FEED] + [('Current blend', D['inputs']['blend'])]
        if st.button('Run scenarios'):
            with st.spinner('Running scenarios...'):
                ss.scenarios = sweep(base, kind, vals, kinetics=KIN, feedstocks=FEED)
        if ss.scenarios is not None:
            df = ss.scenarios
            st.dataframe(df, hide_index=True, use_container_width=True)
            if 'CH4_Nm3_d' in df:
                f = go.Figure()
                f.add_trace(go.Bar(x=df['Scenario'], y=df['CH4_Nm3_d'], name='CH₄ (Nm³/d)'))
                f.add_trace(go.Scatter(x=df['Scenario'], y=df['Net_kW'], name='Net energy (kW)', yaxis='y2', mode='lines+markers'))
                f.update_layout(height=360, yaxis2=dict(overlaying='y', side='right'), margin=dict(t=20))
                st.plotly_chart(f, use_container_width=True)
            st.caption('"No_flags = False" means at least one process flag was raised for that case (see tab 3 logic).')

# ---------------------------------------------------------------- 8 Report
with tabs[7]:
    if not D:
        need_run()
    else:
        st.markdown('**Report contents:** design summary and flags · inputs and blend · kinetics (with references) · '
                    'process performance · agitator sizing · heat load and coil · energy balance · 3D model image · '
                    'scenario comparison (if run on tab 7) · validation and limitations.')
        if st.button('Prepare design report', type='primary'):
            with st.spinner('Building Word report...'):
                tmp = Path(tempfile.mkdtemp()) / 'cad.png'
                render_png(build_parts(cad_params(D)), str(tmp))
                ss.report_bytes = build_report(D, str(tmp), ss.scenarios, project)
        if ss.get('report_bytes'):
            st.download_button('Download design report (.docx)', ss.report_bytes,
                               f"{project.replace(' ', '_')}_design_report.docx",
                               'application/vnd.openxmlformats-officedocument.wordprocessingml.document')
