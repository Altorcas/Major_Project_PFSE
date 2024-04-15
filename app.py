import streamlit as st
import plotly.graph_objects as go
import base_plates as bp
from handcalcs.decorator import handcalc
from math import pi, sqrt

with st.expander("Assumptions", expanded=True):
    st.write("-It is assumed that there is no foundation joint material, so that the base plate is in direct contact with the concrete foundation")
    st.write("-When there is only compression, stiffeners have not been considered")
    st.write("-When bending moment produces tenssion stresses, base plate with stiffeners is considered:")
    st.image("Baseplate_withstiff.jpg", caption='Base plate with stiffeners', width=200)

st.markdown("# Base plate verification")

st.sidebar.subheader("Input data- Cross section")
h = st.sidebar.number_input("Depth (mm)", value=400)
b = st.sidebar.number_input("Width (mm)", value=180)
tw = st.sidebar.number_input("Web thickness (mm)", value=8.6)
tf = st.sidebar.number_input("Flange thickness (mm)", value=13.5)
fyk = st.sidebar.number_input("Steel grade (MPa)", value=275)
phi0 = st.sidebar.number_input("Steel safety factor", value=1.05)
d_edge = st.sidebar.number_input("Flange edge to base plate edge distance (mm)", value=100)

st.sidebar.subheader("Input data- Base plate")
H = st.sidebar.number_input("Base plate depth (mm)", value=600)
B = st.sidebar.number_input("Base plate width (mm)", value=400)
t = st.sidebar.number_input("Thickness (mm)", value=22)
fck = st.sidebar.number_input("Characteristic compression resistante of concrete (MPa)", value=25)
phic = st.sidebar.number_input("Concrete safety factor", value=1.5)

st.sidebar.subheader("Bolts")
D_bolt = st.sidebar.number_input("Bolt diameter (mm)", value=20)
L_bolt = st.sidebar.number_input("Bolt length (mm)", value=80)
d_bolt = st.sidebar.number_input("Distance from the edge of the base plate to the center of the bolt (mm)", value=50)
fyk_bolt = st.sidebar.number_input("Bolt steel grade (MPa)", value=275)
phi_bolt = st.sidebar.number_input("Bolt safety factor (MPa)", value=1.15)
E_bolt = st.sidebar.number_input("Bolt Yound modulus (MPa)", value=210000)


st.sidebar.subheader("Stiffener")
h_stiff = st.sidebar.number_input("Stiffener height (mm)", value=150)
t_stiff = st.sidebar.number_input("Stiffener thickness (mm)", value=10)
fyk_stiff = st.sidebar.number_input("Stiffener grade steel (MPa)", value=275)

st.sidebar.subheader("Design values")
NEd = st.sidebar.number_input("Design axil value (N)", value=600)
MEd = st.sidebar.number_input("Design bending moment value (N·mm)", value=400)


plate = bp.baseplate(h, b, tw, tf, fyk, phi0, d_edge, H, B, t, fck, phic, h_stiff, t_stiff, fyk_stiff, D_bolt, L_bolt, d_bolt, fyk_bolt, phi_bolt, E_bolt)

with st.expander("Base plate characterization", expanded=True):

    c = plate.additional_bearing_width()
    st.write(f"The additional bearing widht (c) has a value of {c: .2f} (mm)")
    c_latex, c_value = bp.calc_additional_bearing_width(plate.t,plate.fyk,plate.phi0,plate.fck,plate.phic)
    st.latex(c_latex)

    if c > (H / 2 - h / 2) and c > (B / 2 - b / 2):
        st.image("CASE1.jpg", caption='c > (H / 2 - h / 2) and c > (B / 2 - b / 2)', width=350)
    
    elif c > (H / 2 - h / 2) and c < (B / 2 - b / 2):
        st.image("CASE2.jpg", caption='c > (H / 2 - h / 2) and c < (B / 2 - b / 2)', width=350)
    
    elif c < (H / 2 - h / 2) and c > (B / 2 - b / 2):
        st.image("CASE3.jpg", caption='c < (H / 2 - h / 2) and c > (B / 2 - b / 2)', width=350)

    elif c < (H / 2 - h / 2) and c < (B / 2 - b / 2):
        st.image("CASE4.jpg", caption='c < (H / 2 - h / 2) and c < (B / 2 - b / 2)', width=350)

    A = plate.effective_area()
    st.write(f"The effective area has a value of {A: .2f} ($mm^2$)")
    aeff_latex, aeff_value = bp.calc_effective_area(plate.h, plate.b, plate.tw, plate.tf, plate.fyk, plate.phi0, plate.H, plate.B, plate.t, plate.fck, plate.phic)
    if aeff_latex:
        st.latex(aeff_latex)

    Wel = plate.Wel_effective_baseplate()
    st.write(f"The effective elastic section modulus of the base plate has a value of {Wel: .2f} ($mm^3$)")

    st.subheader("Design value action ratio VS resistance ratio")

    ratios_latex, ratios = bp.ratios(MEd, NEd, aeff_value, Wel)
    st.latex(ratios_latex)

    st.subheader("Characterization")

    if MEd/NEd <= Wel/aeff_value:
        st.latex("Base\;plate\;only\;in\;compression: M_{Ed}/N_{Ed}<W_p/A_p")
        st.image("CompressionOnly.jpg", caption='Compression only', width=200)

    if MEd/NEd > Wel/aeff_value:
        st.latex("Bending\;moment\;on\;the\;base\;plate\;exists: M_{Ed}/N_{Ed}>W_p/A_p")
        st.image("Compression + Bending.jpg", caption='Compression+Bending moment', width=200)

if MEd/NEd <= Wel/aeff_value:
    with st.expander("SLU base plate: Compression only", expanded=True):

        f_jd_latex, f_jd_value = bp.joint_design_bearing_strength(fck, phic)
        st.latex(f_jd_latex)

        stres_latex, stress_value = bp.check_compression(NEd, A, MEd, Wel)
        st.latex(stres_latex)

        if stress_value <= f_jd_value:
            st.write("Compression ULS for base plate is satisfied")

        if stress_value > f_jd_value:
            st.write("Compression ULS for base plate is NOT satisfied")
            st.write("Advice:")
            st.write("  -Increase base plate thickness")
            st.write("  -Include base plate stiffeners")
        
if MEd/NEd > Wel/aeff_value:
    with st.expander("ULS base plate: Compression + Bending (stiffeners included)", expanded=True):

        st.image("Compression+Bending+Stiffeners.jpg", caption='Compression+Bending+Stiffeners', width=300)
        st.image("Bearing surface in Compression + Bending.jpg", caption='Bearing surface in Compression + Bending', width=400)
        bearing_b_latex, bearing_b_value = bp.bearing_width(t_stiff, b, B, c)
        st.latex(bearing_b_latex)


        st.image("Scheme 1.jpg", caption='Scheme of actions and structural behaviour', width=250)
        y1, y2 = bp.compression_width(H, d_bolt, MEd, NEd, bearing_b_value, fck, phic)

        st.latex(r"\Sigma M_A = 0")
        st.latex(r"M_{Ed} + N_{Ed} \cdot (H/2-d_{bolt})=b_{bearing} \cdot y \cdot f_{jd} \cdot (H-d_{bolt}-y/2)")
        st.latex(r"y^2 + 2 \cdot (d_{bolt}-H) + 2 \cdot [(M_{Ed} + N_{Ed} \cdot (H/2-d_{bolt}))/(b_{bearing} \cdot f_{jd})] = 0")
        st.write("y1: Root 1:", y1)
        st.write("y2: Root 2:", y2)

        selected_root = st.radio("Select root which has geometric sense:", ("Root 1", "Root 2"))

        selected_value = y1 if selected_root == "Root 1" else y2

        st.write("Selected value:", round(selected_value, 2))

        st.write(r"Tenssion force $T_{d}$ is obtained as:")
        st.latex(r"\Sigma F_V = 0 ")
        st.latex(r"T_d + N_{Ed} = b_{bearing} \cdot y \cdot f_{jd}")
        Td = bp.bolt_stress(NEd, bearing_b_value, selected_value, fck, phic)
        st.latex(r"Td = " + str(round(Td, 2)) + r" \, \text{N} = " + str(round(Td/1000, 2)) + r" \, \text{kN}")

        
        
        st.write(r"New gravity center $y_{g}$: stiffener + base plate (mm)")
        st.latex(r"y_g = (\Sigma A_i \cdot y_i)/(\Sigma A_i)")
        new_yg_latex, new_yg_value = bp.bp_stiff_gravity_center(h_stiff, t_stiff, t, B)
        st.latex(new_yg_latex)

        
        st.write(r"New moment of inertia $I_{g}$: stiffener + base plate (mm$^4$)")
        st.latex(r"I_g = \Sigma I_i + \Sigma A_i \cdot (y_i-y_g)^2")
        st.latex(r"I_g = I_{stiff}+I_{bp} + d_1 + d_2")
        new_I_latex, new_I_value = bp.bp_stiff_momentinertia(h_stiff, t_stiff, t, B, new_yg_value)
        st.latex(new_I_latex)     

        st.write(r"New $W_{el}$: stiffener + base plate (mm$^3$)")
        new_Wel_latex, new_Wel_value = bp.bp_stiff_wel(new_I_value, new_yg_value, h_stiff, t)
        st.latex(new_Wel_latex)



        st.write("ULS checking")
        st.image("Scheme 2.jpg", caption='Calculation bending moment', width=250)
        st.write(r"Bending moment at section A")
        st.write(r" if: y <= v")
        st.latex(r"M_{A,Ed} = f_{jd} \cdot b_{bearing} \cdot y \cdot (v - y/2)")
        st.write(r" if: y > v")
        st.latex(r"M_{A,Ed} = f_{jd} \cdot b_{bearing} \cdot v^2/2")
        st.write(r"Bending moment at section B")
        st.latex(r"M_{B,Ed} = T_d \cdot (v - b_{bearing})")
        st.write(r"Calculation bending moment")
        st.latex(r"M = max(M_{A,Ed}, M_{B,Ed})")
        Med = bp.calculation_moment(fck, phic, bearing_b_value, selected_value, d_edge, Td, d_bolt)
        st.latex(r"M_{Ed} = " + str(round(Med, 2)) + r" \, \text{Nmm}= " + str(round(Med/1000000, 2)) + r" \, \text{kNm}")
        
        ratio_latex, ratio_value = bp.uls_compression_bp_stiffener(Med, new_Wel_value, fyk, phi0)
        st.latex(ratio_latex)

        if ratio_value <= 1:
            st.write("Satisfies ULS for compression in the base plate")
        if ratio_value > 1:
            st.write("Does not satisfy ULS for compression in the base plate")            

    with st.expander("ULS bolts: Compression + Bending (stiffeners included)", expanded=True):
        st.write("Threaded bolts are considered")
        bp_stress_latex, bp_stress_value = bp.steel_stress_baseplate(selected_value, H, d_bolt, E_bolt)
        st.latex(bp_stress_latex)

        if bp_stress_value <= fyk_bolt/phi0:
            st.write("4 bolts (diameter 20 mm) simetrically located from the centre of the base plate is enough (2 per side)")

        if bp_stress_value > fyk_bolt/phi0:
            st.write("It is necessary to calculate number of bolts and their diameter.")
            st.write("Maximum 3 bolts per side is considered.")

            st.write("The bolt surface area required is:")
            bolt_area_latex, bolt_area_value = bp.bolt_area(Td, fyk_bolt, phi_bolt)
            st.latex(bolt_area_latex)

            st.write("Iteration with different bolts:")

            num_bolt_20mm_latex, num_bolt_20mm_value = bp.number_bolt_20(bolt_area_value)
            st.latex(r"20 \ mm - > N = " +str(round(num_bolt_20mm_value, 2))+ r"- >" + str(round(int(num_bolt_20mm_value)+1, 2)) + r" \ bolts \ per \ side")

            num_bolt_25mm_latex, num_bolt_25mm_value = bp.number_bolt_25(bolt_area_value)
            st.latex(r"25 \ mm - > N = " +str(round(num_bolt_25mm_value, 2))+ r"- >" + str(round(int(num_bolt_25mm_value)+1, 2)) + r" \ bolts \ per \ side")

            num_bolt_30mm_latex, num_bolt_30mm_value = bp.number_bolt_30(bolt_area_value)
            st.latex(r"30 \ mm - > N = " +str(round(num_bolt_30mm_value, 2))+ r"- >" + str(round(int(num_bolt_30mm_value)+1, 2)) + r" \ bolts \ per \ side")

            st.write("For this case, it is necessary to use:")

            if int(num_bolt_20mm_value)+1 <= 3:
                st.write(f"-{max(int(num_bolt_20mm_value)+1,2)} bolts, diameter 20 (mm) per side")
             
            if int(num_bolt_20mm_value)+1 > 3 and int(num_bolt_25mm_value)+1 <= 3:
                st.write(f"-{max(int(num_bolt_25mm_value)+1,2)} bolts, diameter 25 (mm) per side")
            
            if int(num_bolt_20mm_value)+1 > 3 and int(num_bolt_25mm_value)+1 > 3 and int(num_bolt_30mm_value)+1 <= 3:
                st.write(f"-{max(int(num_bolt_30mm_value)+1,2)} bolts, diameter 30 (mm) per side")
        
                    