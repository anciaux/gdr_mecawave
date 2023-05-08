import streamlit as st
import os

lectures = ["0_Introduction_class",
            "1_Course_Dynamics",
            "2_Course_Fracture_Mechanics",
            "3_Course_Fragmentation",
            "4_Course_fe_dynamics",
            "5_1_Course_Comp_Frac_Cohesive",
            "5_2_Course_Comp_Frac_Phase_Field",
            "6_Course_HPC_Fragmentation"]

for f in lectures:
    fname = f'slides/{f}/{f}.html'
    content = open(fname, 'r').read()
    st.markdown(content, unsafe_allow_html=True)
