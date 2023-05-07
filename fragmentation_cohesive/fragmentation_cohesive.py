#!/usr/bin/env python
# coding: utf-8

################################################################
import main
import numpy as np
import matplotlib.pyplot as plt
from app_helper import params_selector, make_figure, tqdm
import os
import streamlit as st
import solidipes as sp
################################################################

params = {
    'young_modulus':  275.0 * 10**9,  # (Pa)
    'density':  2750.0,  # (kg/m3)
    'fracture_energy':  100.0,  # (N/m)
    'stress_limit':  300.0 * 10**6,  # (Pa)
    'bar_length': 50 * 10**-3,  # (m)
    'number_elements':  500 * 2,  # Total number of triangular elements
    'area':  1.0,  # Cross sectional area (m2) (Equal to element size )
    # Load
    'strain_rate':  10.0**4,  # (s-1)
    # Time
    'time_simulation':  5 * 10**-7,  # Total time of simulation (s)
    'dump_freq': 1,
    'variation_sigma_c': 1e6,  # random(uniform) variation of sigma_c
    'contact': 0
}

dirname = os.path.dirname(__file__)


col1, col2 = st.columns(2)

with col1:
    # f = sp.load_file(os.path.join(dirname, 'plate-hole-2.svg'))
    # f.view()
    pass

with col2:
    params_selector(params, dirname)

button = st.button('Compute!', key='button'+dirname,
                   use_container_width=True, type='primary')


params.update({
    # Left extremitiy x coordinate / 0-initial
    'x0': -0.5 * params['bar_length'],
    # Rigth extremitiy x coordinate / f-final
    'xf':  0.5 * params['bar_length'],
    'paraview': False,
})


if button:
    model, mesh = main.generate_model.create_model(**params)
    results, energies = main.runSimulation(model, **params)

    n = results['n_fragments'].iloc[-1]
    st.markdown(f'### Created {n} fragments')
    main.analyze_results(model, results, energies, bins=10, **params)
