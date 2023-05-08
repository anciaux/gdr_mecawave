#!/usr/bin/env python
# coding: utf-8
################################################################
import matplotlib.pyplot as plt
import tuto
import numpy as np
import akantu as aka
import streamlit as st
import solidipes as sp
import os
################################################################
from app_helper import params_selector, make_figure, tqdm
################################################################
dirname = os.path.dirname(__file__)

params = {
    'L': 10,
    'l': 10,
    'h1': .05,
    'h2': .5,
    'U': .1,
}
col1, col2 = st.columns(2)

with col1:
    f = sp.load_file(os.path.join(dirname, 'schematic.png'))
    f.view()


with col2:
    params_selector(params, dirname)

    compute_variations = st.checkbox('Compute K variations!')
    if compute_variations is False:
        zoom_factor = st.number_input('Zoom factor', value=1)
        zoom_range = params['L']/zoom_factor

crack_length = params['l']

button = st.button('Compute!', use_container_width=True, type='primary')


def on_click(**params):
    if compute_variations:
        variations_view(**params)
    else:
        fields_view(**params)


def fields_view(**params):
    with st.spinner("Creating model and solving"):
        model, mesh = tuto.createModel(**params)
        fig = tuto.plotMesh(mesh, displacement=model.getDisplacement())
        st.pyplot(fig)
        fig = tuto.plotResult(model, displacement=model.getDisplacement(),
                              field='displacement', contour=5,
                              xrange=[-zoom_range+crack_length,
                                      zoom_range+crack_length],
                              yrange=[-zoom_range, zoom_range])
        st.pyplot(fig)

        fig = tuto.plotResult(
            model, displacement=model.getDisplacement(), field='stress',
            xrange=[-zoom_range+crack_length,
                    zoom_range+crack_length],
            yrange=[-zoom_range, zoom_range])
        st.pyplot(fig)

    ################################################################
    with st.spinner("Plot stress"):
        stress = model.getMaterial(0).getStress(aka._triangle_3)
        stress_norm = np.linalg.norm(stress, axis=1)
        stress_max = stress_norm.max()

        c = np.logspace(np.log10(stress_max/7), np.log10(stress_max), 10)
        fig = tuto.plotResult(model, displacement=model.getDisplacement(),
                              field='stress', contour=c, log_color=True,
                              xrange=[-zoom_range+crack_length,
                                      zoom_range+crack_length],
                              yrange=[-zoom_range, zoom_range])
        st.pyplot(fig)
    ################################################################
    with st.spinner(r"Fit williams from $\sigma_11$"):
        curve, K = tuto.extract_stress(
            model, theta=0, r_fit=params['h2'], **params)
        r = curve[:, 0]
        sigma = curve[:, 1]
        st.markdown(f"# Found Stress Intensity Factor $K={K}$")
        with make_figure() as (fig, axe):
            axe.plot(r, sigma, '-', label='FE')
            axe.plot(r, K/np.sqrt(r*2*np.pi), '--',
                     label=f'Williams $K^I = {K}$')
            axe.set_xlabel(r'$\quad[m]$')
            axe.set_ylabel(r'$\sigma_{rr}\quad[Pa]$')
            axe.legend(loc='best')

    ################################################################
    with st.spinner(r"Fit williams from $\sigma_22$"):
        curve, K = tuto.extract_stress(
            model, theta=0, r_fit=.2, shear=True, **params)
        st.markdown(f"# Found Stress Intensity Factor $K={K}$")
        r = curve[:, 0]
        sigma = curve[:, 1]
        with make_figure() as (fig, axe):
            axe.plot(r, sigma, '-', label='FE')
            axe.plot(r, K/np.sqrt(r*2*np.pi), '--',
                     label=f'Williams ($K^I = {K}$)')
            axe.set_xlabel(r'$\quad[m]$')
            axe.set_ylabel(r'$\sigma_{\theta\theta}\quad[Pa]$')
            axe.legend(loc='best')


def variations_view(**params):
    with st.spinner("Creating model and solving"):
        model, mesh = tuto.createModel(**params)

    res = []

    for l in tqdm(np.arange(1, 10), init_text="Varying crack length"):
        params = {
            'L': 10,
            'l': l,
            'h1': .0005,
            'h2': .2,
            'U': .1,
        }
        K = tuto.full_study_extract_K(quiet=True, **params)
        res.append((l, K))

    res = np.array(res)
    with make_figure() as (fig, axe):
        axe.plot(res[:, 0], res[:, 1], 'o-')
        axe.set_xlabel('Crack length [m]')
        _ = axe.set_ylabel('$K^I$')

    ################################################################

    res = []

    for h1 in tqdm([.05, .01, .005, .001, .0005],
                   init_text="Varying refinement"):
        params = {
            'L': 10,
            'l': 10,
            'h1': h1,
            'h2': .2,
            'U': .1,
        }
        K = tuto.full_study_extract_K(quiet=True, **params)
        res.append((h1, K))

    res = np.array(res)
    with make_figure() as (fig, axe):
        axe.plot(res[:, 0], res[:, 1], 'o-')
        axe.set_xlabel('Mesh characteristic size [m]')
        _ = axe.set_ylabel('$K^I$')


################################################################
if button:
    on_click(**params)
