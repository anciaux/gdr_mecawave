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


class tqdm:

    def __init__(self, container, init_text="Please wait"):
        self.container = container
        self.text = init_text
        self.bar = st.progress(0, init_text)

    class _iter:
        def __init__(self, cont):
            self.cont = cont
            self.n = len(cont.container)
            self.current = 0
            self._internal_iter = self.cont.container.__iter__()
            self.cont.bar.progress(float(self.current/self.n),
                                   text=f"{self.current}/{self.n}")

        def __next__(self):
            self.current += 1
            res = self._internal_iter.__next__()
            self.cont.bar.progress(
                float(self.current/self.n),
                text=f"{self.cont.text}: {res} ({self.current}/{self.n})")
            return res

    def __iter__(self):
        return tqdm._iter(self)


################################################################
dirname = os.path.dirname(__file__)

params = {
    'L': 10,
    'crack_length': 10,
    'h1': .05,
    'h2': .5,
    'U': .1,
}
col1, col2 = st.columns(2)

with col1:
    f = sp.load_file(os.path.join(dirname, 'schematic.svg'))
    f.view()


for p, value in params.items():
    _type = type(value)
    val = col2.text_input(p, value=value)
    params[p] = _type(val)

button = st.button('Compute!')

if button:
    with st.spinner("Creating model and solving"):
        model, mesh = tuto.createModel(**params)
        fig = tuto.plotMesh(mesh, displacement=model.getDisplacement())
        st.pyplot(fig)
        fig = tuto.plotResult(model, displacement=model.getDisplacement(),
                              field='displacement', contour=5)
        st.pyplot(fig)

        fig = tuto.plotResult(
            model, displacement=model.getDisplacement(), field='stress')
        st.pyplot(fig)
    ################################################################
    with st.spinner("Plot stress"):
        stress = model.getMaterial(0).getStress(aka._triangle_3)
        stress_norm = np.linalg.norm(stress, axis=1)
        stress_max = stress_norm.max()
        c = np.logspace(np.log10(stress_max/7), np.log10(stress_max), 10)
        crack_length = params['crack_length']
        zoom_range = params['h1']*10
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
        fig = plt.figure()
        axe = fig.add_subplot(111)
        axe.plot(r, sigma, '-', label='FE')
        axe.plot(r, K/np.sqrt(r*2*np.pi), '--',
                 label=f'Williams $K^I = {K}$')
        axe.set_xlabel(r'$\quad[m]$')
        axe.set_ylabel(r'$\sigma_{rr}\quad[Pa]$')
        axe.legend(loc='best')
        st.pyplot(fig)

    ################################################################
    with st.spinner(r"Fit williams from $\sigma_22$"):
        curve, K = tuto.extract_stress(
            model, theta=0, r_fit=.2, shear=True, **params)
        st.markdown(f"# Found Stress Intensity Factor $K={K}$")
        r = curve[:, 0]
        sigma = curve[:, 1]
        fig = plt.figure()
        axe = fig.add_subplot(111)
        axe.plot(r, sigma, '-', label='FE')
        axe.plot(r, K/np.sqrt(r*2*np.pi), '--',
                 label=f'Williams ($K^I = {K}$)')
        axe.set_xlabel(r'$\quad[m]$')
        axe.set_ylabel(r'$\sigma_{\theta\theta}\quad[Pa]$')
        axe.legend(loc='best')
        st.pyplot(fig)

    ################################################################
    res = []

    for l in tqdm(np.arange(1, 10), init_text="Varying crack length"):
        params = {
            'L': 10,
            'crack_length': l,
            'h1': .0005,
            'h2': .2,
            'U': .1,
        }
        K = tuto.full_study_extract_K(quiet=True, **params)
        res.append((l, K))

    res = np.array(res)
    fig = plt.figure()
    axe = fig.add_subplot(111)
    axe.plot(res[:, 0], res[:, 1], 'o-')
    axe.set_xlabel('Crack length [m]')
    _ = axe.set_ylabel('$K^I$')
    st.pyplot(fig)

    ################################################################

    res = []

    for h1 in tqdm([.05, .01, .005, .001, .0005], init_text="Varying refinement"):
        params = {
            'L': 10,
            'crack_length': 10,
            'h1': h1,
            'h2': .2,
            'U': .1,
        }
        K = tuto.full_study_extract_K(quiet=True, **params)
        res.append((h1, K))

    res = np.array(res)
    fig = plt.figure()
    axe = fig.add_subplot(111)
    axe.plot(res[:, 0], res[:, 1], 'o-')
    axe.set_xlabel('Mesh characteristic size [m]')
    _ = axe.set_ylabel('$K^I$')
    st.pyplot(fig)

    ################################################################
