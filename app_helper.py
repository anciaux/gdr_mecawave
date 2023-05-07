import streamlit as st
from contextlib import contextmanager
import matplotlib.pyplot as plt

################################################################


@contextmanager
def make_figure(equal_ratio=False):

    fig = plt.figure()
    axe = fig.add_subplot(111)
    if equal_ratio:
        axe.set_aspect('equal')
    try:
        yield fig, axe
    finally:
        st.pyplot(fig)


################################################################

def params_selector(params, key):
    with st.expander('Model parameters'):
        for p, value in params.items():
            _type = type(value)
            val = st.text_input(p, value=value, key=p+key)
            params[p] = _type(val)


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
