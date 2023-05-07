import os
import sys
import streamlit as st
import matplotlib.pyplot as plt
################################################################
plt.rcParams['figure.figsize'] = [10, 10]


sub_apps = ['stress_intensity_factor',
            'plate_hole']

tabs = st.tabs(sub_apps)

for name, tab in zip(sub_apps, tabs):
    with tab:
        try:
            import importlib
            sys.path.append(name)
            spec = importlib.util.spec_from_file_location(
                name, os.path.join(name, name+'.py'))
            foo = importlib.util.module_from_spec(spec)
            sys.modules[name] = foo
            spec.loader.exec_module(foo)
        except Exception as err:
            st.error(err)
            raise err
