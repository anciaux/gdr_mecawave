import os
import sys
import streamlit as st
################################################################


sub_apps = [e for e in os.listdir() if os.path.isdir(e) and e[0] != '.']

tabs = st.tabs(sub_apps)
_dir = os.getcwd()

for name, tab in zip(sub_apps, tabs):
    with tab:
        try:
            import importlib
            import os
            sys.path.append(name)
            spec = importlib.util.spec_from_file_location(
                name, os.path.join(name, name+'.py'))
            foo = importlib.util.module_from_spec(spec)
            sys.modules[name] = foo
            spec.loader.exec_module(foo)
        except Exception as err:
            st.error(err)
