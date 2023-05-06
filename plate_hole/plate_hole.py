#!/usr/bin/env python
# coding: utf-8

import subprocess
import matplotlib.tri as tri
import matplotlib.pyplot as plt
import numpy as np
import akantu as aka
import solidipes as sp
import os
import streamlit as st

dirname = os.path.dirname(__file__)

params = {
    'w': 10.,     # width (x-axis)
    'l': 5.,      # length (y-axis)
    'h1': 0.05,   # characteristic mesh size at the hole
    'h2': 0.3,    # characteristic mesh size in corners
    'R': 2.      # radius of the hole
}

col1, col2 = st.columns(2)

with col1:
    f = sp.load_file(os.path.join(dirname, 'plate-hole-2.svg'))
    f.view()

for p, value in params.items():
    _type = type(value)
    val = col2.text_input(p, value=value, key=p+dirname)
    params[p] = _type(val)


material_file = """
material elastic [
    name = steel
    rho = 1        # density
    E   = 1        # young's modulus
    nu  = 1        # poisson's ratio
]"""
# writing the material file
open('material.dat', 'w').write(material_file)
# reading the material file
material_file = 'material.dat'
aka.parseInput(material_file)

# geometric parameters
mesh_file = f"""
Point(1) = {{0, 0, 0, {h2} }};
Point(2) = {{ {R}, 0, 0, {h1} }};
Point(3) = {{ {w}, 0, 0, {h2} }};
Point(4) = {{ {w}, {l}, 0, {h2} }};
Point(5) = {{ 0,   {l}, 0, {h2} }};
Point(6) = {{0,    {R}, 0, {h1} }};
"""

mesh_file += """
Circle(1) = {6, 1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line Loop(6) = {1, 2, 3, 4, 5};
Plane Surface(7) = {6};
"""

mesh_file += """
Physical Surface(8) = {7};
Physical Line("XBlocked") = {5};
Physical Line("YBlocked") = {2};
Physical Line("Traction") = {3};
"""


open('plate.geo', 'w').write(mesh_file)

ret = subprocess.run("gmsh -2 -order 1 -o plate.msh plate.geo", shell=True)
if ret.returncode:
    print("Beware, gmsh could not run: mesh is not regenerated")
else:
    print("Mesh generated")


# reading the mesh
spatial_dimension = 2
mesh_file = 'plate.msh'
mesh = aka.Mesh(spatial_dimension)
mesh.read(mesh_file)


# extract the mesh
conn = mesh.getConnectivity(aka._triangle_3)
nodes = mesh.getNodes()
triangles = tri.Triangulation(nodes[:, 0], nodes[:, 1], conn)

# plot the result
plt.axes().set_aspect('equal')
# plots the pristine state
t = plt.triplot(triangles, '--', lw=.8)


# creating the solid mechanics model
model = aka.SolidMechanicsModel(mesh)


# initialize a static solver
model.initFull(_analysis_method=aka._static)


# set the displacement/Dirichlet boundary conditions
model.applyBC(aka.FixedValue(0.0, aka._x), "XBlocked")
model.applyBC(aka.FixedValue(0.0, aka._y), "YBlocked")


# set the force/Neumann boundary conditions
model.getExternalForce()[:] = 0

trac = [10, 0]  # Newtons/m^2

model.applyBC(aka.FromTraction(trac), "Traction")


# configure the linear algebra solver
solver = model.getNonLinearSolver()
solver.set("max_iterations", 2)
solver.set("threshold", 1e-10)
solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

# compute the solution
model.solveStep()

# extract the displacements
u = model.getDisplacement()

# plot the result
plt.axes().set_aspect('equal')
# plots the pristine state
t = plt.triplot(triangles, '--', lw=.8)
# plots an exagerating view of the strained mesh
t = plt.triplot(nodes[:, 0]+u[:, 0]*2e9, nodes[:, 1] +
                u[:, 1]*2e9, triangles=conn)

# plot displacement field
plt.axes().set_aspect('equal')

u_disp = plt.tricontourf(triangles, np.linalg.norm(u, axis=1))
t = plt.triplot(triangles, '--', lw=.8)

cbar = plt.colorbar(u_disp)
cbar.set_label('displacement magnitude [m]')

# plot stress field
plt.axes().set_aspect('equal')
stress_field = model.getMaterial(0).getStress(aka._triangle_3)
stress_disp = plt.tripcolor(triangles, np.linalg.norm(stress_field, axis=1))
cbar = plt.colorbar(stress_disp)
cbar.set_label('stress magnitude [Pa]')


# specify what field to output into paraview files
model.setBaseName("plate")
model.addDumpFieldVector("displacement")
model.addDumpFieldVector("external_force")
model.addDumpField("strain")
model.addDumpField("stress")
model.addDumpField("blocked_dofs")

# generate paraview files
model.dump()
