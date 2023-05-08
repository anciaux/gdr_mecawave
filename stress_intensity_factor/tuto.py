#!/usr/bin/env python
# coding: utf-8

import subprocess
import akantu as aka
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from scipy.optimize import curve_fit
# setting a default image size large enough
plt.rcParams['figure.figsize'] = [5, 5]


def createMesh(l=1, L=10, h1=0.05, h2=.5, quiet=False):

    geometry_file = f"""
h1 = {h1};
h2 = {h2};
L = {L};
l = {l};
Point(1) = {{0, 0, 0, h2}};
Point(2) = {{2*L, 0, 0, h2}};
Point(3) = {{2*L, L, 0, h2}};
Point(4) = {{0, L, 0, h2}};
Point(5) = {{l, 0, 0, h1}};

Point(6) =  {{0, 0, 0, h2}};
Point(7) =  {{2*L, -L, 0, h2}};
Point(8) =  {{0, -L, 0, h2}};


Line(1) = {{1, 5}};
Line(2) = {{4, 1}};
Line(3) = {{3, 4}};
Line(4) = {{2, 3}};
Line(5) = {{5, 2}};

Line Loop(1) = {{2, 3, 4, 5, 1}};
Plane Surface(1) = {{1}};

Line(6) =  {{5, 6}};
Line(7) =  {{6, 8}};
Line(8) =  {{8, 7}};
Line(9) =  {{7, 2}};
Line Loop(2) = {{6, 7, 8, 9, -5}};
Plane Surface(2) = {{2}};


Physical Surface(8) = {{1,2}};
Physical Line("left") = {{2,7}};
Physical Line("bottom") = {{8}};
Physical Line("top") = {{3}};
Physical Line("right") = {{4,9}};

"""

    with open('plate.geo', 'w') as f:
        f.write(geometry_file)

    p = subprocess.Popen(
        "gmsh -2 -order 1 -o plate.msh plate.geo", shell=True,
        stdout=subprocess.PIPE)
    p.wait()
    if p.returncode:
        print("Beware, gmsh could not run: mesh is not regenerated")
    elif quiet is False:
        print("Mesh generated")

    spatial_dimension = 2
    mesh = aka.Mesh(spatial_dimension)
    mesh.read('plate.msh')
    return mesh

################################################################


def createModel(l=1, L=10, h1=0.05, h2=.5, U=0.005, **kwargs):

    mesh = createMesh(l=l, L=L, h1=h1, h2=h2, **kwargs)

    material_file = """
material elastic [
    name = virtual
    rho = 1    # density
    E   = 1    # young's modulus
    nu  = 0.3  # poisson's ratio
    finite_deformation = false
]
"""

    with open('material.dat', 'w') as f:
        f.write(material_file)

    aka.parseInput('material.dat')

    model = aka.SolidMechanicsModel(mesh)
    model.initFull(_analysis_method=aka._static)

    solver = model.getNonLinearSolver('static')
    solver.set('max_iterations', 100)
    solver.set('threshold', 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)

    model.applyBC(aka.FixedValue(U, aka._y), 'top')
    model.applyBC(aka.FixedValue(-U, aka._y), 'bottom')

    model.getExternalForce()[:] = 0

    model.setBaseName('plate')
    model.addDumpFieldVector('displacement')
    model.addDumpFieldVector('external_force')
    model.addDumpField('strain')
    model.addDumpField('stress')
    model.addDumpField('blocked_dofs')

    model.solveStep()
    model.dump()
    return model, mesh


def plotMesh(mesh, displacement=None, field=None):
    conn = mesh.getConnectivity(aka._triangle_3)
    nodes = mesh.getNodes()
    if displacement is None:
        triangles = tri.Triangulation(nodes[:, 0], nodes[:, 1], conn)
    else:
        triangles = tri.Triangulation(nodes[:, 0]+displacement[:, 0],
                                      nodes[:, 1]+displacement[:, 1], conn)

    fig = plt.Figure()
    axe = fig.add_subplot(111)
    axe.set_aspect('equal')
    axe.triplot(triangles, '--', lw=.8)
    return fig


def plotResult(model, displacement=None, field=None,
               log_color=False, contour=None,
               contour_label=True,
               xrange=None, yrange=None):
    mesh = model.getMesh()
    conn = mesh.getConnectivity(aka._triangle_3)
    nodes = mesh.getNodes()

    fig = plt.Figure()
    axe = fig.add_subplot(111)
    axe.set_aspect('equal')

    pos = nodes.copy()

    if displacement is not None:
        pos += displacement

    triangles = tri.Triangulation(pos[:, 0], pos[:, 1], conn)
    contour_fmt = None
    if field == 'stress':
        field = model.getMaterial(0).getStress(aka._triangle_3)
        nodal_field = np.zeros((nodes.shape[0], field.shape[1]))
        for el, c in enumerate(conn):
            for n in c:
                nodal_field[n, :] += field[el]/3
        label = 'stress magnitude [Pa]'
        def stress_fmt(x): return "$" + str(round(x, 2)) + "\\quad [Pa]$"
        contour_fmt = stress_fmt

    elif field == 'displacement':
        field = model.getDisplacement()
        label = 'displacement [m]'
        def disp_fmt(x): return "$" + str(round(x, 2)) + "\\quad [m]$"
        contour_fmt = disp_fmt
        nodal_field = field

    color = np.linalg.norm(field, axis=1)
    if log_color:
        color = np.log(color)
    display = axe.tripcolor(triangles, color)

    if xrange is not None:
        axe.set_xlim(xrange)
    if yrange is not None:
        axe.set_ylim(yrange)

    if contour is not None:
        f = np.linalg.norm(nodal_field, axis=1)
        CS = axe.tricontour(triangles, f, contour,
                            linewidths=2, colors='red')
        if contour_label:
            axe.clabel(CS, CS.levels, fmt=contour_fmt, inline=True)
    else:
        cbar = fig.colorbar(display)
        cbar.set_label(label)

    return fig


def extract_stress(model, theta=0, r_fit=1, angle_threshold=1e-1,
                   l=None, shear=False, **params):
    mesh = model.getMesh()
    quad_coords = aka.ElementTypeMapArrayReal()
    quad_coords.initialize(mesh, nb_component=2)

    fe_engine = model.getFEEngine()
    fe_engine.interpolateOnIntegrationPoints(
        mesh.getNodes(), quad_coords)

    q_points = quad_coords(aka._triangle_3)
    tip = np.array((l, 0))
    curve = []
    stress_field = model.getMaterial(0).getStress(aka._triangle_3)

    for coord, stress in zip(q_points, stress_field):
        stress = stress.reshape((2, 2))
        p = coord-tip
        r = np.linalg.norm(p)
        _theta = np.arctan2(p[1], p[0])
        if np.abs(theta - _theta) > angle_threshold:
            continue
        n = p.T/r
        if shear:
            n = np.array([[0, 1], [-1, 0]])
        sigma = stress@n

        curve.append((r, np.linalg.norm(sigma)))

    curve = np.array(curve)
    curve = curve[curve[:, 0].argsort()]

    def f_williams_I(r, K, theta):
        sigma = K/np.sqrt(2*np.pi*r)*(5/4*np.cos(theta)-1/4*np.cos(3*theta/2))
        return sigma

    fit_range = curve[:, 0] < r_fit
    res = curve_fit(lambda r, K: f_williams_I(
        r, K, theta),
        curve[fit_range, 0],
        curve[fit_range, 1])
    K = res[0][0]
    return curve, K


def full_study_extract_K(**params):
    model, mesh = createModel(**params)
    plotMesh(mesh, displacement=model.getDisplacement())
    curve, K = extract_stress(model, theta=0, r_fit=.2, shear=True, **params)
    epot = model.getEnergy('potential')
    return K, epot
