import akantu as aka
import numpy as np
import time
import generate_model
import plot
from app_helper import tqdm, make_figure
import os
import argparse
import pandas as pd
import matplotlib.pyplot as plt
import streamlit as st
# Input model for src_akantu
################################################################
params = {
    'uniform_mesh': True,
    # Material
    'young_modulus':  275.0 * 1e9,  # (Pa)
    'density':  2750.0,  # (kg/m3)
    'fracture_energy':  100.0,  # (N/m)
    'stress_limit':  300.0 * 1e6,  # (Pa)
    # Geometry
    'bar_length':  50 * 1e-3,  # (m)
}
params.update({
    # Left extremitiy x coordinate / 0-initial
    'x0': -0.5 * params['bar_length'],
    # Rigth extremitiy x coordinate / f-final
    'xf':  0.5 * params['bar_length'],
    'number_elements':  50 * 2,  # Total number of triangular elements
    'area':  1.0,  # Cross sectional area (m2) (Equal to element size )
    # Load
    'strain_rate':  1e4,  # (s-1)
    # Time
    'time_simulation':  4.0 * 1e-7,  # Total time of simulation (s)

    # if there is previous data to continue the simulation
    # from a previous simulation set
    # continue_simulation_from_step = True and give the time
    # to start the simulation
    'initial_step':  0,
    'continue_simulation_from_step':  False,

    'half_bar':  False,
    # if use symmetry we have to add the bc proper
    'dump_freq': 10,
    'paraview': False
})


################################################################


def runSimulation(model, **params):

    inputdata = argparse.Namespace(**params)

    L = inputdata.bar_length
    h = L/inputdata.number_elements

    # Initiation of variables
    work_previous_step = 0.0

    if inputdata.paraview:
        plot.addPlotVtk(model)
    # Set time increment
    dt_crit = model.getStableTimeStep()
    dt = dt_crit * 0.1
    model.setTimeStep(dt)

    n_steps = int(inputdata.time_simulation / dt)
    applied_vel = inputdata.strain_rate * inputdata.bar_length * 0.5

    data_bc = [0, 0]

    results = []
    results_energies = []

    for n in tqdm(range(n_steps), init_text="Integrating time step"):

        # Apply velocity at the boundaries
        generate_model.applyVel(model, applied_vel, n+1, dt)

        if n % inputdata.dump_freq == 0:
            if inputdata.paraview:
                plot.addVtkFiles(model, n)

        # Run simulation
        model.checkCohesiveStress()
        model.solveStep("explicit_lumped")

        u = model.getDisplacement()[:, 0]
        v = model.getVelocity()[:, 0]
        acel = model.getAcceleration()[:, 0]

        fint = model.getInternalForce()[:, 0]
        stress = model.getMaterial(0).getStress(aka._triangle_3)
        stress_xx = stress[:, 0]
        avg_stress_bar = np.mean(stress_xx)

        # Energy balance
        energies = computeEnergies(
            model, work_previous_step, fint, data_bc, applied_vel, dt)

        work_previous_step = energies["external work"]
        results_energies.append([n*dt] + [v for k, v in energies.items()])
        headers_energies = ['time'] + [k for k, v in energies.items()]
        d = getDamageParameter(model)

        # Fragmentation data
        fragment_data = aka.FragmentManager(model)
        fragment_data.computeAllData()

        # Number of fragments
        n_fragments = fragment_data.getNbFragment()
        # Fragments size (assuming uniform mesh)
        frag_lengths = fragment_data.getNbElementsPerFragment() * h / L

        result_step = [
            n*dt,
            n,
            u,
            v,
            acel,
            d,
            stress,
            avg_stress_bar,
            n_fragments,
            frag_lengths
        ]
        results.append(result_step)

    results = pd.DataFrame(
        results,
        columns=["time",
                 "step",
                 "displacement",
                 "velocity",
                 "acceleration",
                 "damage",
                 "stress",
                 "avg_stress_bar",
                 "n_fragments",
                 "frag_lengths"])

    results_energies = pd.DataFrame(
        results_energies, columns=headers_energies)

    st.markdown(
        f"### Computed {dt} * {n_steps} = {inputdata.time_simulation} seconds")
    return results, results_energies
################################################################


def computeEnergies(model, work_previous_step, fint_current_step,
                    data_bc, applied_vel, dt):

    energy_potential = model.getEnergy("potential")
    energy_kinetic = model.getEnergy("kinetic")
    energy_dissipated = model.getEnergy("dissipated")
    energy_reversible = model.getEnergy("reversible")
    energy_contact = model.getEnergy("cohesive contact")

    mesh = model.getMesh()
    # External work
    nodes_bc_left = mesh.getElementGroup("left").getNodeGroup().getNodes()
    nodes_bc_right = mesh.getElementGroup("right").getNodeGroup().getNodes()
    # Internal for at the current time step
    fint_current_step_bcleft = -np.sum(fint_current_step[nodes_bc_left])
    fint_current_step_bcright = -np.sum(fint_current_step[nodes_bc_right])
    # The reaction force is taken as an average between the internal
    # force in the current and previous time-set
    freact_previous_step_bcleft = data_bc[0]
    freact_previous_step_bcright = data_bc[1]

    freact_bcleft = (fint_current_step_bcleft +
                     freact_previous_step_bcleft) * 0.5
    freact_bcright = (fint_current_step_bcright +
                      freact_previous_step_bcright) * 0.5

    external_work = (
        work_previous_step
        + (freact_bcleft * -applied_vel +
           freact_bcright * applied_vel)
        * dt
    )

    data_bc[:] = [freact_bcleft, freact_bcright]

    energies_step = {
        "energy potential": energy_potential,
        "energy kinetic": energy_kinetic,
        "energy dissipated": energy_dissipated,
        "energy reversible": energy_reversible,
        "energy contact": energy_contact,
        "external work": external_work,
    }

    return energies_step

################################################################


def computeVarEnergiesCZM(energies, n_steps, n_elements):
    """Returns the variation of energies between the current time step
       and the time step 0.
    """

    h = 50*1e-3 / n_elements

    var_energies = (energies - energies.iloc[0]) / h
    var_energies.columns = ['var ' + v for v in energies.columns]
    var_energies['var energy total'] = var_energies['var external work'] - (
        var_energies['var energy potential']
        + var_energies['var energy kinetic']
        + var_energies['var energy dissipated']
        + var_energies['var energy reversible']
        + var_energies['var energy contact']
    )
    var_energies['var time'] = energies['time']

    return var_energies
################################################################


def plotResults(results, label_x, label_y):
    with make_figure() as (fig, axe):
        axe.grid(True, which="both")
        axe.axhline(y=0, color="k")
        axe.set_xlabel(label_x)
        axe.set_ylabel(label_y)

        x_values = results[results.columns[0]]
        headers = [c for c in results.columns]
        for label in headers[1:]:
            y_values = results[label]
            axe.plot(x_values, y_values, label=label)
        axe.legend()


################################################################


def plotVarEnergiesCZM(var_energies, title):
    """Plot variation of energy from time t to t0."""

    var_energy_potential = var_energies["var energy potential"] / 1e3
    var_energy_kinetic = var_energies["var energy kinetic"] / 1e3
    var_energy_dissipated = var_energies["var energy dissipated"] / 1e3
    var_energy_reversible = var_energies["var energy reversible"] / 1e3
    var_energy_contact = var_energies["var energy contact"] / 1e3
    var_external_work = var_energies["var external work"] / 1e3
    var_energy_total = var_energies["var energy total"] / 1e3

    with make_figure() as (fig, axe):
        axe.grid(True, which="both")
        axe.axhline(y=0, color="k")
        axe.set_title(title)
        axe.set_xlabel(str("Time (s)"))
        axe.set_ylabel("Variation of energy ($ kJ/ {m^2} $)")

        x = var_energies['var time']
        axe.plot(x, var_energy_potential, label="varEpot")
        axe.plot(x, var_energy_kinetic, label="varEkin")
        axe.plot(x, var_energy_dissipated, label="varEdis")
        axe.plot(x, var_energy_reversible, label="varErev")
        axe.plot(x, var_energy_contact, label="varEcon")
        axe.plot(x, -var_external_work, label="-varWext")
        axe.plot(x, var_energy_total, label="varEtot")
        axe.legend()

################################################################


def plotFragmentSizeHistogram(fragment_sizes, bins=10, **kwargs):
    with make_figure() as (fig, axe):
        axe.grid(True, which="both")
        axe.axhline(y=0, color="k")
        axe.set_title("Fragment size distribution")

        counts, bins = np.histogram(fragment_sizes, bins=bins)
        bin_centers = 0.5*(bins[:-1] + bins[1:])
        non_zero = counts > 0
        counts = counts[non_zero]
        bins = bin_centers[non_zero]
        axe.plot(bins, counts, 'o-')
        axe.set_xlabel('Fragment size s/L')
        axe.set_ylabel("Number of fragments")
        axe.set_xlim(xmin=0, xmax=bins.max())

################################################################


def analyze_results(model, results, energies, **params):

    inputdata = argparse.Namespace(**params)
    dt_crit = model.getStableTimeStep()
    dt = dt_crit * 0.1

    # time_data = readResultsVariable(file_address, "time_data")

    time_simulation = inputdata.time_simulation
    n_steps = int(time_simulation / dt)
    n_files = int(n_steps / 10 + 1)

    plotResults(
        results[['time', 'avg_stress_bar']],
        label_x="time (s)",
        label_y="Average stress at the bar (MPa)",
    )
    plotResults(
        results[['time', 'n_fragments']],
        label_x="time (s)",
        label_y="N",
    )

    var_energies = computeVarEnergiesCZM(energies, n_files, 1250)
    plotVarEnergiesCZM(var_energies, title="CZM: Variation of energy")

    frag_lengths = results['frag_lengths'].iloc[-1]
    plotFragmentSizeHistogram(frag_lengths, **params)

################################################################


def getDamageParameter(model):

    d = model.getMaterial(1).getInternalReal("damage")
    d = d(aka._cohesive_2d_4)
    # mat = model.getMaterial(1)
    # coh_id = model.getMaterial("insertion").getElementFilter()(
    #  aka._cohesive_2d_4
    # )
    # model.getFEEngine().computeIntegrationPointsCoordinate(quad_coords,
    # model.getMaterial("cohesive material").getElementFilter())
    return d

################################################################


def main():
    try:
        os.mkdir("output")
    except FileExistsError:
        pass

    start_time = time.time()
    model, mesh = generate_model.create_model(**params)
    results, energies = runSimulation(model, **params)
    total_time = time.time() - start_time
    analyze_results(model, results, energies, **params)
    plt.show()

    print("--- %s seconds ---" % (total_time))


################################################################
if __name__ == "__main__":
    main()
