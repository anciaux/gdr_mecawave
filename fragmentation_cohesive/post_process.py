import numpy as np
import akantu as aka

import generate_model as DFModel


def getDamageParameter(model):

    d = model.getMaterial(1).getInternalReal("damage")
    d = d(aka._cohesive_2d_4)
    # mat = model.getMaterial(1)
    # coh_id = model.getMaterial("insertion").getElementFilter()(
    #  aka._cohesive_2d_4
    # )
    # model.getFEEngine().computeIntegrationPointsCoordinate(quad_coords, model.getMaterial("cohesive material").getElementFilter())
    return d


def getEnergy(energies, energy_name):
    for i in range(len(energies)):
        if energies[i][0] == energy_name:
            energy = energies[i][1]
            return energy


def computeVariationEnergy(energies):
    """Returns the variation of energies between the current time step and the time step 0."""

    var_energy_potential = np.zeros(DFModel.n_steps)
    var_energy_kinetic = np.zeros(DFModel.n_steps)
    var_energy_dissipated = np.zeros(DFModel.n_steps)
    var_energy_contact = np.zeros(DFModel.n_steps)
    var_energy_reversible = np.zeros(DFModel.n_steps)
    var_energy_total = np.zeros(DFModel.n_steps)
    var_external_work = np.zeros(DFModel.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFModel.n_steps):
        var_energy_potential[n] = energy_potential[n] - energy_potential[0]
        var_energy_kinetic[n] = energy_kinetic[n] - energy_kinetic[0]
        var_energy_dissipated[n] = energy_dissipated[n] - energy_dissipated[0]
        var_energy_reversible[n] = energy_reversible[n] - energy_reversible[0]
        var_energy_contact[n] = energy_contact[n] - energy_contact[0]
        var_external_work[n] = external_work[n] - external_work[0]
        var_energy_total[n] = var_external_work[n] - (
            var_energy_potential[n]
            + var_energy_kinetic[n]
            + var_energy_dissipated[n]
            + var_energy_reversible[n]
            + var_energy_contact[n]
        )

    var_energies = [
        ["var energy potential", var_energy_potential],
        ["var energy kinetic", var_energy_kinetic],
        ["var energy dissipated", var_energy_dissipated],
        ["var energy reversible", var_energy_reversible],
        ["var energy contact", var_energy_contact],
        ["var external work", var_external_work],
        ["var energy total", var_energy_total],
    ]

    return var_energies


def computePower(energies):
    """Returns the variation of energies between two consecutives time steps."""

    power_potential = np.zeros(DFModel.n_steps)
    power_kinetic = np.zeros(DFModel.n_steps)
    power_dissipated = np.zeros(DFModel.n_steps)
    power_reversible = np.zeros(DFModel.n_steps)
    power_contact = np.zeros(DFModel.n_steps)
    power_external_work = np.zeros(DFModel.n_steps)
    power_total = np.zeros(DFModel.n_steps)

    energy_potential = getEnergy(energies, "energy potential")
    energy_kinetic = getEnergy(energies, "energy kinetic")
    energy_dissipated = getEnergy(energies, "energy dissipated")
    energy_reversible = getEnergy(energies, "energy reversible")
    energy_contact = getEnergy(energies, "energy contact")
    external_work = getEnergy(energies, "external work")

    for n in range(1, DFModel.n_steps):
        power_potential[n] = energy_potential[n] - energy_potential[n - 1]
        power_kinetic[n] = energy_kinetic[n] - energy_kinetic[n - 1]
        power_dissipated[n] = energy_dissipated[n] - energy_dissipated[n - 1]
        power_reversible[n] = energy_reversible[n] - energy_reversible[n - 1]
        power_contact[n] = energy_contact[n] - energy_contact[n - 1]
        power_external_work[n] = external_work[n] - external_work[n - 1]
        power_total[n] = power_external_work[n] - (
            power_potential[n]
            + power_kinetic[n]
            + power_dissipated[n]
            + power_reversible[n]
            + power_contact[n]
        )

    power = [
        ["power potential", power_potential],
        ["power kinetic", power_kinetic],
        ["power dissipated", power_dissipated],
        ["power reversible", power_reversible],
        ["power contact", power_contact],
        ["power external work", power_external_work],
        ["power total", power_total],
    ]

    return power
