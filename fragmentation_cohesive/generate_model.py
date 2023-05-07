import akantu as aka
import numpy as np
import generate_mesh
import argparse
################################################################
# Apply BCs


class FixedVelocity(aka.DirichletFunctor):
    """Fixed velocity at the boundaries."""

    def __init__(self, axis, vel):
        super().__init__(axis)
        self.axis = axis
        self.time = 0
        self.vel = vel

    def set_time(self, t):
        self.time = t

    def __call__(self, node, flags, disp, coord):
        flags[int(self.axis)] = True
        disp[int(self.axis)] = self.vel * self.time
################################################################


def create_model(**params):
    spatial_dimension = 2
    mesh = generate_mesh.generate_mesh(**params)
    node_coord = mesh.getNodes()
    inputdata = argparse.Namespace(**params)
    # Assign material properties
    young_modulus = inputdata.young_modulus
    rho = inputdata.density
    fracture_energy = inputdata.fracture_energy
    stress_limit = inputdata.stress_limit

    # Assign load
    strain_rate = inputdata.strain_rate

    # Contact penalty

    alpha = inputdata.contact * (
        stress_limit**2
        + 4.5
        * strain_rate ** (2 / 3)
        * young_modulus
        * fracture_energy ** (2 / 3)
        * rho ** (1 / 3)
    ) / (4.5 * fracture_energy)

    if inputdata.variation_sigma_c > 0:
        var = inputdata.variation_sigma_c
        sigma_c = f"sigma_c = {stress_limit} uniform [-{var}, {var}] # critical stress (Pa)"
    else:
        sigma_c = f"sigma_c = {stress_limit}"
    material_file = f"""
        seed = 1.0
        model solid_mechanics_model_cohesive [

            material elastic [
                name = linear
                rho = {rho} # Density (kg/m3)
                E = {young_modulus}  # Young's module (Pa)
                nu = 0.
                finite_deformation = true
            ]

            material cohesive_linear [
                name = insertion\n"""
    material_file += f"""
                {sigma_c}
    """
    material_file += f"""
                G_c = {fracture_energy} # Fracture energy (N/m)
                beta = 0.
                penalty = {alpha}
            ]
        ]
        """

    with open("material.dat", "w") as f:
        f.write(material_file)

    # Read material file to akantu
    aka.parseInput("material.dat")

    # Create model a model
    model = aka.SolidMechanicsModelCohesive(mesh)
    model.initFull(_analysis_method=aka._static, _is_extrinsic=True)
    # Configure static solver
    solver = model.getNonLinearSolver("static")
    solver.set("max_iterations", 100)
    solver.set("threshold", 1e-10)
    solver.set("convergence_type", aka.SolveConvergenceCriteria.residual)
    # Add solver
    model.initNewSolver(aka._explicit_lumped_mass)

    # Configure cohesive element insertion
    use_automatic_insertion = False
    if use_automatic_insertion is True:
        model.updateAutomaticInsertion()
    else:
        mesh_facets = mesh.getMeshFacets()
        connect_facets = mesh_facets.getConnectivities()
        # facets_coords = mesh_facets.getNodes()

        # Get the cohesive inserter and check facets
        cohesive_inserter = model.getElementInserter()
        check_facets = cohesive_inserter.getCheckFacets()
        up = np.array([0.0, 1.0])
        # Loop by all facet types used in the simulation
        for facet_type in connect_facets.elementTypes(
                dim=(spatial_dimension - 1)):
            conn_facettype = connect_facets(facet_type)
            check_facettype = check_facets(facet_type)
            for el, conn_facettype in enumerate(conn_facettype):
                # Check the direction of the vector
                dir_vec = (
                    node_coord[conn_facettype[1], :]
                    - node_coord[conn_facettype[0], :]
                )
                direction = (dir_vec / np.linalg.norm(dir_vec)).dot(up)
                # If the direction is not 1 it means that is a diagonal facet,
                # then assign False
                if abs(direction) < 0.99:
                    check_facettype[el] = False

    # Apply Dirichlet BC to block dispacements at y direction on
    # top and botton of the elements
    model.applyBC(aka.FixedValue(0.0, aka._y), "Yblocked")

    # Assign load
    strain_rate = inputdata.strain_rate

    model.applyBC(aka.FixedValue(0, aka._x), "left")
    model.applyBC(aka.FixedValue(0, aka._x), "right")

    # Initial values
    v0 = model.getVelocity()
    v0[:, 0] = np.array([inputdata.strain_rate * x for x,
                        y in mesh.getNodes()])
    # model.getVelocity()[:] = v0
    return model, mesh


def applyVel(model, applied_vel, current_time_step, dt):
    time = current_time_step * dt
    model.applyBC(aka.FixedValue(-applied_vel*time, aka._x), "left")
    model.applyBC(aka.FixedValue(+applied_vel*time, aka._x), "right")
