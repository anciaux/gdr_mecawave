import subprocess
import akantu as aka
import argparse
# import input_files.input_data_aka as inputdata


def generate_mesh(**params):
    inputdata = argparse.Namespace(**params)

    # Assign geometry
    bar_length = inputdata.bar_length
    x0 = inputdata.x0
    xf = inputdata.xf

    n_elements = inputdata.number_elements
    # Mesh
    # Compute the size of elements (h) for uniform mesh
    h_uniform = bar_length / (n_elements * 0.5)

    geometry_file = f"""
    Point(1) = {{ {x0}, 0, 0, {h_uniform} }};
    Point(2) = {{ {xf}, 0, 0, {h_uniform} }};
    Point(3) = {{ {xf}, {h_uniform}, 0, {h_uniform} }};
    Point(4) = {{ {x0}, {h_uniform}, 0, {h_uniform} }};

    Line(1) = {{1,2}};
    Line(2) = {{2,3}};
    Line(3) = {{3,4}};
    Line(4) = {{4,1}};

    Line Loop(1) = {{1,2,3,4}};
    Plane Surface(1) = {{1}};
    Physical Surface(1) = {{1}};
    Physical Line("Yblocked") = {{1,3}};
    Physical Line("right") = {{2}};
    Physical Line("left") = {{4}};

    Transfinite Surface {1} Left;
    """

    with open("bar.geo", "w") as f:
        f.write(geometry_file)

    p = subprocess.Popen(
        "gmsh -2 -order 1 -o bar.msh bar.geo", shell=True,
        stdout=subprocess.PIPE
    )
    p.wait()
    if p.returncode:
        print("FATAL    : Gmsh error")
    else:
        print("Info    : Mesh generated")

    # Read mesh
    spatial_dimension = 2
    mesh = aka.Mesh(spatial_dimension)
    mesh.read("bar.msh")
    return mesh
