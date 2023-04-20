
using gmsh

using LinearAlgebra, SparseArrays
using WriteVTK

using BenchmarkTools

include("../constants.jl")
include("../get_mesh_data.jl")
include("../process.jl")
include("../definitions/assemble_Kf.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


############################################################
#
#       Steady state solution
#
# This is the solution using omega = 0.
#
############################################################

gmsh.open(MESH_LOCATION)
mshdata = get_mesh_data();
println("Successfully loaded mesh data.")

# Constants

"Primary peak phase current"
Ip = 0;

"Secondary peak phase current"
Is = 777.62;


Np = 266;
Ns = 6;

omega = 0;  # Frequency

# Calculate current density in the windings
Jp = Np * Ip / Awhv;
Js = Ns * Is / Awlv;


"""
# Source current density J

One term for each of the windings, with a positive and negative part. Note the phase shift between
the phases.
"""
function source_func(group_id)
    Jp * exp(1im * 2pi/3) * (-1 * (group_id==3) + 1 * (group_id==4))
    + Jp * (-1 * (group_id==5) + 1 * (group_id==6))
    + Jp * exp(-1im * 2pi/3) * (-1 * (group_id==7) + 1 * (group_id==8))
    + Js * exp(1im * 2pi/3) * (1 * (group_id==9) - 1 * (group_id==10))
    + Js * (1 * (group_id==11) - 1 * (group_id==12))
    + Js * exp(-1im * 2pi/3) * (1 * (group_id==13) - 1 * (group_id==14));
end
source_per_element = map(source_func, mshdata.e_group);

# Relative permeability model
mu0 = 4e-7 * pi;
mur = 1000;       # Relative permeability of the core
reluctivity_func(group_id) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
reluctivity_per_element = map(reluctivity_func, mshdata.e_group);

# Conductivity
conductivity_func(group_id) = 0;
conductivity_per_element = map(conductivity_func, mshdata.e_group);


println("Assembling steady state...")

# Calculate the vector potential
K, f = assemble_steadystate(mshdata, source_per_element, reluctivity_per_element);

println("System complete. Solving...")
u = K \ f;
println("Done solving.")

# Post-process for magnetic field and current density
B, H, Wm, Jel = process(mshdata, u, source_per_element, reluctivity_per_element, conductivity_per_element, omega);
Bnorm = norm.(sqrt.(B[1].^2 + B[2].^2));


# Define nodes (points) and elements (cells)
points = [mshdata.xnode mshdata.ynode]';
cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mshdata.elements];

# Create VTK file structure using nodes and elements
vtkfile = vtk_grid(string(OUTPUT_LOCATION, "transformer1.vtu"), points, cells);

# Store data in the VTK file
vtkfile["Az", VTKPointData()]   = norm.(u);
vtkfile["imA", VTKPointData()]  = imag.(u);
vtkfile["Bnorm", VTKCellData()] = Bnorm;
vtkfile["Jel", VTKCellData()]   = Jel;

# Save the file
println("Saving result in a file.")
outfiles = vtk_save(vtkfile);

println("Model complete.")
