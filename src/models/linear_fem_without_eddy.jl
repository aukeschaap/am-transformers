
using gmsh

using LinearAlgebra, SparseArrays
using WriteVTK

using BenchmarkTools

include("../constants.jl")
include("../get_mesh_data.jl")
include("../process.jl")
include("../definitions/fem2d.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


############################################################

# Linear FEM without eddy currents

############################################################

gmsh.open(MESH_LOCATION)
mshdata = get_mesh_data();


# Constants


Ip = 0;       # Primary peak phase current
Is = 777.62;  # Secondary peak phase current
Np = 266;
Ns = 6;

omega = 2*pi*50;  # Frequency

# Calculate current density in the windings
Jp = Np * Ip / Awhv;
Js = Ns * Is / Awlv;



# Source current density J
#  One term for each of the windings, with a positive and negative part
#  Note the phase shift between the phases
sourcefunction(group_id) = Jp * exp(1im * 2pi/3) * (-1 * (group_id==3) + 1 * (group_id==4)) + 
                           Jp * (-1 * (group_id==5) + 1 * (group_id==6)) + 
                           Jp * exp(-1im * 2pi/3) * (-1 * (group_id==7) + 1 * (group_id==8)) + 
                           Js * exp(1im * 2pi/3) * (1 * (group_id==9) - 1 * (group_id==10)) +
                           Js * (1 * (group_id==11) - 1 * (group_id==12)) + 
                           Js * exp(-1im * 2pi/3) * (1 * (group_id==13) - 1 * (group_id==14));
sourceperelement = map(sourcefunction, mshdata.e_group);

# Relative permeability model
mu0 = 4e-7 * pi;
mur = 1000;       # Relative permeability of the core
reluctivityfunction(group_id) = (1 / mu0) + (1/(mu0*mur) - 1/mu0) * (group_id == 2)
reluctivityperelement = map(reluctivityfunction, mshdata.e_group);

# Conductivity
conductivityfunction(group_id) = 0;
conductivityperelement = map(conductivityfunction, mshdata.e_group);


# Calculate the vector potential
u = fem2d(mshdata, sourceperelement, reluctivityperelement, conductivityperelement, omega);

# Post-process for magnetic field and current density
B, H, Wm, Jel = process(mshdata, u, sourceperelement, reluctivityperelement, conductivityperelement, omega);
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
outfiles = vtk_save(vtkfile);

