
# Linux uses Gmsh, Windows uses gmsh
try
    using Gmsh: gmsh
catch
    using gmsh
end

using LinearAlgebra, SparseArrays
using Logging

using WriteVTK
using BenchmarkTools

include("../constants.jl")
include("../get_mesh_data.jl")
include("../process.jl")
include("../definitions/source.jl")
include("../definitions/linear_reluctivity.jl")
include("../definitions/assemble_Kf.jl")

const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


"Primary peak phase current"
Ip = 0;

"Secondary peak phase current"
Is = 777.62;

Np = 266;
Ns = 6;

# Calculate current density in the windings
Jp = Np * Ip / Awhv;
Js = Ns * Is / Awlv;

# Vacuum permeability
μ_0 = 4e-7 * pi;

# Relative permeability of the core
μ_r = 1000;

"Frequency"
ω = 0



function main()

    # Build mesh
    gmsh.open(MESH_LOCATION)
    mshdata = get_mesh_data();


    # Calculate source, reluctivity and conductivity
    source_per_element = map(
        id -> source(Jp, Js, id),
        mshdata.e_group
    );
    reluctivity_per_element = map(
        id -> linear_reluctivity(μ_0, μ_r, id),
        mshdata.e_group
    );
    conductivity_per_element = map(
        id -> 0,
        mshdata.e_group
    );


    # Assemble K and f
    K, f = assemble_Kf(
        mshdata,
        source_per_element,
        reluctivity_per_element,
    );

    # Solve the system
    u = K \ f;


    # Post processing
    post(u, mshdata, source_per_element, reluctivity_per_element, conductivity_per_element, ω)

end



function post(u, mshdata, source_per_element, reluctivity_per_element, conductivity_per_element, ω)

    # Post-process for magnetic field and current density
    B, H, Wm, Jel = process(mshdata, u, source_per_element, reluctivity_per_element, conductivity_per_element, ω);

    # Define nodes (points) and elements (cells)
    points = [mshdata.xnode mshdata.ynode]';
    cells = [MeshCell(VTKCellTypes.VTK_TRIANGLE, el) for el in mshdata.elements];

    # Create VTK file structure using nodes and elements
    vtkfile = vtk_grid(string(OUTPUT_LOCATION, "steadystate1.vtu"), points, cells);

    # Store data in the VTK file
    vtkfile["Az", VTKPointData()]   = norm.(u);
    vtkfile["imA", VTKPointData()]  = imag.(u);
    vtkfile["Bnorm", VTKCellData()] = norm.(sqrt.(B[1].^2 + B[2].^2));
    vtkfile["B_vec", VTKCellData()] = real.(B)
    vtkfile["Jel", VTKCellData()]   = Jel;

    # Save the file
    print("Saving result in a file...")
    outfiles = vtk_save(vtkfile);
    println(" Done.")

end



# Execute main function
main()
