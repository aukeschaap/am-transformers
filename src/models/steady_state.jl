
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

include("../FastSparse.jl")

include("../get_mesh_data.jl")
include("../utils/save_vtk.jl")
include("../utils/solution.jl")

include("../definitions/constants.jl")
include("../definitions/general.jl")
include("../definitions/assemble_Kf.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


function main()

    # Build mesh
    gmsh.open(MESH_LOCATION)
    mesh_data = get_mesh_data();
    println("\nMesh built.")


    # Calculate source, reluctivity and conductivity
    print("Evaluating parameters on the elements... ")
    source_per_element = map(
        id -> source(Jp, Js, id),
        mesh_data.e_group
    );
    reluctivity_per_element = map(
        id -> linear_reluctivity(μ_0, μ_r, id),
        mesh_data.e_group
    );
    conductivity_per_element = map(
        id -> conductivity(id),
        mesh_data.e_group
    );
    println("Done.")


    # Assemble linear system
    println("Linear system:")
    print("  ▸ Constructing K and f... \r")
    K, f = assemble_Kf(
        mesh_data,
        source_per_element,
        reluctivity_per_element,
    );
    println("  ✓ Constructed K and f    ")


    # Solve the system
    print("  ▸ Solving...\r")
    u = K \ real.(f);
    println("  ✓ Solved.   ")


    # Post processing
    B, H, Wm, Jel = solution(mesh_data, u, source_per_element, reluctivity_per_element, conductivity_per_element);

    print("Saving result in a file...")
    save_vtk(
        "steadystate3.vtu",
        mesh_data, u, B, H, Wm, Jel
    )
    println(" Done.")
end


# Execute main function
main()
