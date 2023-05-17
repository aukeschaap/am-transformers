using Logging: global_logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

try
    using Gmsh: gmsh
catch
    using gmsh
end


using LinearAlgebra, SparseArrays, WriteVTK, BenchmarkTools, DifferentialEquations, Plots, Printf

include("../FastSparse.jl")

include("../get_mesh_data.jl")
include("../utils/save_vtk_series.jl")
include("../utils/solution.jl")

include("../definitions/constants.jl")
include("../definitions/general.jl")
include("../definitions/assemble_Kf.jl")
include("../definitions/assemble_M.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


function main()

    # Build mesh
    gmsh.open(MESH_LOCATION)
    mesh_data = get_mesh_data();
    println("\nMesh built.")

    # set frequency of source current
    ω = 2π*0.1

    # Calculate source, reluctivity and conductivity (linear element)
    print("Evaluating parameters on the elements...")
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

    # Assemble constant mass matrix
    println("Linear system:")
    print("  ▸ Constructing K and f...\r")
    K, f = assemble_Kf(mesh_data,source_per_element,reluctivity_per_element)
    println("  ✓ Constructed K and f    ")
    print("  ▸ Constructing M...\r")
    M = assemble_M(mesh_data, conductivity_per_element)
    println("  ✓ Constructed M    ")

    # source
    display(f)    

    # Solve linear system
    u = K \ f

    # redfine source function to include frequency
    function construct_Je(c, conductivity, source)
        return norm(source + ω * conductivity * 1/3 * sum(c));
    end

    B, H, Wm, Jel = solution(mesh_data, u, source_per_element, reluctivity_per_element, conductivity_per_element)

    # Save solution
    vtk_path = "frequency_approach.vtu"
    save_vtk(vtk_path, mesh_data, u, B, H, Wm, Jel)

end

main()