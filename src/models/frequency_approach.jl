using Logging: global_logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

try
    using Gmsh: gmsh
catch
    using gmsh
end


using LinearAlgebra, SparseArrays, WriteVTK, BenchmarkTools, Plots, Printf

include("../FastSparse.jl")

include("../get_mesh_data.jl")
include("../utils/save_vtk_series.jl")
include("../utils/solution.jl")

include("../definitions/constants.jl")
include("../definitions/general.jl")
include("../definitions/assemble_Kf.jl")
include("../definitions/assemble_M.jl")


const CLEAR_MESH_DATA = false
const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


# Build mesh
if CLEAR_MESH_DATA == true || (@isdefined mesh_data) == false
    gmsh.finalize()
    gmsh.initialize()
    gmsh.open(MESH_LOCATION)
    mesh_data = get_mesh_data();
    println("\nMesh built.")
else
    println("\nReused built mesh.")
end


function main()

    # set frequency of source current
    freq = 500 # Hz
    ω = 2π*freq

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

    # Solve linear system
    u = (K + 1im*ω*M) \ f

    # redfine source function to include frequency
    function construct_Je(c, conductivity, source)
        return norm(source + ω * conductivity * 1/3 * sum(c));
    end

    B, H, Wm, Jel = solution(mesh_data, u, source_per_element, reluctivity_per_element, conductivity_per_element)

    # Save solution
    vtk_path = "frequency_approach_" * string(round(freq,digits=2)) * ".vtu";
    save_vtk(vtk_path, mesh_data, u, B, Jel, reluctivity_per_element)

end

main()