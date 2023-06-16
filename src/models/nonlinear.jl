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
include("../TimeDependentSolution.jl")

include("../get_mesh_data.jl")
include("../utils/save_vtk_series.jl")
include("../utils/solution.jl")


include("../definitions/constants.jl")
include("../definitions/general.jl")
include("../definitions/assemble_Kf.jl")
include("../definitions/assemble_M.jl")
include("../definitions/backward_euler.jl")

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


    # Frequency of source current
    freq = 50 # Hz
    ω = 2π*freq;

    # Calculate source, reluctivity and conductivity
    print("Evaluating parameters on the elements...")
    source_per_element = map(
        id -> source(Jp, Js, id),
        mesh_data.e_group
    );
    nonlinear_reluctivity_per_element(u_k) = map(
        (tup) -> nonlinear_reluctivity(tup[1], tup[2]),
        zip(mesh_data.e_group, B_norm(mesh_data, u_k))
    );
    conductivity_per_element = map(
        id -> conductivity(id),
        mesh_data.e_group
    );
    println("Done.")

    # Assemble constant mass matrix
    println("Linear system:")
    print("  ▸ Constructing f...\r")
    f = assemble_f(mesh_data, source_per_element)
    println("  ✓ Constructed f                   ")
    print("  ▸ Constructing M...\r")
    M = assemble_M(mesh_data, conductivity_per_element)
    println("  ✓ Constructed M    ")
    println("  ▸ Checking singularity of M...")
    M_size = size(M)
    M_rank = rank(M)
    println("       shape of M: ", M_size)
    println("       rank of M: ", M_rank)
    if M_size[1] == M_rank
        println("  ✓ M is non-singular")
    elseif M_size[1] > M_rank
        println("  ✗ M is singular")
    end


    # Specify time start, end and step
    t_0 = 0.0
    T   = 2(2pi/ω)
    dt  = (T-t_0) / 400
    println("Time discretization:")
    println("  ▸ t_0 = ", t_0)
    println("  ▸ T   = ", T)
    println("  ▸ dt  = ", dt)
    println("  ▸ N   = ", Int(round((T-t_0)/dt)))

    # Initial condition
    u0 = zeros(mesh_data.nnodes)

    # Perform time integration using Backward Euler
    println("Solving ODE...")
    sol = nonlinear_solve(t_0, T, dt, u0, nonlinear_reluctivity_per_element, M, f, ω, mesh_data)
    println("  ✓ ODE solved." *
            "                                    ")
    
    # Save time series
    print("  ▸ Saving time series...\r")
    solution_function(mesh_data, u) = nonlinear_solution(mesh_data, u, source_per_element, nonlinear_reluctivity_per_element, conductivity_per_element)
    folder_name = save_vtk_series(
        "transient",
        mesh_data, 
        sol,
        solution_function
    )
    println("  ✓ Time series saved as '" * string(folder_name) * "'")
end

# Execute main function
main()
