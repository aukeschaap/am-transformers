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
# include("../definitions/assemble_sys.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"

# create solution type
struct sol_type
    t::Array{Float64,1}
    u::Array{Array{ComplexF64,1},1}
end

function main()

    # Build mesh
    gmsh.open(MESH_LOCATION)
    mesh_data = get_mesh_data();
    println("\nMesh built.")

    # set frequency of source current
    freq = 50 # Hz;
    ω = 2π*freq;

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
    K, f = assemble_Kf(mesh_data, source_per_element, reluctivity_per_element)
    println("  ✓ Constructed K and f    ")
    print("  ▸ Constructing M...\r")
    M = assemble_M(mesh_data, conductivity_per_element)
    println("  ✓ Constructed M    ")
    println("  ▸ check singularity of M...")
    M_size = size(M)
    M_rank = rank(M)
    println("       shape of M: ", M_size)
    println("       rank of M: ", M_rank)
    if M_size[1] == M_rank
        println("  ✓ M is non-singular")
    elseif M_size[1] > M_rank
        println("  ✗ M is singular")    
    end

    # set time begin and end
    number_of_periods = 5
    t0 = 0.0
    tf = number_of_periods * (2*pi/ω)
    tspan = (t0, tf)
    dt = (tf-t0)/1000
    tvec = Vector(t0:dt:tf)
    Nt = length(tvec)

     # set initial condition
    u0 = fill(0., mesh_data.nnodes)
    uk = copy(u0)
    u = Array{Array{ComplexF64,1},1}(undef, Nt)
    u[1:Nt] = [copy(u0) for i in 1:Nt]

    # perform time integration using backward Euler  
    println("Solving ODE...")
    for k = 2:Nt
        t = (k-1)*dt
        uk = (M + dt .* K) \ (M * uk + dt .* real.(exp(1im * ω * t) .* f))
        u[k] = copy(uk[:,1])
        message = "  ▸ " * string(round(k/Nt*100,digits=2)) * "%" * " time = " * string(round(10^6*t, digits=2)) * " μs" * "\r"
        print(message)
    end 
    sol = sol_type(tvec, u)
    print("                                                             \r")
    println("  ✓ ODE solved    ")
    
    # check initial solution
    println("Checking initial solution...")
    B, H, Wm, Jel = solution(mesh_data, u[1], source_per_element, reluctivity_per_element, conductivity_per_element);
    Bnorm = sqrt.(B[1].^2 + B[2].^2)
    println(" → number of nodes: ", mesh_data.nnodes)
    println(" → number of xnodes: ", length(mesh_data.xnode))
    println(" → number of ynodes: ", length(mesh_data.ynode))
    println(" → length of Bnorm: ", length(Bnorm))
    println(" → number of elements: ", mesh_data.nelements)
    n_time_steps = length(tvec)
    println(" → number of time steps: ", n_time_steps)

    println("saving solution...")
    # create solution type

    # save time series
    print("  ▸ saving time series...\r")
    save_vtk_series(
        "transient", 
        mesh_data, 
        sol,
        source_per_element,
        reluctivity_per_element,
        conductivity_per_element
    )
    println("  ✓ time series saved    ")
end

# Execute main function
main()
