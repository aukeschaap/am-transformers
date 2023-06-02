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
include("../definitions/assemble_sys.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


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
    # println("Linear system:")
    # print("  ▸ Constructing K and f...\r")
    # K, f = assemble_Kf(mesh_data,source_per_element,reluctivity_per_element)
    # println("  ✓ Constructed K and f    ")
    # print("  ▸ Constructing M...\r")
    # M = assemble_M(mesh_data, conductivity_per_element)
    # println("  ✓ Constructed M    ")
    println("Linear system:")
    print("  ▸ Constructing FE system...\r")
    K, M, f = assemble_sys(mesh_data, source_per_element, reluctivity_per_element, conductivity_per_element)
    println("  ✓ Constructed FE system    ")

    # time stepping: source is causing instability!
    v = Vector{Float64}(undef, mesh_data.nnodes)
    function magneticVectorPotentialEquation!(du,u,p,t)
        f_floats = reinterpret(Float64, exp(1im*ω*t).*f)
        f_real = @view f_floats[1:2:end-1]
        du .= f_real .+ mul!(v,K,u) #real.
        # du .= M \ (f_real .+ mul!(v,K,u)) #real.
    end
    f_func = ODEFunction(magneticVectorPotentialEquation!, mass_matrix = M)
    
    # set initial condition
    u0 = fill(0., mesh_data.nnodes)
                                        
    # set time begin and end
    number_of_periods = 1
    t0 = 0.0
    tf = number_of_periods*(2*pi/ω)
    tspan = (t0, tf)
    dt = (tf-t0)/30
    tvec = Vector(t0:dt:tf)

    # define ODE problem to be solved  +
    println("ODE problem:")
    print("  ▸ Defining ODE...\r")
    prob_magneticVectorPotential = ODEProblem(f_func, u0, tspan)
    println("  ✓ ODE defined                ")
    print("  ▸ solving ODE...\r")
    # sol = solve(prob_magneticVectorPotential, Euler(), dt=dt, force_dtmin = false, progress=true, progress_steps=10);
    sol = solve(prob_magneticVectorPotential, ROS3P(autodiff=false), progress=true, progress_steps=10)
    println("  ✓ ODE solved    ")
    
    # check initial solution
    println("Checking initial solution...")
    B, H, Wm, Jel = solution(mesh_data, sol(0.0), source_per_element, reluctivity_per_element, conductivity_per_element);
    Bnorm = sqrt.(B[1].^2 + B[2].^2)
    println(" → number of nodes: ", mesh_data.nnodes)
    println(" → number of xnodes: ", length(mesh_data.xnode))
    println(" → number of ynodes: ", length(mesh_data.ynode))
    println(" → length of Bnorm: ", length(Bnorm))
    println(" → number of elements: ", mesh_data.nelements)
    n_time_steps = length(sol.t)
    println(" → number of time steps: ", n_time_steps)

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
