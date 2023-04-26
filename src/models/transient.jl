using gmsh
using LinearAlgebra, SparseArrays, WriteVTK, BenchmarkTools
using Logging
include("../constants.jl")
include("../get_mesh_data.jl")
include("../process.jl")
include("../definitions/source.jl")
include("../definitions/linear_reluctivity.jl")
include("../definitions/conductivity.jl")
include("../definitions/assemble_Kf.jl")
include("../definitions/assemble_M.jl")
include("../definitions/assemble_solution.jl")

const MESH_LOCATION = "../../mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "../../out/"
gmsh.open(MESH_LOCATION)
mshdata = get_mesh_data();

# peak phase current
Ip = 0;
Is = 777.62;

# number of windings
Np = 266;
Ns = 6;

# calculate current density in the windings
Jp = Np * Ip / Awhv;
Js = Ns * Is / Awlv;

# vacuum & relative permeability
μ_0 = 4e-7 * pi;
μ_r = 1000;    


print("Evaluating parameters on the elements...")
source_per_element = map(
    id -> source(Jp, Js, id),
    mshdata.e_group
);

reluctivity_per_element = map(
    id -> linear_reluctivity(μ_0, μ_r, id),
    mshdata.e_group
);

conductivity_per_element = map(
    id -> conductivity(id),
    mshdata.e_group
);
print("Finished evaluating parameters.")

print("Constructing linear system...")
K, f = @timev assemble_Kf(mshdata,source_per_element,reluctivity_per_element)
print("Finished construction of K & f.")
M = @timev assemble_M(mshdata, conductivity_per_element)
print("Finished construction of M.")

B, H, Wm, Jel = assemble_solution(mshdata, u, source_per_element, reluctivity_per_element, conductivity_per_element);
Bnorm = norm.(sqrt.(B[1].^2 + B[2].^2));