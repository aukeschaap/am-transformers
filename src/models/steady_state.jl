
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

include("../constants.jl")
include("../get_mesh_data.jl")
include("../utils/process.jl")
include("../utils/save.jl")
include("../utils/result.jl")

include("../definitions/general.jl")
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
    print("Evaluating parameters on the elements... ")
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
    println("Done.")


    # Assemble K and f
    println("Linear system:")
    print("  ▸ Constructing K and f... \r")
    K, f = assemble_Kf(
        mshdata,
        source_per_element,
        reluctivity_per_element,
    );
    println("  ✓ Constructed K and f    ")


    # Solve the system
    print("  ▸ Solving...\r")
    u = K \ f;
    println("  ✓ Solved.   ")


    # Post processing
    B, H, Wm, Jel = process(mshdata, u, source_per_element, reluctivity_per_element, conductivity_per_element, ω);

    save(
        "steadystate2.vtu",
        mshdata, u, B, H, Wm, Jel
    )

end



# Execute main function
main()
