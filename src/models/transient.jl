
try
    using Gmsh: gmsh
catch
    using gmsh
end


using LinearAlgebra, SparseArrays, WriteVTK, BenchmarkTools
using Logging

include("../FastSparse.jl")

include("../get_mesh_data.jl")
include("../utils/process.jl")
include("../utils/save.jl")
include("../utils/result.jl")

include("../definitions/constants.jl")
include("../definitions/general.jl")
include("../definitions/assemble_Kf.jl")
include("../definitions/assemble_M.jl")


const MESH_LOCATION = "./mesh/transformer_stedin.msh"
const OUTPUT_LOCATION = "./out/"


function main()

    # Build mesh
    gmsh.open(MESH_LOCATION)
    mshdata = get_mesh_data();
    println("\nMesh built.")

    # Calculate source, reluctivity and conductivity
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
    println("Done.")

    # Assemble linear system
    println("Linear system:")
    print("  ▸ Constructing K and f...\r")
    K, f = assemble_Kf(mshdata,source_per_element,reluctivity_per_element)
    println("  ✓ Constructed K and f    ")
    print("  ▸ Constructing M...\r")
    M = assemble_M(mshdata, conductivity_per_element)
    println("  ✓ Constructed M    ")



    # etc...


    B, H, Wm, Jel = solution(mshdata, u, source_per_element, reluctivity_per_element, conductivity_per_element);

    save(
        "transient1.vtu",
        mshdata, u, B, H, Wm, Jel
    )

end


# Execute main function
main()
