
include("../definitions/construct_fe.jl")
include("../definitions/construct_Ke.jl")


mutable struct FastSparse
    i_row::Vector{Int}
    i_col::Vector{Int}
    value::Vector{Complex{Float64}}
    FastSparse(nelements::Int) = new(Vector{Int}(undef, 9*nelements), Vector{Int}(undef, 9*nelements), Vector{Complex{Float64}}(undef, 9*nelements))
end


function add!(fsp::FastSparse, id::Int, nodes::Vector{Int}, matrix::Matrix{Float64})
    @debug "add: " id nodes matrix
    for (num, (index, element)) in enumerate(pairs(matrix))
        (i, j) = Tuple(index)
        @debug "    ($i, $j) @ ($((id-1)*9 + num)) = $element)"
        fsp.i_row[(id-1)*9 + num] = nodes[i]
        fsp.i_col[(id-1)*9 + num] = nodes[j]
        fsp.value[(id-1)*9 + num] = element
    end
    return nothing
end



"""
# assemble_steadystate()

Assembles stiffnes matrix, K, and right hand side, f or source, for the steady state solution of the z-component of the 
vector potential, A_z, derived from the Maxwell Equations in a 2D plane for vector potential A = (0,0,A_z). Solving 
the resulting linear system, Ku=f, for u gives the coefficients u_i of the bases functions phi_i in the weighted 
sum that is A_z:

A_z = sum_i{u_i*phi_i}.

Arguments:
- mesh_hdata: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
- sourceperelement: f evaluated at each element in the mesh
- reluctivityperelement: 1/mu evaluated at each element in the mesh

Returns:
- K, stiffnes matrix
- f, source vector
"""
function assemble_Kf(mesh_data, source_per_element, reluctivity_per_element)
    f = zeros(Complex{Float64}, mesh_data.nnodes, 1);
    fsp = FastSparse(mesh_data.nelements);

    for (element_id, nodes) in enumerate(mesh_data.elements)

        # The x and y coordinates of the (triangular) element
        xs(i) = mesh_data.xnode[nodes[i]];
        ys(i) = mesh_data.ynode[nodes[i]];
        
        # Compute the area of the current element
        area = (xs(2) - xs(1))*(ys(3) - ys(1)) - (xs(3) - xs(1))*(ys(2) - ys(1));
        area = abs(area) / 2;

        # Construct local contributions to f and K
        f_loc = construct_fe(area, source_per_element[element_id]);
        K_loc = construct_Ke(xs, ys, area, reluctivity_per_element[element_id]);

        # Add local contribution to K & f
        add!(fsp, element_id, nodes, K_loc);

        f[nodes] += f_loc;

    end

    K = sparse(fsp.i_row, fsp.i_col, fsp.value, mesh_data.nnodes, mesh_data.nnodes, +);

    # Handle the boundary conditions
    bnd_node_ids, _ = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1);
    K[bnd_node_ids,:] .= 0;
    K[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)));
    f[bnd_node_ids] .= 0;

    return K, f
end
