
using .FastSparse

include("../definitions/construct_Me.jl")



"""
# construct_M()

Assembles mass matrix, M, and right hand side, f or source, for the time-dependent solution of the z-component of the 
vector potential, A_z, derived from the Maxwell Equations in a 2D plane for vector potential A = (0,0,A_z). The linear 
systems that the solution vector u must satisfy are

1) sigma*Mu' = (1/mu)*Ku + f, or
2) (sigma*omega*j*M -(1/mu)*K)u = f,

where sigma is the conductivity, mu the permeabality, K is the stiffnes matrix and f the source vector.
Solving for u gives the coefficients u_i of the bases functions phi_i in the weighted sum that is A_z:

A_z = sum_i{u_i*phi_i}.

Arguments:
- mesh_data: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
-conductivity_per_element: 1/mu evaluated at each element in the mesh

Returns:
- M, mass matrix
"""
function assemble_M(mesh_data, conductivity_per_element)
    M = spzeros(Complex{Float64}, mesh_data.nnodes, mesh_data.nnodes)
    fsp = FastSparseMatrix(mesh_data.nelements);

    for (element_id, nodes) in enumerate(mesh_data.elements)

        # The x and y coordinates of the (triangular) element
        xs = [mesh_data.xnode[node] for node in nodes];
        ys = [mesh_data.ynode[node] for node in nodes];
        
        # Compute the area of the current element
        area = (xs[2] - xs[1])*(ys[3] - ys[1]) - (xs[3] - xs[1])*(ys[2] - ys[1]);
        area = abs(area) / 2;

        # Construct local contributions to M
        M_loc = construct_Me(area, conductivity_per_element[element_id])

        # Add local contribution to M
        add!(fsp, element_id, nodes, M_loc);
    end

    M = sparse(fsp.i_row, fsp.i_col, fsp.value, mesh_data.nnodes, mesh_data.nnodes, +);

    # Handle the boundary conditions
    bnd_node_ids, _ = gmsh.model.mesh.getNodesForPhysicalGroup(1, 1);
    M[bnd_node_ids,:] .= 0;
    M[bnd_node_ids,bnd_node_ids] = Diagonal(ones(size(bnd_node_ids)));
    return M

end