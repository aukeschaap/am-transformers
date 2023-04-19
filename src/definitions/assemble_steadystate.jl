"""
# assemble_steadystate()

Assembles stiffnes matrix, K, and right hand side, f or source, for the steady state solution of the z-component of the 
vector potential, A_z, derived from the Maxwell Equations in a 2D plane for vector potential A = (0,0,A_z). Solving 
the resulting linear system, Ku=f, for u gives the coefficients u_i of the bases functions phi_i in the weighted 
sum that is A_z:

A_z = sum_i{u_i*phi_i}.

Arguments:
- mshdata: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
- sourceperelement: f evaluated at each element
- reluctivityperelement: 1/mu evaluated at each element

Returns:
- K, stiffnes matrix
- f, source function
"""
function assemble_steadystate(mshdata, sourceperelement, reluctivityperelement)
    K = spzeros(Complex{Float64}, mshdata.nnodes, mshdata.nnodes)
    f = zeros(Complex{Float64}, mshdata.nnodes, 1)

    xnode = mshdata.xnode;
    ynode = mshdata.ynode;

    # Perform a loop over the elements
    for (element_id, nodes) in enumerate(mshdata.elements)
        # Retrieve global numbering of the local nodes of the current element
        node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];

        # Retrieve the x and y coordinates of the local nodes of the current element
        xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
        ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

        # Compute surface area of the current element
        x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
        y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
        area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

        # Compute local vector contribution floc of the current element
        floc = area_id/3*sourceperelement[element_id] * [1; 1; 1]

        # Compute local matrix contribution Aloc of the current element
        Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
        Emat[3,:] .= 0;
        Bloc = area_id*reluctivityperelement[element_id]*(transpose(Emat)*Emat);

        # Add local contribution to f and K
        f[nodes]       += floc;
        K[nodes,nodes] += Bloc;

        return K, f
    end