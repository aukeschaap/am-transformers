

"""
Post processing.
"""
function process(mesh_data, u, sourceperelement, reluctivityperelement, conductivityperelement, omega)
    Bx = zeros(Complex{Float64}, mesh_data.nelements);
    By = zeros(Complex{Float64}, mesh_data.nelements);
    Bz = zeros(Complex{Float64}, mesh_data.nelements);
    
    Jel = zeros(mesh_data.nelements);
    
    # Perform a loop over the elements
    for (element_id, nodes) in enumerate(mesh_data.elements)
        # Retrieve global numbering of the local nodes of the current element
        node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];
        
        # Get x and y coordinates of the three nodes
        xnode = mesh_data.xnode; ynode = mesh_data.ynode;
        xnode1 = xnode[node1_id]; xnode2 = xnode[node2_id]; xnode3 = xnode[node3_id];
        ynode1 = ynode[node1_id]; ynode2 = ynode[node2_id]; ynode3 = ynode[node3_id];

        # Compute surface area of the current element
        x12 = xnode2 - xnode1; x13 = xnode3-xnode1;
        y12 = ynode2 - ynode1; y13 = ynode3-ynode1;
        area_id = x12*y13 - x13*y12; area_id = abs(area_id)/2

        # Calculate shape function parameters
        Emat = [[xnode1;xnode2;xnode3] [ynode1;ynode2;ynode3] [1;1;1]] \ UniformScaling(1.);
    
        # Calculate Bx and By from the solution coefficients and the shape function parameters
        c = u[[node1_id, node2_id, node3_id]];
        Bx[element_id] = sum(c .* Emat[2,:]);
        By[element_id] = -sum(c .* Emat[1,:]);
        
        # Calculate eddy current loss
        sigma = conductivityperelement[element_id];
        Jel[element_id] = norm(sourceperelement[element_id] + omega * sigma * 1/3 * sum(c));
    end
    
    # H is related to B through the reluctivity
    Hx = reluctivityperelement' .* Bx;
    Hy = reluctivityperelement' .* By;
    
    # Energy is 0.5 * dot(B, H)
    Wm = 0.5 * (Bx .* Hx .+ By .* Hy);
    
    return (Bx,By,Bz), (Hx, Hy), Wm, Jel;
end
