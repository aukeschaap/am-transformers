"""
Constructs the components of the B-field given the solution vector, u.

Arguments:
- mshdata: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
- u: solution vector

Returns:
- Bx, By: x- and y-components of the B-field. 
"""
function process(mshdata, u)
    Bx = zeros(Complex{Float64}, mshdata.nelements);
    By = zeros(Complex{Float64}, mshdata.nelements);
    xnode = mshdata.xnode; ynode = mshdata.ynode;

    # Perform a loop over the elements
    for (element_id, nodes) in enumerate(mshdata.elements)
        # Retrieve global numbering of the local nodes of the current element
        node1_id = nodes[1]; node2_id = nodes[2]; node3_id = nodes[3];
        
        # Get x and y coordinates of the three nodes
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
    return Bx, By