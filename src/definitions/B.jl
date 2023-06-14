

function B_norm(mesh_data, u)
    Bx = zeros(Complex{Float64}, mesh_data.nelements);
    By = zeros(Complex{Float64}, mesh_data.nelements);
    Bz = zeros(Complex{Float64}, mesh_data.nelements);

    for (element_id, nodes) in enumerate(mesh_data.elements)

        # nodal coordinates
        xs(i) = mesh_data.xnode[nodes[i]];
        ys(i) = mesh_data.ynode[nodes[i]];

        # solution coefficients
        c = u[nodes[1:3]];

        # B components
        Bx[element_id], By[element_id] = construct_Be(c, xs, ys)
    end
    
    return real(.√(Bx.^2 + By.^2))
end