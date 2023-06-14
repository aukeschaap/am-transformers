

"""
Constructs the components of the B- & H-fields, the energy and current density given the solution vector, u.

Arguments:
- mesh_data: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
- u: solution vector

Returns:
- Bx, By, Bz: x-, y- and z-components of the B-field. 
"""
function solution(mesh_data, u, source_per_element, reluctivity_per_element, conductivity_per_element)
    Bx = zeros(Complex{Float64}, mesh_data.nelements);
    By = zeros(Complex{Float64}, mesh_data.nelements);
    Bz = zeros(Complex{Float64}, mesh_data.nelements);
    Jel = zeros(mesh_data.nelements);

    for (element_id, nodes) in enumerate(mesh_data.elements)

        # nodal coordinates
        xs(i) = mesh_data.xnode[nodes[i]];
        ys(i) = mesh_data.ynode[nodes[i]];

        # solution coefficients
        c = u[nodes[1:3]];

        # B components
        Bx[element_id], By[element_id] = construct_Be(c, xs, ys)

        # current 
        Jel[element_id] = construct_Je(c, source_per_element[element_id], conductivity_per_element[element_id])
    end

    # H is related to B through the reluctivity
    Hx = reluctivity_per_element' .* Bx;
    Hy = reluctivity_per_element' .* By;
    
    # Energy is 0.5 * dot(B, H)
    Wm = 0.5 * (Bx .* Hx .+ By .* Hy);
    
    return (Bx,By,Bz), (Hx, Hy), Wm, Jel;
    
end


function nonlinear_solution(mesh_data, u, source_per_element, nonlinear_reluctivity_per_element, conductivity_per_element)
    Bx = zeros(Complex{Float64}, mesh_data.nelements);
    By = zeros(Complex{Float64}, mesh_data.nelements);
    Bz = zeros(Complex{Float64}, mesh_data.nelements);
    Jel = zeros(mesh_data.nelements);

    for (element_id, nodes) in enumerate(mesh_data.elements)

        # nodal coordinates
        xs(i) = mesh_data.xnode[nodes[i]];
        ys(i) = mesh_data.ynode[nodes[i]];

        # solution coefficients
        c = u[nodes[1:3]];

        # B components
        Bx[element_id], By[element_id] = construct_Be(c, xs, ys)

        # current 
        Jel[element_id] = construct_Je(c, source_per_element[element_id], conductivity_per_element[element_id])
    end

    # H is related to B through the reluctivity
    B_norm = real(.√(Bx.^2 + By.^2))
    rel = nonlinear_reluctivity_per_element(B_norm)
    perm = 1 ./ rel
    Hx = rel' .* Bx;
    Hy = rel' .* By;
    
    # Energy is 0.5 * dot(B, H)
    Wm = 0.5 * (Bx .* Hx .+ By .* Hy);
    
    return (Bx,By,Bz), (Hx, Hy), Wm, Jel, rel, perm;
    
end


"""
Construct the B components
"""
function construct_Be(c, xs, ys)
    # construct shape function
    Emat = [
        xs(1:3) ys(1:3) [1, 1, 1]
    ] \ UniformScaling(1.);

    # Calculate Bx and By from the solution coefficients and the shape function parameters
    Bx = sum(c .* Emat[2,:]);
    By = -sum(c .* Emat[1,:]);

    return Bx, By
end


"""
Calculate eddy current loss.
"""
function construct_Je(c, conductivity, source)
    return norm(source + conductivity * 1/3 * sum(c));
end
