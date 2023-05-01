include("construct_Be.jl")
include("construct_Je.jl")

"""
Constructs the components of the B-field given the solution vector, u.

Arguments:
- mshdata: contains all mesh information: (number of) elements & nodes (global & local) & coordinates.
- u: solution vector

Returns:
- Bx, By, Bz: x-, y- and z-components of the B-field. 
"""
function assemble_solution(mshdata, u, source_per_element, reluctivity_per_element, conductivityperelement,)
    Bx = zeros(Complex{Float64}, mshdata.nelements);
    By = zeros(Complex{Float64}, mshdata.nelements);
    Bz = zeros(Complex{Float64}, mshdata.nelements);
    Jel = zeros(mshdata.nelements);

    for (element_id, nodes) in enumerate(mshdata.elements)

        # nodal coordinates
        xs(i) = mesh_data.xnode[nodes[i]];
        ys(i) = mesh_data.ynode[nodes[i]];

        # solution coefficients
        c = u[nodes(1:3)];

        # B components
        Bx[element_id],By[element_id] = construct_Be(c, xs, ys)

        # current 
        Je[element_id] = construct_Je(c, sourceperelement[element_id], conductivityperelement[element_id])
    end

    # H is related to B through the reluctivity
    Hx = reluctivityperelement' .* Bx;
    Hy = reluctivityperelement' .* By;
    
    # Energy is 0.5 * dot(B, H)
    Wm = 0.5 * (Bx .* Hx .+ By .* Hy);
    
    return (Bx,By,Bz), (Hx, Hy), Wm, Jel;
    
end