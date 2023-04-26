"""
Constructs local contribution to stiffnes matrix K on element e

Arguments:
- xs/ys: contain nodal x- and y-coordinates of the current element
- area: area of the current element
- reluctivity: 1/mu (with mu being the permeabality) evaluated on the element

Returns:
- Bloc, local contribution to stiffnes matrix
"""
function construct_Ke(xs, ys, area, reluctivity)
    Emat = [
        xs(1:3) ys(1:3) [1, 1, 1]
    ] \ UniformScaling(1.);
    Emat[3,:] .= 0;
    return area*reluctivity*(transpose(Emat)*Emat);
end
