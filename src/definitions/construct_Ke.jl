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
        [xs[1] ys[1] 1]
        [xs[2] ys[2] 1]
        [xs[3] ys[3] 1]
    ] \ UniformScaling(1.);
    Emat[3,:] .= 0;
    Kloc = area*reluctivity*(transpose(Emat)*Emat);
    return Kloc
end