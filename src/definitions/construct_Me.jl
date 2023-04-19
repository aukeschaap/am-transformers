"""
Constructs local contribution to mass matrix M on element e

Arguments:
- area: area of the current element
- conductivity: sigma evaluated on the element

Returns:
- Bloc, local contribution to stiffnes matrix
"""
function construct_Ke(area, conductivity)
    Mloc = area / 3 * conductivity * Diagonal(ones(3));
    return Mloc
end