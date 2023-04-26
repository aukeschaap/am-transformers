"""
Constructs local contribution to source f on element e

Arguments:
- area: area of the element
- source_e: 1/u on evaluated on the element

Returns:
- floc, local contribution to global source
"""
function construct_fe(area, source)
    fe = area/3 * source
    return [fe; fe; fe;]
end