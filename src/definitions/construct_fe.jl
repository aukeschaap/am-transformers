"""
Constructs local contribution to source f on element e

Arguments:
- area: area of the element
- source_e: 1/u on evaluated on the element

Returns:
- floc, local contribution to global source
"""
function construct_fe(area, source)
    floc = area/3 * source * [1; 1; 1]
    return floc
end