"""
# Source current density J

One term for each of the windings, with a positive and negative part. Note the phase shift between
the phases. The group_ids refer to the different physical groups the windings belong to

- Lower voltage windings: group 3 up to and including 8
- High voltage windings: group 9 up to and including 14

Arguments:
- Jp: Primary current density
- Js: Sceondary current density

Returns (implicitly):
- source current density for a specific physic group_id
"""
function source(group_id)
    Jp * exp(1im * 2pi/3) * (-1 * (group_id==3) + 1 * (group_id==4))
    + Jp * (-1 * (group_id==5) + 1 * (group_id==6))
    + Jp * exp(-1im * 2pi/3) * (-1 * (group_id==7) + 1 * (group_id==8))
    + Js * exp(1im * 2pi/3) * (1 * (group_id==9) - 1 * (group_id==10))
    + Js * (1 * (group_id==11) - 1 * (group_id==12))
    + Js * exp(-1im * 2pi/3) * (1 * (group_id==13) - 1 * (group_id==14));
end
