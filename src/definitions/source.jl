
"""
# Source current density J

One term for each of the windings, with a positive and negative part. Note the typical three phase current:
-2*pi/3, 0 , 2*pi/3. The ids refer to the different (physical) groups the windings belong to.

- Lower voltage windings: group 3 up to and including 8
- High voltage windings: group 9 up to and including 14

Furthermore, 
- Lower voltage windings: even/odd ids represent current flowing into/out of the plane
- Higher voltage windings: even/odd ids represent current flowing out of/into the plane

Arguments:
- Jp: Primary current density.
- Js: Sceondary current density.
- id: Group id of the physical group the winding belongs to. 

Returns (implicitly):
- source current density for a specific physic group_id
"""
function source(Jp, Js, id)
    (
        Jp * (
            exp(1im * 2pi/3) * (
                (id==4) - (id==3)
            ) + (
                (id==6) - (id==5)
            ) + exp(-1im * 2pi/3) * (
                (id==8) - (id==7)
            )
        )

        + 

        Js * (
            exp(1im * 2pi/3) * (
                (id==9) - (id==10)
            )
            + (
                (id==11) - (id==12)
            )
            + exp(-1im * 2pi/3) * (
                (id==13) - (id==14)
            )
        )
    );
end
