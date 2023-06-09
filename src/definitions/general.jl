
"""
# Conductivity

The `conductivity` function defines the conductivity of the iron core. By Max's thesis,
the core's conductivity is taken to be constant throughout the entire core and is reduced by a
factor of 10 to model the effect of eddy currents.

Arguments:
- group_id: the physical group id

Returns:
value of the conductivity
"""
function conductivity(id)
    if id==2
        0.1
    else
        1/1000
    end
end



"""
# Linear reluctivity

The linear reluctivity of the different materials.

Arguments:
- μ_0: permeability of free space
- μ_r: relative permeability of the material
- id: the group id of the physical group the element belongs to. For id==2 the permeability of the iron core is returned.
"""
function linear_reluctivity(μ_0, μ_r, id)
    (1/μ_0) + (1/(μ_0*μ_r) - 1/μ_0) * (id == 2)
end



"""
# Nonlinear reluctivity

The nonlinear reluctivity of the different materials.

Arguments:
- μ_0: permeability of free space
- μ_r: relative permeability of the material
- id: the group id of the physical group the element belongs to. For id==2 the permeability of the iron core is returned.
"""
function nonlinear_reluctivity(id, B)
    bh_a = 2.12e-4; 
    bh_b = 7.358;
    bh_c = 1.18e7;

    μ_r_core(B) = 1 / (bh_a + (1 - bh_a) * B^(2*bh_b) / (B^(2*bh_b) + bh_c));
    
    return 1 / (
        μ_0 + 
        μ_0 * (id == 2) * (μ_r_core(B) - 1)
    )
end

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
- Js: Secondary current density.
- id: Group id of the physical group the winding belongs to. 

Returns (implicitly):
- source current density for a specific physic group_id
"""
function source(Jp, Js, id)
    (
        Jp * (
            (exp(1im * 2pi/3) * ((id==4) - (id==3)))
            +
            ((id==6) - (id==5)) 
            + 
            (exp(-1im * 2pi/3) * ((id==8) - (id==7)))
        )

        + 

        Js * (
            (exp(1im * 2pi/3) * ((id==9) - (id==10)))
            + 
            ((id==11) - (id==12)) 
            + 
            (exp(-1im * 2pi/3) * ((id==13) - (id==14)))
        )
    );
end
