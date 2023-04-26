"""
# conductivity(group_id)
Defines the conductivity of the iron core. By Max's thesis, the core's conductivity is taken to be constant throughout 
the entire core and is reduced by a factor of 10 due to eddy currents. 

Arguments:
- group_id: the physical group id

Returns:
value of the conductivity
"""
function conductivity(group_id)
    0.1
end