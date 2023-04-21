
"""
Linear reluctivity of the different materials.

"""
function linear_reluctivity(μ_0, μ_r, id)
    (1 / μ_0) + (1/(μ_0*μ_r) - 1/μ_0) * (id == 2)
end
