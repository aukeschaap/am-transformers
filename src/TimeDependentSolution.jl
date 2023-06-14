module TimeDependentSolution

export Solution

"""
# Solution

A struct that contains the data for a solution of a time dependent problem.
"""
mutable struct Solution
    t::Vector{Float64}
    u::Vector{Array{ComplexF64,1}}
    reluctivity::Vector{Array{ComplexF64,1}}
end


end # module TimeDependentSolution