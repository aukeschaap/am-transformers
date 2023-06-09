

"""
# Backward Euler
Solve our time dependent problem using the backward Euler method.
"""
function solve(t_start, t_end, dt, u0, K, M, f, ω) :: Solution
    tvec = Vector(t_start:dt:t_end)
    u = Vector{Array{ComplexF64,1}}(undef, length(tvec))
    u[1] = u0
    
    # TODO: implement LU decomposition
    # TODO: check julia slack
    start = time_ns()
    for k = 2:length(tvec)
        t = (k-1)*dt
        if (k-1) % 100 == 0
            print("  ▸ " * string(round(k/length(tvec)*100)) * "% (" * string(round((time_ns() - start)/10^9, digits=2)) * " s)" * "      \r")
        end

        # Compute the solution at the next time step
        u[k] = vec(
            (M + dt * K) \ (M * u[k-1] + dt * real(exp(1im * ω * t) .* f))
        )
    end

    return Solution(tvec, u)
end


"""
# Backward Euler
Solve our time dependent problem using the backward Euler method.
"""
function nonlinear_solve(t_start, t_end, dt, u0, K, M, f, ω) :: Solution
    tvec = Vector(t_start:dt:t_end)
    number_of_timesteps = length(tvec)

    u = Vector{Array{ComplexF64,1}}(undef, number_of_timesteps)
    u[1] = u0
    
    # TODO: implement LU decomposition
    # TODO: check julia slack
    start = time_ns()
    for k = 2:number_of_timesteps
        t = (k-1)*dt
        if (k-1) % 100 == 0
            print("  ▸ " * string(round(k/number_of_timesteps*100)) * "% (" * string(round((time_ns() - start)/10^9, digits=2)) * " s)" * "      \r")
        end
        u[k] = vec(
            (M + dt * K) \ (M * u[k-1] + dt * real(exp(1im * ω * t) .* f))
        )
    end

    return Solution(tvec, u)
end