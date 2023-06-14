
include("../definitions/B.jl")

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
Solve our nonlinear problem.
"""
function nonlinear_solve(t_start, t_end, dt, u0, K, M, f, ω, mesh_data) :: Solution
    tvec = Vector(t_start:dt:t_end)
    number_of_timesteps = length(tvec)

    u = Vector{Array{ComplexF64,1}}(undef, number_of_timesteps)
    u[1] = u0

    K_t = K(0)
    α = 0.9
    u_hist = u0
    # TODO: implement LU decomposition
    # TODO: check julia slack
    start = time_ns()
    for k = 2:number_of_timesteps
        u_hist = u_hist * α + u[k-1] * (1-α)

        t = (k-1)*dt
        if (k-1) % 10 == 0
            print("  ▸ " * string(round(k/number_of_timesteps*100, digits=1)) * "% (" * string(round((time_ns() - start)/10^9, digits=2)) * " s)" * "      \r")
        end

        b_norm = B_norm(mesh_data, u_hist)
        print(b_norm)

        K_t = K(b_norm)
        u[k] = vec(
            (M + dt * K_t) \ (M * u_hist + dt * real(exp(1im * ω * t) .* f))
        )
    end

    return Solution(tvec, u)
end



"""
# Backward Euler
Solve our nonlinear problem.
"""
function nonlinear_solve2(t_start, t_end, dt, u0, nonlinear_reluctivity_function, M, f, ω, mesh_data) :: Solution
    tvec = Vector(t_start:dt:t_end)
    number_of_timesteps = length(tvec)

    u = Vector{Array{ComplexF64,1}}(undef, number_of_timesteps)
    u[1] = u0

    reluctivity = Vector{Array{Float64,1}}(undef, number_of_timesteps)
    reluctivity[1] = nonlinear_reluctivity_function(u[1])


    K_t = assemble_K(mesh_data, reluctivity[1])
    α = 0.9
    u_hist = u0

    start = time_ns()
    for k = 2:number_of_timesteps
        u_hist = u_hist * α + u[k-1] * (1-α)

        t = (k-1)*dt
        if (k-1) % 10 == 0
            print("  ▸ " * string(round(k/number_of_timesteps*100, digits=1)) * "% (" * string(round((time_ns() - start)/10^9, digits=2)) * " s)" * "      \r")
        end

        reluctivity[k] = nonlinear_reluctivity_function(u_hist)
        K_t = assemble_K(mesh_data, reluctivity[k])
        u[k] = vec(
            (M + dt * K_t) \ (M * u_hist + dt * real(exp(1im * ω * t) .* f))
        )
    end

    return Solution(tvec, u, reluctivity)
end