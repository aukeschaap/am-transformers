
include("../definitions/B.jl")

"""
# Backward Euler
Solve our time dependent problem using the backward Euler method.
"""
function linear_solve(t_start, t_end, dt, u0, K, M, f, ω) :: Solution
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

TODO: implement LU decomposition
TODO: check julia slack
"""
function nonlinear_solve_wrong(t_start, t_end, dt, u0, K, M, f, ω, mesh_data) :: Solution
    tvec = Vector(t_start:dt:t_end)
    number_of_timesteps = length(tvec)

    u = Vector{Array{ComplexF64,1}}(undef, number_of_timesteps)
    u[1] = u0

    K_t = K(0)
    α = 0.9
    u_hist = u0
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

TODO: implement LU decomposition
TODO: check julia slack
"""
function nonlinear_solve(t_start, t_end, dt, u0, nonlinear_reluctivity_function, M, f, ω, mesh_data) :: Solution
    tvec = Vector(t_start:dt:t_end)
    number_of_timesteps = length(tvec)

    u = Vector{Array{ComplexF64,1}}(undef, number_of_timesteps)
    u[1] = u0


    start = time_ns()
    for k = 2:number_of_timesteps
        _display_progress(k, number_of_timesteps, start)

        u[k] = _solve(u[k-1], tvec[k], dt, M, f, ω, nonlinear_reluctivity_function)
    end

    return Solution(tvec, u)
end


function _solve(u, t, dt, M, f, ω, nonlinear_reluctivity_function)
    NUM_ITERATIONS = 4
    α = 0.9

    reluctivity = nonlinear_reluctivity_function(u) # prev reluctivity
    K_t = assemble_K(mesh_data, reluctivity)
    u = vec(
        (M + dt * K_t) \ (M * u + dt * real(exp(1im * ω * t) .* f))
    )

    u_hist = u
    u_prev = u

    for i = 1:NUM_ITERATIONS
        u_prev = u;
        u_hist = u_hist * α + u_prev * (1-α)

        reluctivity = nonlinear_reluctivity_function(u_hist)
        K_t = assemble_K(mesh_data, reluctivity)
        u = vec(
            (M + dt * K_t) \ (M * u_hist + dt * real(exp(1im * ω * t) .* f))
        )
    end
    # println("\n acc: " * string(norm(u-u_prev)) * "\n")

    return u
end


function _display_progress(k, number_of_timesteps, start)
    if (k-1) % 10 == 0
        progress = round(k/number_of_timesteps*100, digits=1)
        elapsed = round((time_ns() - start)/10^9, digits=2)
        estimated = round(number_of_timesteps * elapsed / k, digits=2)
        print(
            "  ▸ " * string(progress) * 
            "% (" * string(elapsed) * 
            " of est. " * string(estimated) *
            " s)        \r"
        )
    end
end