using Logging: global_logger
using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using DifferentialEquations, Plots

function lorenz!(du, u, p, t)
    du[1] = 10.0 * (u[2] - u[1])
    du[2] = u[1] * (28.0 - u[3]) - u[2]
    du[3] = u[1] * u[2] - (8 / 3) * u[3]
end

u0 = [1.0; 0.0; 0.0]
tspan = (0.0, 100.0)
prob = ODEProblem(lorenz!, u0, tspan, progress=true, progress_steps=1)
@timev solve(prob, Tsit5())

sol = solve(ODEProblem((u, p, t) -> (sleep(0.01); -u), 1.0, nothing),
    Tsit5();
    dt = 0.5,
    tspan = (0.0, 1000.0),
    progress = true,
    progress_steps = 1)

sol.t

