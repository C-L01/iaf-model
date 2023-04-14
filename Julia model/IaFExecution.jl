includet("IaFMechanics.jl")
# includet("IaFVisualization.jl")

using Parameters, DifferentialEquations, Plots, Statistics
using .IaFMechanics, .IaFVisualization

# Get parameters
parameters = IaFParameters{Float64}(Iext=(t -> 6), N=100, w0=1)

u0 = genu0(3, parameters)         # uniformly spaced with some radius

tspan = (0, 20)

prob = ODEProblem(genfLeaky(parameters), u0, tspan, parameters)
sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(parameters.N))

display(plot(sol))
# plot(sol, vars=1, continuity=:right)

plot(sol.t, var.(sol.u))

# png("uplot")

# animate(sol)
