includet("IaFMechanics.jl")

using .IaFMechanics, Parameters, DifferentialEquations, Plots

# Get parameters
parameters = IaFParameters{Float64}(Iext=(t -> 6), N=100, w0=1)

u0 = genu0(3, parameters)         # uniformly spaced with some radius
f = genf(parameters; exponential=false)

tspan = (0, 20)

prob = ODEProblem(f, u0, tspan)
sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(parameters))

plot(sol)
# plot(sol, idxs=1)
