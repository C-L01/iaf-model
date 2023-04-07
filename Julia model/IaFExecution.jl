using .IaFMechanics, Parameters, DifferentialEquations, ParameterizedFunctions, Plots

# TODO: maybe create such a function in the future
# function solveiaf(u0, tStart, tEnd)#, params::Para)
#     # @unpack f, p, options = params

#     tspan = (tStart, tEnd)
    
#     # Solve until the first spike (or tEnd if there are none)
#     prob = ODEProblem(f, u0, tspan, p)
#     sol = solve(prob; options...)

#     return sol
# end




# @unpack f, p, options = params

tspan = (tStart, tEnd)

# Solve until the first spike (or tEnd if there are none)
prob = ODEProblem(f, u0, tspan, p)
sol = solve(prob; options...)





#(V) Starting potentials
r::Float64 = 5                        #(V) (Expected) radius of the initial potentials
u0::Vector{Float64} = range(-r, r, N)
# u0 = 2*r*(rand(N,1)-0.5)   # ~Unif[-r,r]
# u0 = r*randn(N,1)          # ~N(0,r^2)

sol = solveiaf(u0, 0, 5)#, params)

# plot(sol, idxs = (1,2))
# plot(sol.t, sol.u[1])
