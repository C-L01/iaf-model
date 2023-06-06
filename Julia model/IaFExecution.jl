includet("IaFMechanics.jl")
includet("IaFVisualizations.jl")

using Parameters, DifferentialEquations, Random#, Plots#, Statistics, LaTeXStrings, Printf
using .IaFMechanics, .IaFVisualizations


### Parameter setting

Random.seed!(0)     # for replicability

para = IaFParameters{Float64}(
    leaky = false,
    N = 300,
    tend = 40,
    u0distr = :uniform,
    r = 0.4,
    wdistr = :gaussian,
    w0 = 5 / 15,
    sig1 = 0.2,
    sig2 = 0.5
    )

# Update parameters that (usually) depend on other parameters

I0 = para.leaky ? 0.4 : 0.3   # constant input sufficient to generate spikes, depends on model

para = IaFParameters(para,
    Iext = (t,x) -> 0.9*I0 + 0.3*I0*(x[1] < 0.2),
    # Iext = (t,x) -> 0.95*I0 + sin(t*π/10)*0.5*I0*(x[1] < 0.3)*(x[2] < 0.3) - sin(t*π/10)*0.5*I0*(x[1] > 0.7)*(x[2] > 0.7),       # sin(t), exp(t/10)
    V_rest = para.leaky ? 10 / 15 : 10 / 17,
    )

# Wbal = false

# if Wbal
#     para = IaFParameters(para, W=(para.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
# end


### Solving

# u0 = [3.8, 1.5, 1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]     # somewhat pathological counter-syncing example
# sol = solveiaf(para, u0; savepotentials=true)

sol = solveiaf(para; savepotentials=true)

println("System solved")



### Plotting

save = false

uvaplot(sol, para; save=save)


### Animating

fps = 10

# udensityanim(sol, para; fps=fps, playspeed=2, save=save)
# utorusanim(sol, para; fps=fps, playspeed=1.5, save=save)
# uspatialanim(sol, para; fps=fps, playspeed=1, save=save)
# Aspatialanim(para; spatialbinsize=0.1, timebinsize=0.4, playspeed=1, save=save)
