includet("IaFMechanics.jl")
includet("IaFVisualizations.jl")

using Parameters, DifferentialEquations, Random#, Plots#, Statistics, LaTeXStrings, Printf
using .IaFMechanics, .IaFVisualizations


### Parameter setting

Random.seed!(0)     # for replicability

para = IaFParameters{Float64}(
    leaky = true,
    N = 100,
    tend = 50,
    wdistr = :constant,
    w0 = 1,
    sig1 = 0.2,
    sig2 = 0.5,
    u0distr = :uniform,
    r = 5
    )

# Update parameters that (usually) depend on other parameters

para = IaFParameters(para,
    Iext = para.leaky ? (t -> 6) : (t -> 4.75),       # sin(t), exp(t/10)
    V_F = para.leaky ? 5 : 7,
    )

# Wbal = false

# if Wbal
#     para = IaFParameters(para, W=(para.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
# end

u0 = genu0(para)      # not wrapped in solveiaf for now, maybe later
# u0 = [3, 2, -5.5, -9]
updateW(para)         # not wrapped in solveiaf for now, maybe later


### Solving

sol = solveiaf(para, u0)

println("System solved")



### Plotting

save = false

uvaplot(sol, para; save=save)


### Animating

fps = 10

udensityanim(sol, para; fps=fps, playspeed=2, save=save)
utorusanim(sol, para; fps=fps, playspeed=1, save=save)
uspatialanim(sol, para; fps=fps, playspeed=1, save=save)

