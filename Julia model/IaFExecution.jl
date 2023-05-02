includet("IaFMechanics.jl")
includet("IaFVisualizations.jl")

using Parameters, DifferentialEquations, Random#, Plots, Statistics, LaTeXStrings, Printf
using .IaFMechanics, .IaFVisualizations


### Parameter setting

Random.seed!(0)     # in case of random u0

para = IaFParameters{Float64}(
    leaky = true,
    N = 50,
    w0 = 1,
    tend = 50,
    u0distr = :Uniform,
    r = 10
    )

# Update parameters that (usually) depend on other parameters

para = IaFParameters(para,
    Iext = para.leaky ? (t -> 6) : (t -> 4.75),       # sin(t), exp(t/10)
    V_F = para.leaky ? 5 : 7
    )

# Wbal = false

# if Wbal
#     para = IaFParameters(para, W=(para.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
# end

u0 = genu0(para)      # not wrapped in solveiaf for now, maybe later



### Solving

sol = solveiaf(para, u0)

println("System solved")



### Plotting

uvaplot(sol, para; save=false)


### Animating

fps = 5

udensityanim(sol, para, fps=fps; save=false)
torusanim(sol, para, fps=fps; save=false)
