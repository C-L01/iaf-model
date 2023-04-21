includet("IaFMechanics.jl")

using Parameters, DifferentialEquations, Plots, Statistics, LaTeXStrings, Printf
using .IaFMechanics

### Solving

Wbal = false
leaky = true

# Set parameters
Iext = leaky ? (t -> 6) : (t -> 4.75)
V_F = leaky ? 5 : 7

r = 3                         # TODO: maybe put in params struct, especially if solveiaf is revived
tstart, tend = 0, 30

params = IaFParameters{Float64}(Iext=Iext, V_F=V_F, N=5000, w0=-1)

@unpack N = params



if Wbal
    params = IaFParameters(params, W=(params.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
end

u0 = genu0(r, params)         # uniformly spaced with some radius

f = leaky ? genfleaky(params) : genfexp(params)

prob = ODEProblem(f, u0, (tstart,tend), params)
sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(N), dense=false)

println("System solved")


### Plotting

# plotlyjs()
gr()

pu = plot(sol, idxs=[round(Int, i) for i in range(1, N, length=3)],
            title="Potential evolution of a subset of 3 neurons", xlabel="", ylabel="Potential (V)")

# plot(sol, idxs=1, continuity=:right)

pv = plot(sol.t, var.(sol.u), title="Variance", ylabel=L"Var$(t)$ (V$^2$)")

pa = histogram(params.spikes.t, bins=tstart:0.2:tend, weights=1/N*params.spikes.cnt, normalize=:density,
                title="Activity", xlabel=L"$t$ (s)", ylabel=L"$A(t)$ (#/s)")

# Combine plots of potential (u), variance (v) and activity (a)
uvaplot = plot(pu, pv, pa, layout=(3, 1), link=:x, legend=false, thickness_scaling=1)


### Animating

fps = 15
bins = range(params.V_R, params.V_F, step=1)
# bins = :rice

udensityanim = @animate for t in range(tstart, tend, fps*(tend-tstart))
    histogram(sol(t), bins=bins, normalize=:pdf,
                xlims=(params.V_R, params.V_F), ylims=(0,1), title="Potential density at t = $(@sprintf("%.2f", t))", legend=false)          
end


### Saving
     
# Assumes Iext is constant
filesuffix = "$(leaky ? "leaky" : "exp");N=$(params.N);w0=$(params.w0);$(Wbal ? "Wbal;" : "")I0=$(params.Iext(0));r=$r;tend=$tend"

png(uvaplot, "images/uva_" * filesuffix)

gif(udensityanim, "animations/udensity_" * filesuffix * ".avi", fps=fps)

display(uvaplot)
gif(udensityanim, fps=fps)
