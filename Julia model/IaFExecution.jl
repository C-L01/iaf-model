includet("IaFMechanics.jl")

using Parameters, DifferentialEquations, Plots, Statistics, LaTeXStrings
using .IaFMechanics

### Solving

# Get parameters
params = IaFParameters{Float64}(Iext=(t -> 6), N=100, w0=-1)
@unpack N = params

u0 = genu0(3, params)         # uniformly spaced with some radius

tstart, tend = 0, 40

prob = ODEProblem(genfleaky(params), u0, (tstart,tend), params)
sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(N), dense=false)

println("System solved")


### Plotting

# plotlyjs()
gr()

pu = plot(sol, idxs=[round(Int, i) for i in range(1, N, length=3)],
            title="Potential evolution of a subset of 3 neurons", xlabel="", ylabel="Potential (V)")

# plot(sol, idxs=1, continuity=:right)

pv = plot(sol.t, var.(sol.u), title="Variance", ylabel=L"Var$(t)$ (V$^2$)")

pa = histogram(params.spikes.t, bins=tstart:0.1:tend, weights=1/N * params.spikes.cnt, normalize=:density,
                title="Activity", xlabel=L"$t$ (s)", ylabel=L"$A(t)$ (#/s)")

plot(pu, pv, pa, layout=(3, 1), plot_title="Metrics for $N coupled neurons", link=:x, legend=false, thickness_scaling=1)

# png("uplot")

numframes = 500
numbins::Int = (params.V_F - params.V_R)

densityanim = @animate for t in range(tstart, tend, numframes)
    histogram(sol(t), bins=numbins, normalize=:density,         # TODO: find out if there is a more appropriate normalization, maybe scale by N
                xlims=(params.V_R, params.V_F), ylims=(0,40), label="t = $(round(t, digits=2))", legend=:topleft)
end

gif(densityanim, fps=numframes/(tend-tstart))
# gif(densityanim, "densityanim.avi", fps=numframes/(tend-tstart))
