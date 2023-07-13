includet("IaFMechanics.jl")
includet("IaFVisualizations.jl")

using Parameters, DifferentialEquations, Random, Plots, LaTeXStrings, Printf
using .IaFMechanics, .IaFVisualizations


#########################
### Parameter setting ###
#########################

Random.seed!(0)     # for replicability

### Choice of model and performance impacting settings

para = IaFParameters{Float64}(
    leaky = false,
    N = 200,
    tend = 20
    )


### Initial conditions (potentials & weights)

para = IaFParameters(para,
    u0distr = :uniform,
    r = 0.85 - (para.leaky ? 0.37 : 0.17),
    w0distr = :gaussian,
    w0 = 0.1,
    sig1 = 0.2,
    sig2 = 0.5
    )


### Resting potential & external stimulus

# I0 = para.leaky ? 0.4 : 0.3   # constant input sufficient to generate spikes, depends on model
I0 = para.leaky ? 0.65 : 0.2   # constant input sufficient for firing period of approximately 4 seconds

para = IaFParameters(para,
    V_rest = para.leaky ? 0.37 : 0.17,
    # Iext = (t,x) -> I0
    # Iext = (t,x) -> I0*(0.8 + 0.5*x[1])
    Iext = (t,x) -> (0.5*I0 + 0.3*I0*(x[1] < 0.2))
    # Iext = function(t,x)
    #     (
    #     (0.9*I0 + 0.3*I0*(x[1] < 0.2)) * (t < para.tend/2)
    #     +
    #     I0 * (t >= para.tend/2)
    #     )
    # end
    # Iext = (t,x) -> 0.95*I0 + sin(t*π/10)*0.5*I0*(x[1] < 0.3)*(x[2] < 0.3) - sin(t*π/10)*0.5*I0*(x[1] > 0.7)*(x[2] > 0.7),       # sin(t), exp(t/10)
    )


### Learning rule

@unpack w0, N, tend = para

para = IaFParameters(para,
    learning = true,
    # Fp = (w,Δt) -> abs(w)/w0 * (w0 - w) * exp(-2*Δt),
    # Fm = (w,Δt) -> -abs(w)/w0 * (w0 + w) * exp(-4*Δt),
    Fp = (w,Δt) -> 0.1*0.5w0 * (w0 - w) * exp(-2*Δt),
    Fm = (w,Δt) -> 0,
    Gp = w -> -0.1w,
    Gm = w -> 0
    )



# Wbal = false

# if Wbal
#     para = IaFParameters(para, W=(para.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
# end



###############
### Solving ###
###############

# u0 = [3.8, 1.5, 1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]     # somewhat pathological syncing counter-example
# sol = solveiaf(para, u0; savepotentials=true)


if para.learning
    sol, savedweights = solveiaf(para; savepotentials=true, saveweightsperiod=para.tend/2)
else
    sol = solveiaf(para; savepotentials=true)
end

# sol = solveiaf(para, [0.4]; savepotentials=true)
sol = solveiaf(para; savepotentials=true)

println("System solved")


############################
### Plotting & Animating ###
############################

save = false

if N == 1
    pu = plot(sol, idxs=1, continuity=:right, grid=:none, xlabel=L"$t$ (s)", ylabel=L"$u$ (V)")
    display(p)
    png(p, "images/u_" * genfilesuffix(para))
else
    uvaplot(sol, para; activitybinsize=0.2, u=true, v=true, save=save)
end

fps = 10

# udensityanim(sol, para; binsize=50/N, fps=fps, playspeed=2, save=true)
# utorusanim(sol, para; fps=fps, playspeed=1.5, save=save)
# uspatialanim(sol, para; fps=fps, playspeed=0.5, save=save)
Aspatialanim(para; spatialbinsize=0.2, timebinsize=0.5, playspeed=0.75, save=save)

if para.learning
    Wplot(para, savedweights; binsize=0.2, save=save)
end


# Separate activity
# pa = histogram(para.spikes.t, bins=0:0.5:para.tend, weights=1/N*length.(para.spikes.neurons), normalize=:density,
#                     xlabel=L"$t$ (s)", ylabel=L"$A(t)$ (#/s)", label="Numerical")
# plot!([0, tend], [0.454462, 0.454462], label="Theoretical", lw=2)
# plot!(legend=:topright)
# png(pa, "images/a_" * genfilesuffix(para))


# Potential density
# pu = histogram(sol(tend), bins=range(0, 1, step=25/N), normalize=:pdf, xlabel=L"$u$",
#                     xlims=(0,1), ylims=(0,N/250), label="Numerical")

# pexact = u -> 0.454462 / (1.12455 - u)
# plot!(pexact, 0:0.001:1, lw=2, label="Theoretical")
# plot!(legend=:topleft)
# png(pu, "images/udensity_" * genfilesuffix(para))



# Spatial still
# spatplots = []

# for t in [0, para.spikes.t[126],  para.spikes.t[176]]
#     spatplot = scatter(para.X, marker_z=sol(t), color=cgrad([:blue, :orange], [0.5, 0.95]),
#                 xlabel=L"$x^1$", ylabel=L"$x^2$",
#                 cbar=true, cbarlims=(0,1), cbar_title=L"$u$", colorbar_titlefontrotation=270.0,
#                 xlims=(-0.1,1.1), ylims=(-0.1,1.1))#, title="t = $(@sprintf("%.2f", t))")
#     if t in para.spikes.t
#         scatter!(para.X[para.spikes[findfirst(==(t), para.spikes.t), :neurons]], mc=:red, ms=8)
#     end
#     display(spatplot)

#     # png(spatplot, "images/spat;t=$(@sprintf("%.2f", t))_" * genfilesuffix(para))
# end
