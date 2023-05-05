"""
Contains visualizations for the solution of an integrate-and-fire model.
"""
module IaFVisualizations

# include("IaFMechanics.jl")

export uvaplot, udensityanim, utorusanim, uspatialanim

using Parameters, DifferentialEquations, Plots, Statistics, LaTeXStrings, Printf
# using .IaFMechanics: IaFParameters

# plotlyjs()
gr()
default(legend=false)

"""
Generate file suffix used for saving plots. Assumes Iext is constant.
"""
function genfilesuffix(para)
    @unpack leaky, N, w0, Iext, tend, r, u0distr, wdistr = para
    return "$(leaky ? "leaky" : "exp");N=$N;w0=$w0;I0=$(Iext(0));u0distr=$u0distr;r=$r;wdistr=$wdistr;tend=$tend"
    # $(Wbal ? "Wbal;" : "")
end


"""
Create a stacked plot of 3 potentials, the variance and the activity over time
"""
function uvaplot(sol::ODESolution, para; activitybinsize::Float64 = 0.2, save::Bool = false)
    @unpack N, tend = para

    pu = plot(sol, idxs=[round(Int, i) for i in range(1, N, length=3)], continuity=:right,
                title="Potential evolution of a subset of 3 neurons", xlabel="", ylabel="Potential (V)")

    pv = plot(sol.t, var.(sol.u), title="Variance", ylabel=L"Var$(t)$ (V$^2$)")

    pa = histogram(para.spikes.t, bins=0:activitybinsize:tend, weights=1/N*para.spikes.cnt, normalize=:density,
                    title="Activity", xlabel=L"$t$ (s)", ylabel=L"$A(t)$ (#/s)")

    # Stack the plots
    uvaplot = plot(pu, pv, pa, layout=(3, 1), link=:x, thickness_scaling=1)

    if save
        filesuffix = genfilesuffix(para)
        png(uvaplot, "images/uva_" * filesuffix)
    end

    display(uvaplot)
end


"""
Create an animation of the potential density evolution.
"""
function udensityanim(sol::ODESolution, para; fps::Int = 5, binsize::Float64 = 0.5, save::Bool = false)
    @unpack tend, V_R, V_F = para

    timesteps = range(0, tend, round(Int, fps*tend))       # NOTE: tstart = 0
    bins = range(V_R, V_F, step=binsize)
    # bins = :rice
    ymax = 1/binsize    # because of pdf normalization

    udensityanim = @animate for t in timesteps
        histogram(sol(t), bins=bins, normalize=:pdf,
                    xlims=(V_R,V_F), ylims=(0,ymax), title="Potential density at t = $(@sprintf("%.2f", t))")
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(udensityanim, "animations/udensity_" * filesuffix * ".avi", fps=fps)
    end

    display(gif(udensityanim, fps=fps))
end


"""
Create an animation of the evolution of the potentials mapped onto a torus from V_R to V_F.
"""
function utorusanim(sol::ODESolution, para; fps::Int = 5, save::Bool = false)
    @unpack V_R, V_F, tend = para

    timesteps = range(0, tend, round(Int, fps*tend))       # NOTE: tstart = 0

    # Potential TODO
    # @userplot TorusPlot
    # @recipe function g(tp::TorusPlot)
    #     u = tp.args
    #     angles = 2π * (u .- V_R) / (V_F - V_R) + π/2
    #     cos.(angles), sin.(angles)
    # end

    # Background torus
    τ = range(0, 2π, length=100)
    torusx = cos.(τ)
    torusy = sin.(τ)

    # Reset location
    resetx, resety = [0], [1]

    torusanim = @animate for t in timesteps
        # torusplot(sol(t))
        potentialangles = 2π * (sol(t) .- V_R) / (V_F - V_R) .+ π/2
        scatter(cos.(potentialangles), sin.(potentialangles), mc=:red, ms=8, ma=0.3, framestyle=:none,
                xlims=(-1.5,1.5), ylims=(-1.1,1.1), title="Potential torus at t = $(@sprintf("%.2f", t))")
        plot!(torusx, torusy, lc=:black)
        scatter!(resetx, resety, shape=:vline, mc=:black, ms=12)
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(torusanim, "animations/utorus_" * filesuffix * ".avi", fps=fps/4)
    end

    display(gif(torusanim, fps=fps))
end


"""
Create an animation of the evolution of the potentials visualized per location.
"""
function uspatialanim(sol::ODESolution, para; fps::Int = 10, save::Bool = false)
    @unpack tend, X, spikes = para

    timesteps::Vector{Float64} = range(0, tend, round(Int, fps*tend))       # NOTE: tstart = 0
    append!(timesteps, spikes.t)                                            # add spike times to ensure those are not stepped over
    sort!(timesteps)

    uspatialanim = @animate for t in timesteps
        scatter(X, marker_z=sol(t), color=cgrad(:blues, rev=false),
                    xlims=(-0.1,1.1), ylims=(-0.1,1.1), title="Spatial potentials at t = $(@sprintf("%.2f", t))")
        if t in spikes.t
            scatter!(X[spikes[findfirst(==(t), spikes.t), :neurons]], mc=:red, ms=8)
        end
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(uspatialanim, "animations/uspatial_" * filesuffix * ".avi", fps=fps)
    end

    display(gif(uspatialanim, fps=fps))
end

end
