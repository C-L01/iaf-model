"""
Contains visualizations for the solution of an integrate-and-fire model.
"""
module IaFVisualizations

# include("IaFMechanics.jl")

export genfilesuffix, uvaplot, udensityanim, utorusanim, uspatialanim, Aspatialanim, Wplot

using Parameters, DifferentialEquations, Plots, Statistics, LaTeXStrings, Printf
# using .IaFMechanics: IaFParameters

# plotlyjs()
gr()
default(legend=false)

"""
Generate file suffix used for saving plots. Assumes Iext is constant.
"""
function genfilesuffix(para)
    @unpack leaky, N, w0, Iext, tend, r, u0distr, w0distr, sig1, sig2, learning = para
    includesigmas = (w0distr != :constant)       # the sigmas are only relevant with non-constant initial weights
    return (
        "$(leaky ? "leaky" : "exp");N=$N;"
        * "I0=$(round(Iext(0,(0,0)), digits=2));tend=$tend"
        * (N > 1 ? 
        (";w0=$(round(w0, digits=2));"
        * "u0distr=$u0distr;r=$r;w0distr=$w0distr;"
        * "$(includesigmas ? "sig1=$sig1;sig2=$sig2;" : "")"
        * "$(para.learning ? "learning" : "")")
        : "")
    )
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

    pa = histogram(para.spikes.t, bins=0:activitybinsize:tend, weights=1/N*length.(para.spikes.neurons), normalize=:density,
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
function udensityanim(sol::ODESolution, para; binsize::Real = 0.5, fps::Int = 5, playspeed::Real = 1, save::Bool = false)
    @unpack tend = para

    timesteps = range(0, tend, round(Int, fps*tend))       # NOTE: tstart = 0
    bins = range(0, 1, step=binsize)
    ymax = 1/binsize    # because of pdf normalization

    # TODO: maybe add vertical line at V_rest
    udensityanim = @animate for t in timesteps
        histogram(sol(t), bins=bins, normalize=:pdf,
                    xlims=(0,1), ylims=(0,ymax), title="Potential density at t = $(@sprintf("%.2f", t))")
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(udensityanim, "animations/udensity_" * filesuffix * ".avi", fps=playspeed*fps)
    end

    display(gif(udensityanim, fps=playspeed*fps))
end


# TODO: maybe replace timebinsize with fps? But then fps should not be too high, else you lose smoothness
"""
Create a spatial animation of the activity evolution. The fps depends on the timebinsize, so it is not a parameter.
"""
function Aspatialanim(para; spatialbinsize::Float64 = 0.1, timebinsize::Float64 = 0.2,
                            playspeed::Real = 1, save::Bool = false)
    @unpack N, tend, X, spikes = para

    timesteps = range(0, tend, step=timebinsize)       # NOTE: tstart = 0
    nrbins = ceil(1/spatialbinsize)^2

    # Determine for each timestep which (if any) spike index occurred after it first
    spikeindicesedges = [findfirst(>=(t), spikes.t) for t in timesteps]

    Aspatialanim = @animate for i = 1:(length(timesteps)-1)
        if isnothing(spikeindicesedges[i])
            spikingneurons = []
        elseif isnothing(spikeindicesedges[i+1])
            spikingneurons = vcat(spikes[spikeindicesedges[i]:end, :neurons]...)
        else
            spikingneurons = vcat(spikes[spikeindicesedges[i]:(spikeindicesedges[i+1]-1), :neurons]...)
        end
        
        spikelocations = map(j -> X[j], spikingneurons)
        histogram2d(spikelocations, bins=0:spatialbinsize:1, color=cgrad(:YlOrRd_9), cbarlims=(0,N/nrbins),
                    title="Activity on [$(timesteps[i]), $(timesteps[i+1])) (s)", xlabel=L"$x_1$", ylabel=L"$x_2$")
        
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(Aspatialanim, "animations/Aspatialanim_" * filesuffix * ".avi", fps=playspeed/timebinsize)
    end

    display(gif(Aspatialanim, fps=playspeed/timebinsize))
end


"""
Create an animation of the evolution of the potentials mapped onto a torus from 0 to 1.
"""
function utorusanim(sol::ODESolution, para; fps::Int = 5, playspeed::Real = 1, save::Bool = false)
    @unpack tend = para

    timesteps = range(0, tend, round(Int, fps*tend))       # NOTE: tstart = 0

    # Potential TODO
    # @userplot TorusPlot
    # @recipe function g(tp::TorusPlot)
    #     u = tp.args
    #     angles = 2π * u + π/2
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
        potentialangles = 2π * sol(t) .+ π/2
        scatter(cos.(potentialangles), sin.(potentialangles), mc=:red, ms=8, ma=0.3, framestyle=:none,
                xlims=(-1.5,1.5), ylims=(-1.1,1.1), title="Potential torus at t = $(@sprintf("%.2f", t))")
        plot!(torusx, torusy, lc=:black)
        scatter!(resetx, resety, shape=:vline, mc=:black, ms=12)
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(torusanim, "animations/utorus_" * filesuffix * ".avi", fps=playspeed*fps)
    end

    display(gif(torusanim, fps=playspeed*fps))
end


"""
Create an animation of the evolution of the potentials visualized per location.
"""
function uspatialanim(sol::ODESolution, para; fps::Int = 10, playspeed::Real = 1, save::Bool = false)
    @unpack tend, X, spikes = para

    timesteps::Vector{Float64} = copy(spikes.t)                             # stop at all spikes
    
    numtargetframes = round(Int, fps*tend)
    if numtargetframes > length(timesteps)
        additionaltimesteps::Vector{Float64} = range(0, tend, numtargetframes - length(timesteps) + 1)
        append!(timesteps, additionaltimesteps)
        sort!(timesteps)
        realfps = fps
    else
        realfps = round(Int, length(timesteps)/tend)
        println("Spikes alone account for all frames, real fps is $realfps.")
    end

    uspatialanim = @animate for t in timesteps
        scatter(X, marker_z=sol(t), color=cgrad(:blues, rev=false), framestyle=:none,
                    xlims=(-0.1,1.1), ylims=(-0.1,1.1), title="Spatial potentials at t = $(@sprintf("%.2f", t))")
        if t in spikes.t
            scatter!(X[spikes[findfirst(==(t), spikes.t), :neurons]], mc=:red, ms=8)
        end
    end

    if save
        filesuffix = genfilesuffix(para)
        gif(uspatialanim, "animations/uspatial_" * filesuffix * ".avi", fps=playspeed*realfps)
    end

    display(gif(uspatialanim, fps=playspeed*realfps))
end


"""
Create a plot of the average synaptic strength between spatial bins of neurons.
"""
function Wplot(para, savedweights::SavedValues; binsize::Float64 = 0.1, save::Bool = false)
    @unpack X, w0, N = para

    preplots = []
    postplots = []

    nrbins = ceil(1/binsize)^2

    for n in eachindex(savedweights.t)

        W = savedweights.saveval[n]

        total_presynaptic = vec(sum(W, dims=2))         # total incoming weight strength per neuron
        total_postsynaptic = vec(sum(W, dims=1))        # total outgoing weight strength per neuron

        push!(preplots, histogram2d(X, bins=0:binsize:1, weights=total_presynaptic, color=cgrad(:roma, rev=true), cbarlims=(-N*w0/nrbins,N*w0/nrbins),
                        title="Pre on $(savedweights.t[n]) (s)", xlabel=L"$x_1$", ylabel=L"$x_2$"))

        push!(postplots, histogram2d(X, bins=0:binsize:1, weights=total_postsynaptic, color=cgrad(:roma, rev=true), cbarlims=(-N*w0/nrbins,N*w0/nrbins),
                        title="Post on $(savedweights.t[n]) (s)", xlabel=L"$x_1$", ylabel=L"$x_2$"))

        # for j in spikingneurons
        #     postneurons = 1:N .!= j
            # quiver!(repeat(X[j,:], N-1), quiver=(X[postneurons,:]),
            #         line_z=W[postneurons, j], color=cgrad(:redgreensplit, rev=false), lw=3)
        #     plot!([([X[j,1], X[i,1]], [X[j,2], X[i,2]]) for i = 1:N if i != j],
        #             color=cw.(W[postneurons, j]))
        # end
    end

    Wplot = plot(preplots..., postplots..., layout=(2, length(preplots)), link=:x, thickness_scaling=1)
    
    if save
        filesuffix = genfilesuffix(para)
        png(Wplot, "images/W_" * filesuffix)
    end

    display(Wplot)
end

end
