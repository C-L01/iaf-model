"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""
module IaFMechanics

export IaFParameters, solveiaf, genu0, updateW

using Parameters, DataFrames, Distributions, DifferentialEquations, LinearAlgebra, SparseArrays#, ParameterizedFunctions


"""
Parameters related to the evolution of the potential of a single neuron, the shape of the network, and neuron interactions.
Also contains logging of the spikes, so new instances should be created for every new simulation.
"""
@with_kw struct IaFParameters{T<:Real} @deftype T
    # Driving force type
    leaky::Bool = true                      # Whether to use leaky IaF (alternative is exponential)

    # Driving force parameters (only relevant for exponential)
    delta_T = 0.5 / (7 - (-10))             #( ) Sharpness parameter
    theta_rh = (5 - (-10)) / (7 - (-10))    #(V) Rheobase threshold

    Iext::Function = (t,x) -> 0             #(A) External current

    V_rest = 2/3                            #(V) Normalized resting potential

    # Network parameters
    N::Int; @assert N >= 1                              #(#) Number of neurons
    X::Vector{Tuple{T,T}} = tuple.(rand(N), rand(N))    #(x) Coordinates of neurons (if empty, model is non-spatial)
    wdistr::Symbol = :constant                          # Type of coupling/weight function
    w0 = 1 / (5-(-10))                                  # (Maximum) synaptic weight
    sig1 = 0.2                                          # Standard deviation for Gaussian, and first parameter for Mexican-hat
    sig2 = 0.5                                          # Second Mexican-hat parameter
    W::Matrix{T} = zeros(N,N)                           # Synaptic weights

    # Initial conditions
    r = 3; @assert r >= 0                   #(V) Potential radius around 0
    u0distr::Symbol = :range                # Type of distribution for u0

    # Termination
    tend = 30                               #(s) End time of simulation (NOTE: tstart is normalized to 0)

    # Logging variables
    spikes::DataFrame = DataFrame(t = T[], neurons = Vector{Int}[])
end


"""
Update the weight matrix W, i.e. matrix of w_ij's. Note that in all cases the values on the diagonal do not matter.
Three options for weight distribution:
    - constant
    - gaussian
    - mexican
"""
function updateW(para::IaFParameters)
    @unpack W, w0, wdistr, sig1, sig2, X, N = para

    if wdistr == :constant                                         # standard full-connectivity
        w = x -> w0

    elseif wdistr == :gaussian                                     # Gaussian based on distance
        w = x -> w0/sqrt(2π * sig1) * exp(-x^2/(2sig1^2))        
        
    elseif wdistr == :mexican                                      # "Mexican-hat" based on distance
        @assert sig2 > sig1
        w = x -> w0 * (sig2 * exp(-x^2/(2sig1^2)) - sig1 * exp(-x^2/(2sig2^2))) / (sig2 - sig1)

    else
        error("Invalid type of weight distribution: $wdistr")
    end

    for i = 1:N
        for j = (i+1):N
            W[i,j] = 1/N * w(norm(X[i] .- X[j]))
            W[j,i] = W[i,j]     # because of symmetry of the metric
        end
    end
end



"""
Create the initial condition u0, i.e. a potential vector of length N.
Three options for u0 distribution:
    - range
    - uniform
    - normal
Always returns a sorted u0, from high to low.
"""
function genu0(para::IaFParameters)::Vector{Float64}
    @unpack r, V_rest, u0distr, N = para

    @assert r >= 0 "Radius must be nonnegative."

    if u0distr == :range                                            # deterministic uniform
        u0 = N > 1 ? range(min(0.99, V_rest+r), V_rest-r, N) : [0]

    elseif u0distr == :uniform                                      #~Unif[-r,r]
        d = Uniform(V_rest-r, min(0.99, V_rest+r))
        u0 = sort(rand(d, N), rev=true)

    elseif u0distr == :normal                                       #~N(0,r^2) (censored)
        d = truncated(Normal(V_rest,r), upper=0.99)
        u0 = sort(rand(d, N), rev=true)

    else
        error("Invalid type of u0 distribution: $u0distr")
    end

    return u0
end


"""
Generate the driving force function f. Takes into account the sparsity pattern of the Jacobian.
"""
function genf(para::IaFParameters)
    @unpack leaky, V_rest, delta_T, theta_rh, Iext, N, X = para

    Iextperneuron = t -> map(Ix -> Ix(t), [s -> Iext(s,x) for x in X])

    # For some reason the fact that Iextperneuron is now a vector breaks usage of @.

    function fleaky!(du, u, p, t)
        du .= -(u .- V_rest) .+ Iextperneuron(t)
    end

    # Exponential driving force
    # fexp = @ode_def ExponentialIaF begin
    #     du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R * Iextperneuron(t)) / tau
    # end V_rest delta_T theta_rh R tau

    function fexp!(du, u, p, t)
        du .= -(u .- V_rest) .+ delta_T * exp.((u .- theta_rh) / delta_T) .+ Iextperneuron(t)
    end

    return ODEFunction(leaky ? fleaky! : fexp!, jac_prototype=sparse(I, N, N))
end



"""
Check whether any neurons have fired. Used for Callback.
"""
function fireCondition!(out, u::Vector{Float64}, t::Real, integrator)
    @. out = u - 1
end


"""
Perform post-spike potential updates, where firingneuron is the neuron that fires (first). Used for Callback.

Assumes that the firingneuron it is passed has the highest potential of all neurons (could be tied though).
Assumes that weights are sufficiently small such that potentials will not be elevated all the way from 0 to 1 by other firings alone.
"""
function fire!(integrator, firingneuron::Int)
    # @unpack spikes = integrator.p
    
    # Search for neurons that have exactly the same potential as the firingneuron (they should also fire)
    firingneurons::Vector{Int} = findall(u -> u == integrator.u[firingneuron], integrator.u)

    firedneurons::Vector{Int} = []    # to track neurons that have already fired (at this moment in time)

    while !isempty(firingneurons)

        # Print if there are more then one firingneurons, as this is exceptional
        # TODO: if this never happens, maybe optimize by removing the findalls
        if length(firingneurons) > 1
            println("Neurons $firingneurons fired *exactly* at the same fire iteration")
        end

        # Increment potential of all neurons based on weights
        integrator.u .+= integrator.p.W[:, firingneurons]       # TODO: could exclude firingneurons here because of reset

        # Reset potential of the firing neurons (NOTE: potential(s) can be changed later by new firing neurons)
        integrator.u[firingneurons] .= integrator.p.V_R
        #                         norm(integrator.p.X[j] .- integrator.p.X[i]))
        #                     for sigma in spikes.t
        #                         if j in spikes[findfirst(==(sigma), spikes.t), :neurons])
        #         end
        #     end
        # end


        # Log the neurons that fired in this iteration
        append!(firedneurons, firingneurons)

        # Emtpy firingneurons (to exit while loop eventually)
        empty!(firingneurons)



        ### The arrival of the spike(s) might have caused new neurons to fire

        # Find largest potential
        maxpotential = maximum(integrator.u)

        if maxpotential >= 1
            # Find all neurons with that same potential, which are the new firingneurons
            firingneurons = findall(u -> u == maxpotential, integrator.u)
            @assert !any(i in firingneurons for i in firedneurons) "A neuron refired instantaneously, weights too strong"
        end

        # sortedNeurons = sortperm(integrator.u; alg=QuickSort, rev=true)
    end

    # Log all fired neurons
    push!(integrator.p.spikes, [integrator.t, firedneurons])
end


# TODO: external pulses could be implemented as an integration event; make sure to tstop at the pulse time
# tPulse = 1                     #(s) time of pulse
# dt = 1e-2                      #(s) time interval of pulse input
# dV = 5.5                         #(V) voltage increase due to pulse
# Iext(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))

"""
Create the VectorContinuousCallback used by the ODE solver to check for and handle spikes.
"""
function gencallback(N::Int)::VectorContinuousCallback
    # fireCondition! should never detect a downcrossing root, hence the nothing argument
    return VectorContinuousCallback(fireCondition!, fire!, nothing, N)
end



"""
Wrapper for generating all input for the solver and subsequently solving the integrate-and-fire model for given parameters.
If savepotentials is false, the ODE solver does not save any potential values. Spikes are still saved in para.
"""
function solveiaf(para::IaFParameters, u0::Vector{Float64}; savepotentials::Bool = true)
    @unpack leaky, tend, N, spikes = para

    empty!(spikes)      # clear spike data from a possible previous run

    f = genf(para)
    tspan = (0, tend)

    prob = ODEProblem(f, u0, tspan, para)
    sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(N), save_everystep=savepotentials)

    return sol
end

end
