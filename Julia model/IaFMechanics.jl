"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""
module IaFMechanics

export IaFParameters, solveiaf, genu0

using Parameters, DataFrames, Distributions, DifferentialEquations, LinearAlgebra, SparseArrays#, ParameterizedFunctions


"""
Parameters related to the evolution of the potential of a single neuron, the shape of the network, and neuron interactions.
Also contains logging of the spikes, so new instances should be created for every new simulation.
"""
@with_kw struct IaFParameters{T<:Real} @deftype T
    # Driving force type
    leaky::Bool = true                      # Whether to use leaky IaF (alternative is exponential)

    # Driving force parameters
    delta_T = 0.5                           #( ) Sharpness parameter
    theta_rh = 5                            #(V) Rheobase threshold

    Iext::Function = t -> 0                 #(A) External current

    # Reset parameters
    V_F = 5                                 #(V) Firing threshold
    V_R = -10                               #(V) Potential after reset

    # Network parameters
    N::Int = 100                            #(#) Number of neurons
    w0 = 1
    W::Matrix{T} = (w0 / N) * ones(N, N)    # Synaptic weights

    # Initial conditions
    r = 3                                   #(V) Potential radius around 0
    u0distr::Symbol = :range                # Type of distribution for u0

    # Termination
    tend = 30                               #(s) End time of simulation (NOTE: tstart is normalized to 0)

    # Logging variables
    spikes::DataFrame = DataFrame(t = T[], cnt = Int[])
end


"""
Create the initial condition u0 as a uniformly spaced vector of length N.
"""
function genu0(para::IaFParameters)::Vector{Float64}
    @unpack r, u0distr, V_F, N = para

    @assert r >= 0

    if u0distr == :range                                            # deterministic uniform
        return N > 1 ? range(min((N-1)/N * V_F, r), -r, N) : [0]
    elseif u0distr == :Uniform                                      #~Unif[-r,r]
        d = Uniform(-r,r)
        return rand(d, N)
    elseif u0distr == :Normal                                       #~N(0,r^2) (censored)
        d = censored(Normal(0,r), upper=(N-1)/N * V_F)
        return rand(d, N)
    else
        error("Invalid type of u0 distribution: $u0distr")
    end
end


# Leaky driving force
function genfleaky(para::IaFParameters)
    @unpack Iext, N = para

    function fleaky!(du, u, p, t)
        @. du = -u + Iext(t)
    end

    return ODEFunction(fleaky!, jac_prototype=sparse(I, N, N))
end


# Exponential driving force
# fexp = @ode_def ExponentialIaF begin
#     du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R * Iext(t)) / tau
# end V_rest delta_T theta_rh R tau

function genfexp(para::IaFParameters)
    @unpack delta_T, theta_rh, Iext, N = para

    function fexp!(du, u, p, t)
        @. du = -u + delta_T * exp((u - theta_rh) / delta_T) + Iext(t)
    end

    return ODEFunction(fexp!, jac_prototype=sparse(I, N, N))
end



"""
Check whether any neurons have fired. Used for Callback.
"""
function fireCondition!(out, u::Vector{Float64}, t::Real, integrator)
    @. out = u - integrator.p.V_F
end


"""
Perform post-spike potential updates, where firingNeuron is the neuron that fires (first). Used for Callback.
"""
function fire!(integrator, firingNeuron::Int)
    
    firedNeurons::Vector{Int} = [firingNeuron]

    while firingNeuron != 0

        # Increment potential of all neurons based on weights
        integrator.u .+= integrator.p.W[:, firingNeuron]      # TODO: only doing this for unfired neurons might be better (or not)

        # The arrival of the spike might have caused new neurons to fire
        # Order of processing spikes shouldn't matter (because we reset at end), so we do a linear search and take the first one

        firingNeuron = 0    # keeps this value if there are no more firing neurons

        for i = 1:integrator.p.N
            if integrator.u[i] >= integrator.p.V_F && !(i in firedNeurons)
                firingNeuron = i
                push!(firedNeurons, firingNeuron)
                break
            end
        end
    end

    # Reset potentials of fired neurons
    integrator.u[firedNeurons] .= integrator.p.V_R

    # Log number of fired neurons
    push!(integrator.p.spikes, [integrator.t, length(firedNeurons)])
end


# TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
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
"""
function solveiaf(para::IaFParameters, u0::Vector{Float64})
    @unpack leaky, tend, N, spikes = para

    empty!(spikes)      # clear data from a possible previous run

    f = leaky ? genfleaky(para) : genfexp(para)
    tspan = (0, tend)

    prob = ODEProblem(f, u0, tspan, para)
    sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(N), dense=false)

    return sol
end

end
