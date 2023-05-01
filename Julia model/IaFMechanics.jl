"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""
module IaFMechanics

export IaFParameters, genu0, genfleaky, genfexp, gencallback

using Parameters, DataFrames, DifferentialEquations, LinearAlgebra, SparseArrays#, ParameterizedFunctions


"""
Parameters related to the evolution of the potential of a single neuron, the shape of the network, and neuron interactions.
Also contains logging of the spikes, so new instances should be created for every new simulation.
"""
@with_kw struct IaFParameters{T<:Real} @deftype T
    # Driving force parameters
    delta_T = 0.5       #( ) Sharpness parameter
    theta_rh = 5        #(V) Rheobase threshold

    Iext::Function = t -> 0    #(A) External current
    # Iext = t -> exp(t/10)

    # Reset parameters
    V_F = 5             #(V) Firing threshold
    V_R = -10           #(V) Potential after reset

    # Network parameters
    N::Int = 100                            #(#) number of neurons
    w0 = 1
    W::Matrix{T} = (w0 / N) * ones(N, N)    # synaptic weights

    # Logging variables
    spikes::DataFrame = DataFrame(t = T[], cnt = Int[])     # TODO: make cb initialize this? Or recreate this struct everytime
end


"""
Create the initial condition u0 as a uniformly spaced vector of length N.
"""
function genu0(r::Real, parameters::IaFParameters)::Vector{Float64}
    @unpack V_F, N = parameters

    return N > 1 ? range(min(V_F, r), -r, N) : [0]

    # TODO:
    # u0 = 2*r*(rand(N,1)-0.5)   # ~Unif[-r,r]
    # u0 = r*randn(N,1)          # ~N(0,r^2)
end


# Leaky driving force
function genfleaky(p::IaFParameters)
    @unpack Iext, N = p

    function fleaky!(du, u, p, t)
        @. du = -u + Iext(t)
    end

    return ODEFunction(fleaky!, jac_prototype=sparse(I, N, N))
end


# Exponential driving force
# fexp = @ode_def ExponentialIaF begin
#     du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R * Iext(t)) / tau
# end V_rest delta_T theta_rh R tau

function genfexp(p::IaFParameters)
    @unpack delta_T, theta_rh, Iext, N = p

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


# TODO: this logic currently does not warrant its own function
# function solveiaf(f::Function, u0::Vector{Float64}, tStart::Real, tEnd::Real, options::NamedTuple)
#     tspan = (tStart, tEnd)

#     prob = ODEProblem(f, u0, tspan)
#     sol = solve(prob; options...)

#     return sol
# end

end
