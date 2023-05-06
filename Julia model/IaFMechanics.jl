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
    delta_T = 0.5                           #( ) Sharpness parameter
    theta_rh = 5                            #(V) Rheobase threshold

    Iext::Function = t -> 0                 #(A) External current

    # Reset parameters
    V_F = 5                                 #(V) Firing threshold
    V_R = -10                               #(V) Potential after reset

    # Network parameters
    N::Int; @assert N >= 1                              #(#) Number of neurons
    X::Vector{Tuple{T,T}} = tuple.(rand(N), rand(N))    #(x) Coordinates of neurons (if empty, model is non-spatial)
    wdistr::Symbol = :constant                          # Type of coupling/weight function
    w0 = 1                                              # (Maximum) synaptic weight
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
"""
function genu0(para::IaFParameters)::Vector{Float64}
    @unpack r, u0distr, V_F, N = para

    @assert r >= 0 "Radius must be nonnegative."

    if u0distr == :range                                            # deterministic uniform
        u0 = N > 1 ? range(min(0.99 * V_F, r), -r, N) : [0]

    elseif u0distr == :uniform                                      #~Unif[-r,r]
        d = Uniform(-r, min(0.99 * V_F, r))
        u0 = rand(d, N)

    elseif u0distr == :normal                                       #~N(0,r^2) (censored)
        d = censored(Normal(0,r), upper=0.99 * V_F)
        u0 = rand(d, N)

    else
        error("Invalid type of u0 distribution: $u0distr")
    end

    return u0
end

# TODO: add spatial dependence to Iext
"""
Generate the driving force function f. Takes into account the sparsity pattern of the Jacobian.
"""
function genf(para::IaFParameters)
    @unpack leaky, delta_T, theta_rh, Iext, N = para

    function fleaky!(du, u, p, t)
        @. du = -u + Iext(t)
    end

    # Exponential driving force
    # fexp = @ode_def ExponentialIaF begin
    #     du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R * Iext(t)) / tau
    # end V_rest delta_T theta_rh R tau

    function fexp!(du, u, p, t)
        @. du = -u + delta_T * exp((u - theta_rh) / delta_T) + Iext(t)
    end

    return ODEFunction(leaky ? fleaky! : fexp!, jac_prototype=sparse(I, N, N))
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
    
    firedNeurons::Vector{Int} = [firingNeuron]          # track indices of neurons that have already fired at this moment in time

    while firingNeuron != 0

        # Increment potential of all neurons based on weights
        integrator.u .+= integrator.p.W[:, firingNeuron]      # TODO: only doing this for unfired neurons might be better (or not)

        # The arrival of the spike might have caused new neurons to fire
        # Order of processing spikes shouldn't matter (because we reset at end), so we do a linear search and take the first one

        firingNeuron = 0    # placeholder value for if there are no more firing neurons

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
    push!(integrator.p.spikes, [integrator.t, firedNeurons])
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
"""
function solveiaf(para::IaFParameters, u0::Vector{Float64})
    @unpack leaky, tend, N, spikes = para

    empty!(spikes)      # clear spike data from a possible previous run

    f = genf(para)
    tspan = (0, tend)

    prob = ODEProblem(f, u0, tspan, para)
    sol = solve(prob; reltol=1e-4, abstol=1e-6, callback=gencallback(N), dense=false)

    return sol
end

end
