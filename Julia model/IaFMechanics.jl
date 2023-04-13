"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""
module IaFMechanics

export IaFParameters, genu0, genf, gencallback

using Parameters, DifferentialEquations#, ParameterizedFunctions


"""
Parameters related to the evolution of the potential of a single neuron, the shape of the network, and neuron interactions.
"""
@with_kw struct IaFParameters{T<:Real} @deftype T
    # TODO: maybe wrap in @consts begin end

    # Driving force parameters
    tau = 1             #(s) Time constant
    V_rest = 0          #(V) Resting potential
    delta_T = 0.5       #( ) Sharpness parameter
    theta_rh = 5        #(V) Rheobase threshold
    R = 1               #(Ω) Resistance

    Iext::Function = t -> 0    #(A) External current
    # Iext = t -> exp(t/10)

    # TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
    # tPulse = 1                     #(s) time of pulse
    # dt = 1e-2                      #(s) time interval of pulse input
    # dV = 5.5                         #(V) voltage increase due to pulse
    # Iext(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))

    # Reset parameters
    V_F = 5             #(V) Firing threshold
    V_R = -10           #(V) Potential after reset

    # Network parameters
    N::Int = 100                            #(#) number of neurons
    w0 = 1
    W::Matrix{T} = (w0 / N) * ones(N, N)    # synaptic weights
end


"""
Create the initial condition u0 as a uniformly spaced vector of length N.
"""
function genu0(r::Real, parameters::IaFParameters)::Vector{Float64}
    @unpack V_rest, N = parameters

    return range(V_rest + r, V_rest - r, N)

    # TODO:
    # u0 = 2*r*(rand(N,1)-0.5)   # ~Unif[-r,r]
    # u0 = r*randn(N,1)          # ~N(0,r^2)
end


"""
Create the driving force function f based on given parameters.
"""
function genf(parameters::IaFParameters; exponential::Bool = false)::Function
    @unpack tau, V_rest, delta_T, theta_rh, R, Iext = parameters

    # Exponential driving force
    if exponential
        # fExp = @ode_def ExponentialIaF begin
        #     du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R*Iext(t)) / tau
        # end V_rest delta_T theta_rh R tau

        return function fExp(du, u, p, t)
            @. du = ( -(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R*Iext(t) ) / tau
        end

    # Leaky driving force
    else
        return function fLeaky(du, u, p, t)
            @. du = ( -(u - V_rest) + R*Iext(t) ) / tau
            # mul!(du, UniformScaling(-1/tau), u)
        end
    end
end


"""
Create the VectorContinuousCallback used by the ODE solver to check for and handle spikes.
"""
function gencallback(parameters::IaFParameters)::VectorContinuousCallback
    @unpack V_R, V_F, N, W = parameters

    """
    Check whether any neurons have fired. Used for Callback.
    """
    function fireCondition!(out, u::Vector{Float64}, t::Real, integrator)
        @. out = u - V_F
    end


    """
    Perform post-spike potential updates, where firingNeuron is the neuron that fires (first). Used for Callback.
    """
    function fire!(integrator, firingNeuron::Int)
        
        firedNeurons::Vector{Int} = [firingNeuron]

        while firingNeuron != 0

            # Increment potential of all neurons based on weights
            integrator.u .+= W[:,firingNeuron]          # TODO: only doing this for unfired neurons might be better (or not)

            # The arrival of the spike might have caused new neurons to fire
            # Order of processing spikes shouldn't matter (because we reset at end), so we do a linear search and take the first one

            firingNeuron = 0    # keeps this value if there are no more firing neurons

            for i = 1:N
                if !(i in firedNeurons) && integrator.u[i] >= V_F
                    firingNeuron = i
                    push!(firedNeurons, firingNeuron)
                    break
                end
            end
        end

        # Reset potentials of fired neurons
        integrator.u[firedNeurons] .= V_R
    end

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