module IaFMechanics
"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""

export DrivingForceParameters

using Parameters, DifferentialEquations, ParameterizedFunctions

@with_kw struct DrivingForceParameters{R<:Real} @deftype R
    tau = 1          #(s) Time constant
    V_rest = 0       #(V) Resting potential
    delta_T = 0.5    #( ) Sharpness parameter
    theta_rh = 5     #(V) Rheobase threshold
    R = 1            #(Ω) Resistance
end

#(A) External stimulus
I = t -> 4.75
# I(t) = exp(t/10)

# TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
# tPulse = 1                     #(s) time of pulse
# dt = 1e-2                      #(s) time interval of pulse input
# dV = 5.5                         #(V) voltage increase due to pulse
# I(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))


# Exponential driving force
fExp = @ode_def ExponentialIaF begin
    du = (-(u - V_rest) + delta_T * exp((u-theta_rh) / delta_T) + R*I(t))/tau
end V_rest delta_T theta_rh R tau I

# Leaky driving force
fLeaky = @ode_def ExponentialIaF begin
    du = (-(u - V_rest) + R*I(t))/tau
end V_rest R tau I


## Network settings

N::Int = 1000                  # number of neurons
w0::Float64 = 1
W = (w0 / N) * ones(N, N)        # synaptic weights
# W = w0*hilb(N) / N            # synaptic weights


## Reset conditions

V_F::Float64 = 6             #(V) Firing threshold
V_R::Float64 = -10           #(V) Potential after reset

function fireCondition!(out, u, t, integrator)
    out = u .- V_F
end

function fire!(integrator, firingNeuron::Int)
    """
    Perform post-spike potential updates, where firingNeuron is the neuron that fires (first).
    """
    firedNeurons::Vector{Int} = [firingNeuron]

    while firingNeuron != 0

        # Increment potential of all neurons based on weights
        integrator.u .+= W[:,firingNeuron]                      # TODO: only doing this for unfired neurons might be better (or not)

        # The arrival of the spike might have caused new neurons to fire
        # Order of processing spikes shouldn't matter (because we reset at end), so we just do a linear search and take the first one
        for i = 1:N
            if !(i in firedNeurons) && integrator.u[i] >= V_F
                firingNeuron = i
                push!(firedNeurons, firingNeuron)
                break
            end
            firingNeuron = 0    # no more firing neurons
        end
    end

    # Reset potentials of fired neurons
    integrator.u[firedNeurons] .= V_R
end

cb = VectorContinuousCallback(fireCondition!, fire!, nothing, N)


## ode solver options

options = (reltol = 1e-4,
           abstol = 1e-6,
           callback = cb)


end
