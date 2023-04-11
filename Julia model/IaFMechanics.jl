"""
Contains the mechanics needed to numerically solve an integrate-and-fire model.
"""
module IaFMechanics

export DrivingForceParameters, generate_u0, options

using Parameters, DifferentialEquations, ParameterizedFunctions


@with_kw struct DrivingForceParameters{T<:Real} @deftype T
    tau = 1          #(s) Time constant
    V_rest = 0       #(V) Resting potential
    delta_T = 0.5    #( ) Sharpness parameter
    theta_rh = 5     #(V) Rheobase threshold
    R = 1            #(Ω) Resistance
end


## Network settings
# TODO: turn these into settable parameters
const N::Int = 100              # number of neurons
const w0::Float64 = -1
const W = (w0 / N) * ones(N, N)  # synaptic weights
# W = w0*hilb(N) / N


## Initial conditions
"Create the initial condition u0 as a uniformly spaced vector of length N."
function generate_u0(r::Real, V_rest::Real, N::Int)::Vector{Float64}
    return range(V_rest + r, V_rest - r, N)
end

# TODO:
# u0 = 2*r*(rand(N,1)-0.5)   # ~Unif[-r,r]
# u0 = r*randn(N,1)          # ~N(0,r^2)



## Reset conditions

const V_F::Float64 = 5             #(V) Firing threshold
const V_R::Float64 = -10           #(V) Potential after reset

"""
Check whether any neurons have fired. Used for Callback.
"""
function fireCondition!(out, u::Vector{Float64}, t::Real, integrator)
    out .= u .- V_F
end


"""
Perform post-spike potential updates, where firingNeuron is the neuron that fires (first). Used for Callback.
"""
function fire!(integrator, firingNeuron::Int)
    
    firedNeurons::Vector{Int} = [firingNeuron]

    while firingNeuron != 0

        # Increment potential of all neurons based on weights
        integrator.u .+= W[:,firingNeuron]                      # TODO: only doing this for unfired neurons might be better (or not)

        # The arrival of the spike might have caused new neurons to fire
        # Order of processing spikes shouldn't matter (because we reset at end), so we just do a linear search and take the first one
        for i = 1:length(integrator.u)
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
# affect_neg! is nothing, as it is called when fireCondition! is found to be 0
# and the cross is a downcrossing (from positive to negative)
# This is nothing because this should not happen



## ode solver options

options = (reltol = 1e-4,
           abstol = 1e-6,
           callback = cb)

end
