includet("IaFMechanics.jl")

using .IaFMechanics, Parameters, DifferentialEquations, ParameterizedFunctions, LinearAlgebra, Plots


#(A) External stimulus
Iext = t -> 6
# Iext(t) = exp(t/10)

# TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
# tPulse = 1                     #(s) time of pulse
# dt = 1e-2                      #(s) time interval of pulse input
# dV = 5.5                         #(V) voltage increase due to pulse
# Iext(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))


# Exponential driving force
fExp = @ode_def ExponentialIaF begin
    du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R*Iext(t)) / tau
end V_rest delta_T theta_rh R tau

# Leaky driving force
function fLeaky(du, u, p, t)
    # du = (-(u - V_rest) + R*Iext(t)) / tau
    @. du = ( -(u - p[1]) + p[2]*Iext(t) ) / p[3]
    # mul!(du, UniformScaling(-1/tau), u)
end


function solveiaf(u0::Vector{Float64}, tStart::Real, tEnd::Real, params::DrivingForceParameters)
    # @unpack f, p, options = params
    @unpack V_rest, R, tau = params

    tspan = (tStart, tEnd)
    
    # Solve until the first spike (or tEnd if there are none)
    prob = ODEProblem(fLeaky, u0, tspan, (V_rest, R, tau))
    sol = solve(prob; options...)

    return sol
end

N::Int = IaFMechanics.N

params = DrivingForceParameters{Float64}()
u0 = generate_u0(3, params.V_rest, 100)         # uniformly spaced around 0

sol = solveiaf(u0, 0, 20, params)

plot(sol)
# plot(sol.t, sol.u[1])
