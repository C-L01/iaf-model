include("IaFMechanics.jl")
using .IaFMechanics, Parameters, DifferentialEquations, ParameterizedFunctions, Plots


#(A) External stimulus
I = t -> 5
# I(t) = exp(t/10)

# TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
# tPulse = 1                     #(s) time of pulse
# dt = 1e-2                      #(s) time interval of pulse input
# dV = 5.5                         #(V) voltage increase due to pulse
# I(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))


# Exponential driving force
fExp = @ode_def ExponentialIaF begin
    du = (-(u - V_rest) + delta_T * exp((u - theta_rh) / delta_T) + R*I(t)) / tau
end V_rest delta_T theta_rh R tau

# Leaky driving force
fLeaky = @ode_def LeakyIaF begin
    du = (-(u - V_rest) + R*I(t)) / tau
end V_rest R tau


function solveiaf(u0, tStart, tEnd, params::DrivingForceParameters)
    # @unpack f, p, options = params
    @unpack V_rest, R, tau = params

    tspan = (tStart, tEnd)
    
    # Solve until the first spike (or tEnd if there are none)
    prob = ODEProblem(fLeaky, u0, tspan, (V_rest, R, tau))
    sol = solve(prob; options...)

    return sol
end

params = DrivingForceParameters{Float64}()
u0 = generate_u0(3)         # uniformly spaced around 0

sol = solveiaf(u0, 0, 5, params)

plot(sol, idxs = [1])
# plot(sol.t, sol.u[1])
