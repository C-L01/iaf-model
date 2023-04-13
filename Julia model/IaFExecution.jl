includet("IaFMechanics.jl")

using .IaFMechanics, Parameters, DifferentialEquations, Plots


#(A) External stimulus
Iext = t -> 6
# Iext(t) = exp(t/10)

# TODO: external pulses should be implemented as an integration event; make sure to tstop at the pulse time
# tPulse = 1                     #(s) time of pulse
# dt = 1e-2                      #(s) time interval of pulse input
# dV = 5.5                         #(V) voltage increase due to pulse
# Iext(t) = (dV/dt) * (heaviside(t-tPulse) - heaviside(t-tPulse-dt))


# Get parameters
parameters = IaFParameters{Float64}()

## ode solver options
options = (reltol = 1e-4,
           abstol = 1e-6,
           callback = gencallback(parameters))


u0 = genu0(3, parameters)         # uniformly spaced with some radius
f = genf(Iext, parameters)

sol = solveiaf(f, u0, 0, 20, options)

plot(sol)
# plot(sol, idxs=1)
