includet("IaFMechanics.jl")
includet("IaFVisualizations.jl")

using Parameters, DifferentialEquations, Random#, Plots#, Statistics, LaTeXStrings, Printf
using .IaFMechanics, .IaFVisualizations


#########################
### Parameter setting ###
#########################

Random.seed!(0)     # for replicability

### Choice of model and performance impacting settings

para = IaFParameters{Float64}(
    leaky = false,
    N = 200,
    tend = 20
    )


### Initial conditions (potentials & weights)

para = IaFParameters(para,
    u0distr = :uniform,
    r = 0.4,
    w0distr = :gaussian,
    w0 = 5 / 15,
    sig1 = 0.2,
    sig2 = 0.5
    )


### Resting potential & external stimulus

I0 = para.leaky ? 0.4 : 0.3   # constant input sufficient to generate spikes, depends on model
# I0 = 0.5

para = IaFParameters(para,
    V_rest = para.leaky ? 10 / 15 : 10 / 17,
    Iext = function(t,x)
        (
        (0.9*I0 + 0.3*I0*(x[1] < 0.2)) * (t < para.tend/2)
        +
        I0 * (t >= para.tend/2)
        )
    end
    # Iext = (t,x) -> 0.95*I0 + sin(t*π/10)*0.5*I0*(x[1] < 0.3)*(x[2] < 0.3) - sin(t*π/10)*0.5*I0*(x[1] > 0.7)*(x[2] > 0.7),       # sin(t), exp(t/10)
    )


### Learning rule

@unpack w0, N = para

para = IaFParameters(para,
    learning = true,
    # TODO: does this inequality asymmetry make sense? Especially Δt == 0 requires attention
    F = function(w,Δt,d)
        (
        abs(w)/(2*w0/N) * (2*w0/N - w) * exp(-2*Δt) * (Δt >= 0)
        -
        abs(w)/(2*w0/N) * (2*w0/N + w) * exp(4*Δt) * (Δt < 0)
        )
    end
    )



# Wbal = false

# if Wbal
#     para = IaFParameters(para, W=(para.w0 / N) * hcat(-ones(N, N÷2), ones(N, N÷2)))   # TODO: does not work for uneven N
# end



###############
### Solving ###
###############

# u0 = [3.8, 1.5, 1, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]     # somewhat pathological counter-syncing example
# sol = solveiaf(para, u0; savepotentials=true)


if para.learning
    sol, savedweights = solveiaf(para; savepotentials=true, saveweightsperiod=para.tend/2)
else
    sol = solveiaf(para; savepotentials=true)
end

println("System solved")


############################
### Plotting & Animating ###
############################

save = false

uvaplot(sol, para; save=save)


fps = 10

# udensityanim(sol, para; fps=fps, playspeed=2, save=save)
# utorusanim(sol, para; fps=fps, playspeed=1.5, save=save)
# uspatialanim(sol, para; fps=fps, playspeed=1, save=save)
Aspatialanim(para; spatialbinsize=0.2, timebinsize=0.2, playspeed=0.5, save=save)

if para.learning
    Wplot(para, savedweights; binsize=0.2, save=save)
end
