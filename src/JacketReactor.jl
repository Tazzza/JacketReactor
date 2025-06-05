# Import packages
using DifferentialEquations # For solving the system of ODEs
using Plots                 # For plotting
using Fresa                 # For optimisation
using Statistics            # For ordering arrays of objective function values
using Printf                # Saving and formatting output

#------------------------------------
# Constants in differential equations
#------------------------------------
global const c_f :: Float64 = 0.02
global const T_f :: Float64 = 340.0
global const T_min :: Float64 = 280.0
global const T_max :: Float64 = 400.0
global const T_wmin :: Float64 = 280.0
global const T_wmax :: Float64 = 400.0
global const v :: Float64 = 0.1
global const α :: Float64 = 0.0582
global const β :: Float64 = 0.2
global const γ :: Float64 = 16.659
global const δ :: Float64 = 0.25
global const L :: Float64 = 1.0
global const K :: Float64 = 250000.0

function reactor(x, p, z)
    Tw = input(p, z)
    term1 = (α/v) * (1 - x[1]) * exp(γ * x[2] / (1 + x[2]))
    term2 = (β/v) * ((Tw - T_f) / T_f - x[2])

    return [
        term1,
        term1 * δ + term2
    ]
end

function simulation(profile :: TemperatureProfile)
    x0 = [0, 0]  # Initial conditions
    length = (0.0, L)
    prob = ODEProblem(reactor, x0, length, profile)
    results = DifferentialEquations.solve(prob)
    return results
end