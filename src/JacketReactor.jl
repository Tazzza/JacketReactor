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

# PWL constants
global const N :: UInt64 = 50
global const ΔTw_max :: Float64 = 1.0

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

abstract type TemperatureProfile end

mutable struct PiecewiseLinearProfile <: TemperatureProfile
    fz :: Vector{Float64}
    fT :: Vector{Float64}
    T0 :: Float64
    function PiecewiseLinearProfile(fz, fT, T0)
        if length(fz) != length(fT)
            error("Length mismatch between fz and fT")
        end
        if T0 < T_wmin || T0 > T_wmax
            error("T0=$(T0) must be in [$T_wmin,$T_wmax]")
        end
        new(fz, fT, T0)
    end
end

function input(profile :: PiecewiseLinearProfile, z :: Float64) :: Float64
    Δτ = L / N
    τ = 0.0
    Tw = profile.T0

    for i in 1:N
        ΔTw = profile.fT[i] * ΔTw_max
        if z < τ + Δτ
            Tw += ΔTw * (z - τ) / Δτ
            Tw = clamp(Tw, T_wmin, T_wmax)
            break
        else
            Tw += ΔTw
            Tw = clamp(Tw, T_wmin, T_wmax)
        end
        τ += Δτ
    end

    return Tw
end

function objective(profile :: TemperatureProfile)
    results = simulation(profile)

    x_vals = results.u

    x_end = x_vals[end]
    J1 = c_f * (1 - x_end[1])
    J2 = (T_f * x_end[2]^2) / K

    g = 0.0
    for i in 1:length(results.u)
        Tw = input(profile, results.t[i])
        x = x_vals[i]
        x2 = x[2]

        g = max(
            g,
            Tw - T_wmax,
            T_wmin - Tw,
            x2 - (T_max / T_f - 1),
            (T_min / T_f - 1) - x2
        )
    end

    return ([J1, J2], g)
end

function Fresa.neighbour(x :: PiecewiseLinearProfile, f :: Float64, domain)
    a = domain.lower(x)
    b = domain.upper(x)
    fz = Fresa.neighbour(x.fz, f, Fresa.Domain(x->a.fz, x->b.fz))
    fT = Fresa.neighbour(x.fT, f, Fresa.Domain(x->a.fT, x->b.fT))
    T0 = Fresa.neighbour(x.T0, f, Fresa.Domain(x->a.T0, x->b.T0))
    return PiecewiseLinearProfile(fz, fT, T0)
end


function solve(p0, domain :: Fresa.Domain)
    Fresa.solve(
        # the first 4 arguments are required
        objective,               # the objective function
        p0;                      # an initial point in the design space
        domain = domain,         # search domain for the decision variables
        # the rest are option arguments for Fresa
        archiveelite = false,    # save thinned out elite members
        elite = true,            # elitism by default
        ϵ = 0.001,               # tolerance for similarity detection
        fitnesstype = :hadamard, # how to rank solutions in multi-objective case
        multithreading = true,   # use multiple threads, if available
        ngen = 200,              # number of generations
        #nfmax = 100000,
        np = (40, 100),          # propagation size
        nrmax = 5,               # number of runners maximum
        ns = 100,)
        #plotvectors = true,
        #populationoutput = true)  
                               # number of stable solutions for stopping      
end

domain_pwl = Fresa.Domain(
    x -> PiecewiseLinearProfile(zeros(length(x.fz)), -ones(length(x.fT)), T_wmin),
    x -> PiecewiseLinearProfile(ones(length(x.fz)), ones(length(x.fT)), T_wmax)
)

p0_pwl = [Fresa.Point(PiecewiseLinearProfile(fill(0.5, N), fill(0.0, N), T_f), objective)]

nondominated, population = solve(p0_pwl, domain_pwl)
