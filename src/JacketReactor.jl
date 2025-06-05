module JacketReactor

using DifferentialEquations # For solving the system of ODEs
using Fresa                 # For optimisation

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

# PWP constants
global const T_w0min :: Float64 = T_wmin
global const T_w0max :: Float64 = T_wmax

global const T_w1min :: Float64 = T_wmin
global const T_w1max :: Float64 = T_wmax

global const T_w2min :: Float64 = T_wmin
global const T_w2max :: Float64 = T_wmax

global const T_wfmin :: Float64 = T_wmin
global const T_wfmax :: Float64 = T_wmax

global const T_w1prime_min = -1.0
global const T_w1prime_max = 1.0

global const T_w2prime_min = -1.0
global const T_w2prime_max = 1.0
 
global const z1_min :: Float64 = 0.0 * L
global const z1_max :: Float64 = 0.33 * L
global const z2_min :: Float64 = 1/3 * L
global const z2_max :: Float64 = 2/3 * L
global const z3_min :: Float64 = 2/3 * L
global const z3_max :: Float64 = 1.0 * L

abstract type TemperatureProfile end

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

#------------------------------------
# PWL
#------------------------------------
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

#------------------------------------
# PWP
#------------------------------------
mutable struct PiecewisePolynomialProfile <: TemperatureProfile
    T_w0 :: Float64
    T_w1 :: Float64
    T_w2 :: Float64
    T_wf :: Float64
    T_w1prime :: Float64
    T_w2prime :: Float64
    z1 :: Float64
    z2 :: Float64
    z3 :: Float64
    alpha1 :: Float64
    beta1 :: Float64
    gamma1 :: Float64
    zeta1 :: Float64

    alpha2 :: Float64
    beta2 :: Float64
    gamma2 :: Float64
    zeta2 :: Float64

    alpha3 :: Float64
    beta3 :: Float64
    gamma3 :: Float64
    zeta3 :: Float64

    function PiecewisePolynomialProfile(T_w0, T_w1, T_w2, T_wf, T_w1prime, T_w2prime, z1, z2, z3)
        # Define constants
        alpha1 = T_w0
        beta1 = 0
        gamma1 = -(T_w1prime * z1 - 3 * T_w1 + 3 * T_w0) / z1^2
        zeta1 = (T_w1prime * z1 - 2 * T_w1 + 2 * T_w0) / z1^3

        denominator2 = -z2^3 + 3 * z1 * z2^2 - 3 * z1^2 * z2 + z1^3
        alpha2 = -(z1 * (-T_w1prime * z2^3 - 3 * T_w1 * z2^2) + T_w1 * z2^3 + z1^2 * (-T_w2prime * z2^2 + T_w1prime * z2^2 + 3 * T_w2 * z2) + z1^3 * (T_w2prime * z2 - T_w2)) / denominator2
        beta2 = (-T_w1prime * z2^3 + z1 * (-2 * T_w2prime * z2^2 - T_w1prime * z2^2 + 6 * T_w2 * z2 - 6 * T_w1 * z2) + z1^2 * (T_w2prime * z2 + 2 * T_w1prime * z2) + T_w2prime * z1^3) / denominator2
        gamma2 = -(-T_w2prime * z2^2 - 2 * T_w1prime * z2^2 + z1 * (-T_w2prime * z2 + T_w1prime * z2 + 3 * T_w2 - 3 * T_w1) + 3 * T_w2 * z2 - 3 * T_w1 * z2 + (2 * T_w2prime + T_w1prime) * z1^2) / denominator2
        zeta2 = (-T_w2prime * z2 - T_w1prime * z2 + (T_w2prime + T_w1prime) * z1 + 2 * T_w2 - 2 * T_w1) / denominator2

        denominator3 = -z3^3 + 3 * z2 * z3^2 - 3 * z2^2 * z3 + z2^3
        alpha3 = (z2 * (T_w2prime * z3^3 + 3 * T_w2 * z3^2) - T_w2 * z3^3 + z2^2 * (-T_w2prime * z3^2 - 3 * T_wf * z3) + T_wf * z2^3) / denominator3
        beta3 = (-T_w2prime * z3^3 + z2 * (-T_w2prime * z3^2 + 6 * T_wf * z3 - 6 * T_w2 * z3) + 2 * T_w2prime * z2^2 * z3) / denominator3
        gamma3 = -(-2 * T_w2prime * z3^2 + z2 * (T_w2prime * z3 + 3 * T_wf - 3 * T_w2) + 3 * T_wf * z3 - 3 * T_w2 * z3 + T_w2prime * z2^2) / denominator3
        zeta3 = (-T_w2prime * z3 + T_w2prime * z2 + 2 * T_wf - 2 * T_w2) / denominator3


        new(T_w0, T_w1, T_w2, T_wf, T_w1prime, T_w2prime, z1, z2, z3, alpha1, beta1, gamma1, zeta1,
            alpha2, beta2, gamma2, zeta2, alpha3, beta3, gamma3, zeta3)
    end
end

function input(profile :: PiecewisePolynomialProfile, z :: Float64) :: Float64
    if z <= profile.z1
        Tw = profile.alpha1 + profile.beta1 * z + profile.gamma1 * z^2 + profile.zeta1 * z^3
    elseif z > profile.z1 && z <= profile.z2
        Tw = profile.alpha2 + profile.beta2 * z + profile.gamma2 * z^2 + profile.zeta2 * z^3
    else
        Tw = profile.alpha3 + profile.beta3 * z + profile.gamma3 * z^2 + profile.zeta3 * z^3
    end

    if Tw < T_wmin 
        Tw = T_wmin
    elseif Tw > T_wmax
        Tw = T_wmin
    end

    return Tw
end

#------------------------------------
# QSP
#------------------------------------
mutable struct QuadraticSplineProfile <: TemperatureProfile
    T_w0 :: Float64               
    T_wf :: Float64               
    z1 :: Float64               
    a :: Float64                
    b :: Float64
    c :: Float64
    d :: Float64
    e :: Float64
    f :: Float64
    function QuadraticSplineProfile(T_w0, T_wf, z1)
        a = T_w0
        b = 0
        c = (T_wf-T_w0)/z1
        d = (T_wf*z1 - T_w0) / (z1 - 1)
        e = (2*T_w0 - 2*T_wf) / (z1 - 1)
        f = (T_wf - T_w0) / (z1 - 1)
        new(T_w0, T_wf, z1, a, b, c, d, e, f)
    end
end

function input(profile :: QuadraticSplineProfile, z :: Float64) :: Float64
    if z <= profile.z1
        Tw = profile.a + profile.b*z + profile.c*z^2
    else
        Tw = profile.d + profile.e*z + profile.f*z^2
    end
    return Tw
end

function objective_J1_J2(profile :: TemperatureProfile)
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

function objective_J2_J1(profile :: TemperatureProfile)
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

    return ([J2, J1], g)
end


function Fresa.neighbour(x :: PiecewiseLinearProfile, f :: Float64, domain)
    a = domain.lower(x)
    b = domain.upper(x)
    fz = Fresa.neighbour(x.fz, f, Fresa.Domain(x -> a.fz, x -> b.fz))
    fT = Fresa.neighbour(x.fT, f, Fresa.Domain(x -> a.fT, x -> b.fT))
    T0 = Fresa.neighbour(x.T0, f, Fresa.Domain(x -> a.T0, x -> b.T0))
    return PiecewiseLinearProfile(fz, fT, T0)
end

function Fresa.neighbour(x :: PiecewisePolynomialProfile, f :: Float64, domain)
    a = domain.lower(x)
    b = domain.upper(x)
    T_w0 = Fresa.neighbour(x.T_w0, f, Fresa.Domain(x -> a.T_w0, x -> b.T_w0))
    T_w1 = Fresa.neighbour(x.T_w1, f, Fresa.Domain(x -> a.T_w1, x -> b.T_w1))
    T_w2 = Fresa.neighbour(x.T_w2, f, Fresa.Domain(x -> a.T_w2, x -> b.T_w2))
    T_wf = Fresa.neighbour(x.T_wf, f, Fresa.Domain(x -> a.T_wf, x -> b.T_wf))
    T_w1prime = Fresa.neighbour(x.T_w1prime, f, Fresa.Domain(x -> a.T_w1prime, x -> b.T_w1prime))
    T_w2prime = Fresa.neighbour(x.T_w2prime, f, Fresa.Domain(x -> a.T_w2prime, x -> b.T_w2prime))
    z1 = Fresa.neighbour(x.z1, f, Fresa.Domain(x -> a.z1, x -> b.z1))
    z2 = Fresa.neighbour(x.z2, f, Fresa.Domain(x -> a.z2, x -> b.z2))
    z3 = Fresa.neighbour(x.z3, f, Fresa.Domain(x -> a.z3, x -> b.z3))
    PiecewisePolynomialProfile(T_w0, T_w1, T_w2, T_wf, T_w1prime, T_w2prime, z1, z2, z3)
end

function Fresa.neighbour(x :: QuadraticSplineProfile, f, domain)
    a = domain.lower(x)
    b = domain.upper(x)
    T_w0 = Fresa.neighbour(x.T_w0, f, Fresa.Domain(x -> a.T_w0, x -> b.T_w0))
    T_wf = Fresa.neighbour(x.T_wf, f, Fresa.Domain(x -> a.T_wf, x -> b.T_wf))
    z1 = Fresa.neighbour(x.z1, f, Fresa.Domain(x -> a.z1, x -> b.z1))
    QuadraticSplineProfile(T_w0, T_wf, z1)
end

function solve(p0, domain :: Fresa.Domain)
    Fresa.solve(
        # the first 2 arguments are required
        objective,               # the objective function
        p0;                      # an initial point in the design space
        # the rest are option arguments for Fresa
        domain = domain,         # search domain for the decision variables
        archiveelite = false,    # save thinned out elite members
        elite = true,            # elitism by default
        ϵ = 0.001,               # tolerance for similarity detection
        fitnesstype = :hadamard, # how to rank solutions in multi-objective case
        multithreading = true,   # use multiple threads, if available
        ngen = 200,              # number of generations
        #nfmax = 100000,         # number of function evaluations
        np = (40, 100),          # propagation size
        nrmax = 5,               # number of runners maximum
        ns = 100,)               # number of stable solutions for stopping
        #plotvectors = true,
        #populationoutput = true)  
                                     
end

end