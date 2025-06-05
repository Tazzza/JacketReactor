using Fresa

function lower(profile :: JacketReactor.TemperatureProfile)
    if typeof(profile) == JacketReactor.PiecewiseLinearProfile
        JacketReactor.PiecewiseLinearProfile(zeros(length(profile.fz)), -ones(length(profile.fT)), JacketReactor.T_wmin)
    else
        JacketReactor.PiecewisepolynomialProfile(
        JacketReactor.T_w0min, 
        JacketReactor.T_w1min, 
        JacketReactor.T_w2min, 
        JacketReactor.T_wfmin, 
        JacketReactor.T_w1prime_min, 
        JacketReactor.T_w2prime_min, 
        JacketReactor.z1_min, 
        JacketReactor.z2_min, 
        JacketReactor.z3_min)
    end
end

function upper(profile :: JacketReactor.TemperatureProfile)
    if typeof(profile) == JacketReactor.PiecewiseLinearProfile
        JacketReactor.PiecewiseLinearProfile(ones(length(profile.fz)), ones(length(profile.fT)), JacketReactor.T_wmax)
    else
        PiecewisepolynomialProfile(
        JacketReactor.T_w0max, 
        JacketReactor.T_w1max, 
        JacketReactor.T_w2max, 
        JacketReactor.T_wfmax, 
        JacketReactor.T_w1prime_max, 
        JacketReactor.T_w2prime_min, 
        JacketReactor.z1_max, 
        JacketReactor.z2_max, 
        JacketReactor.z3_max)
    end
end

domain = Fresa.Domain(lower, upper)

p0 = [Fresa.Point(JacketReactor.PiecewiseLinearProfile(fill(0.5, JacketReactor.N), fill(0.0, JacketReactor.N), JacketReactor.T_f), 
      JacketReactor.objective),
      Fresa.Point(JacketReactor.PiecewisepolynomialProfile(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25 * JacketReactor.L, 0.5 * JacketReactor.L, 0.75 * JacketReactor.L),
    JacketReactor.objective)]

nondominated, population = Fresa.solve(
        # the first 3 arguments are required
        JacketReactor.objective, # the objective function
        p0;                      # an initial point in the design space
        domain,         # search domain for the decision variables
        # the rest are option arguments for Fresa
        archiveelite = false,    # save thinned out elite members
        elite = true,            # elitism by default
        ϵ = 0.001,               # tolerance for similarity detection
        fitnesstype = :hadamard, # how to rank solutions in multi-objective case
        multithreading = true,   # use multiple threads, if available
        ngen = 1000,             # number of generations
        #nfmax = 100000,         # number of function evaluations
        np = (40, 100),          # propagation size
        nrmax = 5,               # number of runners maximum
        ns = 100,)

