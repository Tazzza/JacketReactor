using Fresa

domain_pwp = Fresa.Domain(
    x -> JacketReactor.PiecewisePolynomialProfile(
        JacketReactor.T_w0min, 
        JacketReactor.T_w1min, 
        JacketReactor.T_w2min, 
        JacketReactor.T_wfmin, 
        JacketReactor.T_w1prime_min, 
        JacketReactor.T_w2prime_min, 
        JacketReactor.z1_min, 
        JacketReactor.z2_min, 
        JacketReactor.z3_min),
    x -> JacketReactor.PiecewisePolynomialProfile(
        JacketReactor.T_w0max, 
        JacketReactor.T_w1max, 
        JacketReactor.T_w2max, 
        JacketReactor.T_wfmax, 
        JacketReactor.T_w1prime_max, 
        JacketReactor.T_w2prime_min, 
        JacketReactor.z1_max, 
        JacketReactor.z2_max, 
        JacketReactor.z3_max)
)

p0_pwp = [JacketReactor.Fresa.Point(JacketReactor.PiecewisePolynomialProfile(
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.25 * JacketReactor.L, 0.5 * JacketReactor.L, 0.75 * JacketReactor.L),
    JacketReactor.objective_J1_J2)]

nondominated, population = JacketReactor.solve(p0_pwp, domain_pwp)
