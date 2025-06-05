using Fresa

domain_pwl = Fresa.Domain(
    x -> JacketReactor.PiecewiseLinearProfile(zeros(length(x.fz)), -ones(length(x.fT)), JacketReactor.T_wmin),
    x -> JacketReactor.PiecewiseLinearProfile(ones(length(x.fz)), ones(length(x.fT)), JacketReactor.T_wmax)
)

p0_pwl = [Fresa.Point(JacketReactor.PiecewiseLinearProfile(fill(0.5, JacketReactor.N), fill(0.0, JacketReactor.N), JacketReactor.T_f), JacketReactor.objective)]

nondominated, population = JacketReactor.solve(p0_pwl, domain_pwl)
