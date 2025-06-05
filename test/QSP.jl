using Fresa

domain = Fresa.Domain(
    x -> JacketReactor.QuadraticSplineProfile(JacketReactor.T_wmin, JacketReactor.T_wmax, 0.25),
    x -> JacketReactor.QuadraticSplineProfile(JacketReactor.T_wmin, JacketReactor.T_wmax, 0.75)
)

p0 = [Fresa.Point(JacketReactor.QuadraticSplineProfile(JacketReactor.T_f, JacketReactor.T_f, 0.5),
                  JacketReactor.objective_J1_J2)]

nondominated, population = JacketReactor.solve(p0, domain)
