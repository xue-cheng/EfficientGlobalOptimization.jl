import EfficientGlobalOptimization
import NLopt

EGO = EfficientGlobalOptimization
@testset "optimizers" begin
    @testset "stochastic" begin
        f = EGO.Branin_Forrester()
        of = (x)->f(x)
        lb = EGO.lowerbounds(f)
        ub = EGO.upperbounds(f)
        fo, xo = EGO.optimum(f)
        for opt in (EGO.CRSOptimizer, EGO.ISRESOptimizer)
            o = opt(maxeval=10000, population=100)
            fx, x, ret = EGO.minimize(o, of,:Auto, lb, ub, lb)
            @test isapprox(fx, fo, rtol=1e-3)
        end
    end
    @testset "local" begin
        f = EGO.SumSphere(3)
        of = (x)->f(x)
        lb = EGO.lowerbounds(f)
        ub = EGO.upperbounds(f)
        fo, xo = EGO.optimum(f)
        for opt in (EGO.BOBYQAOptimizer, EGO.LBFGSOptimizer)
            o = opt(maxeval=1000)
            fx, x, ret = EGO.minimize(o, of, :Auto, lb, ub, lb)
            @test isapprox(fx, fo, atol=1e-3)
        end
    end
end
EGO = nothing
