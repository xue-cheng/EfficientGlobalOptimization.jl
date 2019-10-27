using Test
import EfficientGlobalOptimization: Branin_Forrester, Ackley, Rosenbrock
import EfficientGlobalOptimization: optimum, lowerbounds, upperbounds

@testset "testfunctions" begin
    function test_optim(test_f, args...)
        f = test_f(args...)
        y, x = optimum(f)
        abs(f(x)-y) < 1e-6
    end
    function test_inbounds(test_f, args...)
        f = test_f(args...)
        y, x = optimum(f)
        lb = lowerbounds(f)
        ub = upperbounds(f)
        all(lb .<=x .<= ub)
    end
    @test test_inbounds(Branin_Forrester)
    @test test_optim(Branin_Forrester)

    for i = 2:6
        @test test_optim(Ackley, i)
        @test test_optim(Rosenbrock, i)
        @test test_inbounds(Ackley, i)
        @test test_inbounds(Rosenbrock, i)
    end
end