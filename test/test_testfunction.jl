using MAT

@testset "testfunctions" begin
@testset "single-obj" begin
    function test_optim(test_f, args...)
        f = test_f(args...)
        y, x = EfficientGlobalOptimization.optimum(f)
        abs(f(x)-y) < 1e-6
    end
    function test_inbounds(test_f, args...)
        f = test_f(args...)
        y, x = EfficientGlobalOptimization.optimum(f)
        lb = EfficientGlobalOptimization.lowerbounds(f)
        ub = EfficientGlobalOptimization.upperbounds(f)
        all(lb .<=x .<= ub)
    end
    @test test_inbounds(EfficientGlobalOptimization.Branin_Forrester)
    @test test_optim(EfficientGlobalOptimization.Branin_Forrester)

    for i = 2:6
        @test test_optim(EfficientGlobalOptimization.Ackley, i)
        @test test_optim(EfficientGlobalOptimization.Rosenbrock, i)
        @test test_optim(EfficientGlobalOptimization.SumSphere, i)
        @test test_inbounds(EfficientGlobalOptimization.Ackley, i)
        @test test_inbounds(EfficientGlobalOptimization.Rosenbrock, i)
        @test test_inbounds(EfficientGlobalOptimization.SumSphere, i)
    end
end

@testset "mult-obj" begin
    x = range(0.0, 1.0; length=11)
    zdt1 = EfficientGlobalOptimization.ZDT(1,3)
    zdt2 = EfficientGlobalOptimization.ZDT(2,3)
    zdt3 = EfficientGlobalOptimization.ZDT(3,3)
    nx = length(x)
    y1 = zeros(nx,nx,nx,2)
    y2 = zeros(nx,nx,nx,2)
    y3 = zeros(nx,nx,nx,2)
    for k=1:nx
        for j=1:nx
            for i=1:nx
                xx = [x[i],x[j],x[k]]
                y1[i,j,k,:] = zdt1(xx)
                y2[i,j,k,:] = zdt2(xx)
                y3[i,j,k,:] = zdt3(xx)
            end
        end
    end
    f = matopen("ZDT.mat")
    y1t = read(f, "zdt1")
    y2t = read(f, "zdt2")
    y3t = read(f, "zdt3")
    @test y1 ≈ y1t
    @test y2 ≈ y2t
    @test y3 ≈ y3t
end
end
