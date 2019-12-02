using MAT

@testset "testfunctions" begin
@testset "single-obj" begin
    function test_optim(test_f, args...)
        f = test_f(args...)
        y, x = EGO.optimum(f)
        abs(f(x)-y) < 1e-6
    end
    function test_inbounds(test_f, args...)
        f = test_f(args...)
        y, x = EGO.optimum(f)
        lb = EGO.lowerbounds(f)
        ub = EGO.upperbounds(f)
        all(lb .<=x .<= ub)
    end
    @test test_inbounds(EGO.Branin_Forrester)
    @test test_optim(EGO.Branin_Forrester)

    for i = 2:6
        @test test_optim(EGO.Ackley, i)
        @test test_optim(EGO.Rosenbrock, i)
        @test test_optim(EGO.SumSphere, i)
        @test test_inbounds(EGO.Ackley, i)
        @test test_inbounds(EGO.Rosenbrock, i)
        @test test_inbounds(EGO.SumSphere, i)
    end
end

@testset "mult-obj" begin
    x = range(0.0, 1.0; length=11)
    zdt1 = EGO.ZDT(1,3)
    zdt2 = EGO.ZDT(2,3)
    zdt3 = EGO.ZDT(3,3)
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
