using EfficientGlobalOptimization

@testset "serial ego" begin
    # EI Min
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    ego = EGO(ExpectedImprovement(:Min), lb, ub)
    x₀ = init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # PI Min
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    ego = EGO(ProbabilityOfImprovement(:Min), lb, ub)
    x₀ = init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # EI Max
    bf = EfficientGlobalOptimization.Branin_Forrester()
    f(x) = -bf(x)
    yₒ, xₒ = optimum(bf)
    yₒ = -yₒ
    lb = lowerbounds(bf)
    ub = upperbounds(bf)
    ego = EGO(ExpectedImprovement(:Max), lb, ub)
    x₀ = init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # PI Max
    bf = EfficientGlobalOptimization.Branin_Forrester()
    f(x) = -bf(x)
    yₒ, xₒ = optimum(bf)
    yₒ = -yₒ
    lb = lowerbounds(bf)
    ub = upperbounds(bf)
    ego = EGO(ProbabilityOfImprovement(:Max), lb, ub)
    x₀ = init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.001
end


@testset "parallel ego" begin
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    ego = EGO(ExpectedImprovement(:Min), lb, ub)
    x₀ = init_sampling(ego, 7)
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    newx = Matrix{Float64}(undef, ndims(ego), 5)
    for j = 1:6
        for i = 1:5
            newx[:, i] = acquire(ego, newx[:, 1:i-1])
        end
        newy = map(f, eachcol(newx))
        append!(ego, newx, newy)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.001
end