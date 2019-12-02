@testset "serial ego" begin
    # EI Min
    f = EGO.Branin_Forrester()
    yₒ, xₒ = EGO.optimum(f)
    lb = EGO.lowerbounds(f)
    ub = EGO.upperbounds(f)
    ego = EGO.EGO(EGO.ExpectedImprovement(:Min), lb, ub)
    EGO.meanbounds!(ego, [-20., 300.])
    x₀ = EGO.init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = EGO.acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = EGO.optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # PI Min
    f = EGO.Branin_Forrester()
    yₒ, xₒ = EGO.optimum(f)
    lb = EGO.lowerbounds(f)
    ub = EGO.upperbounds(f)
    ego = EGO.EGO(EGO.ProbabilityOfImprovement(:Min), lb, ub)
    EGO.meanbounds!(ego, [-20., 20.])
    x₀ = EGO.init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = EGO.acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = EGO.optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # EI Max
    bf = EGO.Branin_Forrester()
    f(x) = -bf(x)
    yₒ, xₒ = EGO.optimum(bf)
    yₒ = -yₒ
    lb = EGO.lowerbounds(bf)
    ub = EGO.upperbounds(bf)
    ego = EGO.EGO(EGO.ExpectedImprovement(:Max), lb, ub)
    EGO.meanbounds!(ego, [-20., 20.])
    x₀ = EGO.init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = EGO.acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = EGO.optimum(ego)
    @test yo ≈ yₒ rtol=0.001
    # PI Max
    bf = EGO.Branin_Forrester()
    f(x) = -bf(x)
    yₒ, xₒ = EGO.optimum(bf)
    yₒ = -yₒ
    lb = EGO.lowerbounds(bf)
    ub = EGO.upperbounds(bf)
    ego = EGO.EGO(EGO.ProbabilityOfImprovement(:Max), lb, ub)
    EGO.meanbounds!(ego, [-20., 20.])
    x₀ = EGO.init_sampling(ego, 7) # 3*N+1
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:30
        x = EGO.acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = EGO.optimum(ego)
    @test yo ≈ yₒ rtol=0.001
end


@testset "parallel ego" begin
    f = EGO.Branin_Forrester()
    yₒ, xₒ = EGO.optimum(f)
    lb = EGO.lowerbounds(f)
    ub = EGO.upperbounds(f)
    ego = EGO.EGO(EGO.ExpectedImprovement(:Min), lb, ub)
    EGO.meanbounds!(ego, [-20., 20.])
    x₀ = EGO.init_sampling(ego, 7)
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    newx = Matrix{Float64}(undef, ndims(ego,1), 5)
    for j = 1:6
        for i = 1:5
            newx[:, i] = EGO.acquire(ego, newx[:, 1:i-1])
        end
        newy = map(f, eachcol(newx))
        append!(ego, newx, newy)
    end
    yo, xo = EGO.optimum(ego)
    @test yo ≈ yₒ rtol=0.001
end