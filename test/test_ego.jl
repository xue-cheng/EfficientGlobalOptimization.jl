using EfficientGlobalOptimization

@testset "serial ego" begin
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    ego = EGO(ExpectedImprovement(:Min), lb, ub)
    x₀ = init_sampling(ego, 4)
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    for i = 1:28
        x = acquire(ego)
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.0001
end


@testset "parallel ego" begin
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    ego = EGO(ExpectedImprovement(:Min), lb, ub)
    x₀ = init_sampling(ego, 3)
    y₀ = map(f, eachcol(x₀))
    append!(ego, x₀, y₀)
    newx = Matrix{Float64}(undef, ndims(ego), 4)
    for j = 1:7
        for i = 1:4
            newx[:, i] = acquire(ego, newx[:, 1:i-1])
        end
        newy = map(f, eachcol(newx))
        append!(ego, newx, newy)
    end
    yo, xo = optimum(ego)
    @test yo ≈ yₒ rtol=0.0001
end