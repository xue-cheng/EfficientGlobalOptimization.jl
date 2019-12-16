@testset "serial ego" begin
    # EI Min
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    x0 = sampling(:LHS, lb, ub, 7)
    y0 = transpose(map(f, eachcol(x0)))
    krg = Kriging(x0, y0, MeanConst(0.0), SEArd([0.0,0.0], 3.0))
    kbd = (fill(-5.0, 3), fill(5.0, 3))
    mbd = ([-1.0], [1.0])
    tuner = KrigingTuner(kbd, mbd; global_every = 8, local_every = 2)
    ego = EGO(krg, tuner, lb, ub)
    for i = 1:30
        x, fx = acquire(ego)
        @info "$i: $x (EI=$fx)"
        if fx < eps()
            @info "EI ≈ 0, stop!"
            break
        end
        y = f(x)
        append!(ego, x, y)
    end
    yo, xo = optimum(ego)
    @test yo[1] ≈ yₒ rtol = 0.001
end


@testset "parallel ego" begin

    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    x0 = sampling(:LHS, lb, ub, 7)
    y0 = transpose(map(f, eachcol(x0)))
    krg = Kriging(x0, y0, MeanConst(0.0), SEArd([0.0,0.0], 3.0))
    kbd = (fill(-5.0, 3), fill(5.0, 3))
    mbd = ([-1.0], [1.0])
    tuner = KrigingTuner(kbd, mbd; global_every = 8, local_every = 2)
    ego = EGO(krg, tuner, lb, ub)
    for j = 1:6
        newx = []
        for i = 1:5
            x, fx = acquire(ego, newx...)
            @info "$j, $i: $x (EI=$fx)"
            if fx < eps()
                @info "EI ≈ 0, stop!"
                break
            end
            push!(newx, x)
        end
        if isempty(newx)
            break
        end
        newy = map(f, newx)
        append!(ego, hcat(newx...), newy)
    end
    yo, xo = optimum(ego)
    @test yo[1] ≈ yₒ rtol = 0.001
end