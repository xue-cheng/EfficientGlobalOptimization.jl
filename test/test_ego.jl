# File: EfficientGlobalOptimization.jl/test/test_ego.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

@testset "serial ego" begin
    # EI Min
    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    x0 = sampling(lb, ub, 7)
    y0 = transpose(map(f, eachcol(x0)))
    krg = Kriging(x0, y0, MeanConst(0.0), SEArd([0.0,0.0], 3.0))
    kbd = (fill(-5.0, 3), fill(5.0, 3))
    mbd = ([-1.0], [1.0])
    tuner = KrigingTuner(kbd, mbd; global_every = 8, local_every = 2)
    ego = EGO(krg, tuner, lb, ub)
    try
        for i = 1:30
            x, fx = acquire(ego)
            y = f(x)
            append!(ego, x, y)
        end
    catch e
        if isa(e, NoFeasibleInfill)
            @info e.msg
        else
            rethrow(e)
        end
    end
    yo, xo = optimum(ego)
    @test yo[1] ≈ yₒ rtol = 0.001
end


@testset "parallel ego" begin

    f = EfficientGlobalOptimization.Branin_Forrester()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    x0 = sampling(lb, ub, 7)
    y0 = transpose(map(f, eachcol(x0)))
    krg = Kriging(x0, y0, MeanConst(0.0), SEArd([0.0,0.0], 3.0))
    kbd = (fill(-5.0, 3), fill(5.0, 3))
    mbd = ([-1.0], [1.0])
    tuner = KrigingTuner(kbd, mbd; global_every = 8, local_every = 2)
    ego = EGO(krg, tuner, lb, ub)
    try
        for j = 1:6
            newx = []
            for i = 1:5
                x, fx = acquire(ego, newx...)
                push!(newx, x)
            end
            if isempty(newx)
                break
            end
            newy = map(f, newx)
            append!(ego, hcat(newx...), newy)
        end
    catch e
        if isa(e, NoFeasibleInfill)
            @info e.msg
        else
            rethrow(e)
        end
    end
    yo, xo = optimum(ego)
    @test yo[1] ≈ yₒ rtol = 0.001
end