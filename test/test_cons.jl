@testset "Constrained" begin 
    # EI Min
    f = EfficientGlobalOptimization.PressureVesselDesign()
    yₒ, xₒ = optimum(f)
    lb = lowerbounds(f)
    ub = upperbounds(f)
    cons = constrains(f)
    xs = sampling(lb, ub, 9, cons...)
    x0 = Vector{Float64}[]
    y0 = Vector{Float64}[]
    maxg(x) =  maximum(g->g(x) ,cons)
    for x in eachcol(xs)
        if maxg(x) > eps()
            continue
        end
        push!(x0, x)
        push!(y0, [f(x)])
    end
    @info "$(length(y0)) initial samples"
    krg = Kriging(hcat(x0...), hcat(y0...), MeanConst(5e3), Mat32Ard(zeros(4), 3.0))
    kbd = (fill(-10.0, 5), fill(10.0, 5))
    mbd = ([1e2], [1e4])
    tuner = KrigingTuner(kbd, mbd; global_every = 8, local_every = 2)
    ego = EGO(krg, tuner, lb, ub; constrains=cons)
    try
        for i = 1:64
            x, fx = acquire(ego)
            y = f(x)
            append!(ego, x, y)
        end
    catch e
        if isa(e, NoFeasibleInfill)
        else
            rethrow(e)
        end
    end
    yo, xo = optimum(ego)
    @test yo[1] ≈ f(xo)
    @test yo[1]/yₒ < 1.001

end