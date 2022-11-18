
mutable struct KrigingTuner
    counter::Int
    g_inter::Int
    l_inter::Int
    g_opt::GlobalOptimizer
    l_opt::LocalOptimizer
    kernbounds
    meanbounds
end

function KrigingTuner(kernbounds,
    meanbounds;
    global_every::Integer = 50,
    global_opt::GlobalOptimizer = ISRESOptimizer(maxeval = 10000),
    local_every::Integer = 5,
    local_opt::LocalOptimizer = LBFGSOptimizer(maxeval = 1000))
    global_every < 0 && throw(ArgumentError("`global_every` must be non-negative"))
    local_every < 0 && throw(ArgumentError("`local_every` must be non-negative"))
    if local_every==global_every==0
        @warn "Both `global_every` and `local_every` are zeros, never tune hyperparameters!"
    end
    if (local_every != 0 && global_every % local_every != 0)
        throw(ArgumentError("`global_every` must be a multiple of `local_every`"))
    end
    KrigingTuner(0, global_every, local_every, global_opt, local_opt,kernbounds, meanbounds)
end

function tunekrg!(krg::Kriging, tuner::KrigingTuner)
    nobs = length(krg.gps[1].y)
    if tuner.g_inter > 0 && ceil(nobs/tuner.g_inter) > ceil(tuner.counter/tuner.g_inter)
        tune!(krg, tuner.kernbounds, tuner.meanbounds, tuner.g_opt)
    end
    if tuner.l_inter > 0 && ceil(nobs/tuner.l_inter) > ceil(tuner.counter/tuner.l_inter)
        tune!(krg, tuner.kernbounds, tuner.meanbounds, tuner.l_opt)
    end
    tuner.counter = nobs
end

function forcetune!(krg::Kriging, tuner::KrigingTuner)
    tune!(krg, tuner.kernbounds, tuner.meanbounds, tuner.g_opt)
    tune!(krg, tuner.kernbounds, tuner.meanbounds, tuner.l_opt)
    tuner.counter = length(krg.gps[1].y)
end