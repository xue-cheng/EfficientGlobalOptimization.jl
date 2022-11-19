# File: EfficientGlobalOptimization.jl/src/sampling/sampling.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function sampling(lb::AbstractVector{T}, ub::AbstractVector{T}, ns::Int, gs::Function...) where T
    @assert length(lb) == length(ub)
    @assert all(lb .< ub)
    ndim = length(lb)
    s = SobolSeq(ndim)
    x = Matrix{T}(undef, ndim, ns)
    scalar = MinMaxScaler(lb, ub)
    for i in 1:ns
        while true
            xx = inverse!(scalar,  next!(s))
            if all(map(g->g(xx)<=0, gs))
                x[:, i] .= xx
                break
            end
        end
    end
    return x
end
