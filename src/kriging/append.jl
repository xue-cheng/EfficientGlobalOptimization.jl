# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function Base.append!(
    krg::Kriging{N,M},
    x::AbstractMatrix,
    y::AbstractMatrix,
) where {N,M}
    xs = krg.xscaler * x
    ys = krg.yscaler * y
    for i in 1:M
        append!(krg.gps[i], xs, ys[i, :])
    end
end
