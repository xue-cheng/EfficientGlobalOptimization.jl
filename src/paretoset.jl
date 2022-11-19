# File: EfficientGlobalOptimization.jl/src/paretoset.jl
# Copyright (c) 2019-2022 XUE Cheng
# 
# This software is released under the MIT License.
# https://opensource.org/licenses/MIT

function ParetoSub(X, perm)
    Z= X.-X[:,1]
    nondominated=reshape(any(Z.<0,dims=1),:)
    ispareto = all(any(Z[:,nondominated].>0, dims=1))
    X=X[:, nondominated]
    nondominated = perm[nondominated]
    return ispareto, nondominated, X
end

"""
    index = ParetoSet(X[, sense])
Get the Pareto set from a given set of oberservations `X`. 
# Argument
- `X::AbstractMatrix{Real}` size(X)=(n_objectives, n_samples)
- `sense::AbstractVector{Union{Symbol, EGOSense}}=repeat(:Min, n_objectives)` sense of each objective, size(S)=(n_objectives,)
"""
function ParetoSet(X::AbstractMatrix)
    n,m = size(X)
    Xmin = minimum(X, dims=2)
    Xt = X.-Xmin
    Xm = mean(Xt, dims=2)
    Xn = maximum(Xt ./ (Xm.+maximum(Xm)), dims=1)
    perm = sortperm(reshape(Xn,:))
    Y = X[:, perm]
    membership = falses(m)
    while length(perm)>1
        k = perm[1]
        membership[k], perm, Y = ParetoSub(Y, perm)
    end
    membership[perm] .= true
    membership
end

ParetoSet(X, sense::AbstractVector{EGOSense}) = ParetoSet(X.*Int.(sense))

ParetoSet(X, sense::AbstractVector{Symbol}) = ParetoSet(X, sym2sense.(sense))
