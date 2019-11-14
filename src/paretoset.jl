import Statistics:mean

function ParetoSub(X, perm)
    Z= X.-X[:,1]
    nondominated=reshape(any(Z.<0,dims=1),:)
    ispareto = all(any(Z[:,nondominated].>0, dims=1))
    X=X[:, nondominated]
    nondominated = perm[nondominated]
    return ispareto, nondominated, X
end

"""
    index = ParetoSet(X)
Get the Pareto set from a given set of oberservations `X`.
# Argument
- `X::AbstractMatrix` n_objectives * n_samples

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
