using Random:randperm

@enum SamplingMethods LHS GRID 
_sym2sm = Dict(Symbol(i) => i for i in instances(SamplingMethods))
function sym2smethod(s::Symbol)
    try
        return _sym2sm[s]
    catch
        error("invalid sampling method `$s`")
    end
end

sampling(s::Symbol, args...) = sampling(sym2smethod(s), args...)

sampling(s::SamplingMethods, lb::AbstractVector{T}, ub::AbstractVector{T}, ns) where {T} = sampling(Val(s), lb, ub, ns)

function sampling(::Val{LHS},
    lb::AbstractVector{T},
    ub::AbstractVector{T},
    ns::Int) where {T<:AbstractFloat} 
    ndim = length(lb)
    length(ub) != ndim && throw(ArgumentError("length of `lb` and `ub` must be same"))
    s = similar(lb, ndim, ns)
    @inbounds for i in 1:ndim
        s[i,:] .=  (rand(ns) .+ randperm(ns) .- 1)./ns
    end
    scaler = MinMaxScaler(lb,ub) 
    KrigingModel.inverse!(scaler, s)
    s
end

function sampling(::Val{GRID},
    lb::AbstractVector{T},
    ub::AbstractVector{T},
    ns::NTuple{N,Int}) where {N,T}
    length(lb)!=N && throw(ArgumentError("length of `lb` must be $N(length of `ns`)"))
    length(ub)!=N && throw(ArgumentError("length of `ub` must be $N(length of `ns`)"))
    any(ns.<1) && throw(ArgumentError("`ns`=$ns: at least 1 sample in each dimension"))
    ntot = prod(ns)
    s = similar(lb, N, ntot)
    d = (ub-lb)./ns
    iters = Iterators.product(((lb[i]+d[i]/2):d[i]:ub[i] for i in 1:N)...)
    @inbounds @views for (i,p) in zip(1:ntot, iters)
        s[:,i] .= p 
    end
    s
end