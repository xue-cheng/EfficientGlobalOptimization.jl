struct ConstrainedAquisition{A<:EGOAquisition,N} <: EGOAquisition
    base::A
    cons::NTuple{N,Function}
end

update_parameters!(a::ConstrainedAquisition, args...) = update_parameters!(a.base, args...)

function objective(a::ConstrainedAquisition, model::Kriging, x)
    for g in a.cons
        if g(x) > eps()
            return 0.0
        end
    end
    return objective(a.base, model, x)
end
