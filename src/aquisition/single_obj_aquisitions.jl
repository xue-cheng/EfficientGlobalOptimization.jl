

abstract type SingleObjAquisition{S}<:EGOAquisition end

update_parameters!(a::SingleObjAquisition, gp::GP.GPBase) = nothing

include("single_obj_improvement.jl")
