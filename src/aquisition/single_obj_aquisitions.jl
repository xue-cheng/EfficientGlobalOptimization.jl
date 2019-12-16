

abstract type SingleObjAquisition{S}<:EGOAquisition end

update_parameters!(a::SingleObjAquisition, krg::Kriging) = nothing

include("single_obj_improvement.jl")
