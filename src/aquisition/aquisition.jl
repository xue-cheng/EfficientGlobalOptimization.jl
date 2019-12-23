abstract type EGOAquisition end

update_parameters!(a::EGOAquisition, args...) = error("NOT IMPLEMENTED")

import StatsFuns

const φ = StatsFuns.normpdf
const Φ = StatsFuns.normcdf

include("single_obj_aquisitions.jl")
#include("multi_obj_aquisitions.jl")
