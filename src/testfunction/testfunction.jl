import Base: ndims

abstract type TestFunction{N,M} end
ndims(::TestFunction{N,M}) where{N,M} = (N,M)
include("single_obj_functions.jl")
