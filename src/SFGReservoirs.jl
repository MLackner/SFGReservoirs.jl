module SFGReservoirs

using SparseArrays
using LinearAlgebra
using DifferentialEquations

include("sfspectra.jl")
include("reservoirs.jl")

export sfspectrum
export model

end # module
