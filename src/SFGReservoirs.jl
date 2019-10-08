module SFGReservoirs

using SparseArrays

include("sfspectra.jl")
include("reservoirs.jl")

export  sfspectrum,
        model_A

end # module
