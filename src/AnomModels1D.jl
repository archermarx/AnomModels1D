module AnomModels1D

using HallThruster

include("models/lafleur.jl")
include("models/datadriven.jl")

export LafleurModel
export DataDriven1


end # module AnomModels1D
