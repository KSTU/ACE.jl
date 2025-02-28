module ACE

export EnergyStream

export MaterialStream
export mstream_TQA
export mstream_TpA
export mstream_pQA

export mstream_H_T

export heater_W
export heater_Q

export splitter_G
export splitter_dol

export mixer

using Clapeyron
using Roots
using NLsolve
using Integrals


include("structs.jl")
include("material.jl")
include("properties/enthalpy.jl")
include("properties/density.jl")
include("heat.jl")
include("mixer.jl")
include("splitter.jl")




end # module ACE
