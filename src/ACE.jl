module ACE

export MaterialStream
export mstream_TQA
export mstream_TpA
export mstream_pQA

export mstream_H_T

export heater_W
export heater_Q

using Clapeyron
using Roots
using NLsolve
using Integrals


include("structs.jl")
include("material.jl")
include("enthalpy.jl")
include("heat.jl")




end # module ACE
