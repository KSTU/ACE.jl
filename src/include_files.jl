include("structs.jl")
include("support.jl")
include("errors.jl")

include("propmodel/viscosity.jl")
include("propmodel/thermal.jl")

include("properties/material_stream/pQA.jl")
include("properties/material_stream/TpA.jl")

#TESTS
include("testing/td_tests.jl")

#include("errors.jl")

#include("properties/enthalpy.jl")
#include("properties/density.jl")
#include("properties/intenergy.jl")

#include("apps/heat.jl")
#include("apps/mixer.jl")
#include("apps/splitter.jl")
#include("apps/pressure_change.jl")
