

"""
Flow properties
G -- mass flow [kg/s]
T -- temperature [K]

"""
mutable struct MaterialStream
    G::Float64          #mass flow
    T::Float64          #temperature
    p::Float64          #pressure
    x::AbstractArray    #mol frac liquid
    y::AbstractArray    #mol frac vapor
    Q::Float64          #vapor mole frac
    model::EoSModel     #thermodynamic model
end

"material stream copy"
function material_stream_copy(s::MaterialStream)
    return  Stream(s.G, s.T, s.p, s.x, s.y, s.Q, s.model)
end

