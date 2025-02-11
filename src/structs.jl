"""
Фаза
надо будет подумать
"""
mutable struct Phase

end


"""
Flow properties
G -- mass flow [kg/s]
T -- temperature [K]

"""
mutable struct MaterialStream
    G::Float64          #mass flow
    T::Float64          #temperature
    p::Float64          #pressure
    x::AbstractArray    #mol frac in total flow
    y::AbstractArray    #mol frac vapor
    Q::Float64          #vapor mole frac
    model::EoSModel     #thermodynamic model
end

"material stream copy"
function material_stream_copy(s::MaterialStream)
    return  MaterialStream(s.G, s.T, s.p, s.x, s.y, s.Q, s.model)
end

"""

"""
mutable struct EnergyStream
    Q::Float64  #тепловой поток [Вт]
end

