abstract type ViscModel end
abstract type TermModel end

"""
Properties of stream
"""
mutable struct Props
    td::EoSModel
    visc::ViscModel
    term::TermModel
end

"""
Flow properties:
    N - mole flow [mol/s]
    T - temperature [K]
    p - pressure [Pa]
    x - liquid composition [mol frac]
    y - vapor composition [mol frac]
    Q - vapor moe fraction
    model - thermodynamic model
    name - strem name
"""
mutable struct MaterialStream
    N::Float64          #mole flow  [mol/s]
    T::Float64          #temperature  [K]
    p::Float64          #pressure  [Pa]
    x::AbstractArray    #liquid mole fraction
    y::AbstractArray    #vapor mole fraction
    Q::Float64          #vapor mole frac
    prop::Props     #thermodynamic model
    name::String        #name of stream
end

"Make copy of material stream"
function material_stream_copy(s::MaterialStream)
    return  MaterialStream(s.N, s.T, s.p, s.x, s.y, s.Q, s.model, s.name)
end

"""
Energy stream
"""
mutable struct EnergyStream
    Q::Float64  #Enegry stream [W]
end



