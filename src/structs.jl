"""
Flow properties
G -- mass flow [kg/s]
T -- temperature [K]

"""
mutable struct Stream
    G::Float64          #mass flow
    T::Float64          #temperature
    p::Float64          #pressure
    x::AbstractArray    #mol frac liquid
    y::AbstractArray    #mol frac vapor
    Q::Float64          #vapor mole frac
    model::TDmodel      #thermodynamic model
end