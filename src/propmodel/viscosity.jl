
mutable struct ViscModelDefault <: ViscModel
    p::Float64
end

function viscosity_empty()
    return ViscModelDefault(0.1)
end

function μ_vap()

end

function μ_liq()
    
end
