"""
calculates density [kg/m^3] of stream
"""
function mstream_D(s::MaterialStream)
    if 0 < s.Q < 1
        ρᵥ = mass_density(s.model, s.p, s.T, s.y)
        ρₗ = mass_density(s.model, s.p, s.T, s.x)
        return 1.0 / ((1.0 - s.Q) / ρₗ + s.Q / ρᵥ)
    elseif s.Q == 1.0
        return mass_density(s.model, s.p, s.T, s.y)
    elseif s.Q == 0.0
        return mass_density(s.model, s.p, s.T, s.x)
    end
end