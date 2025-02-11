"""
calculates density [kg/m^3] of stream
"""
function mstream_D_(s::MaterialStream)
    if 0 < s.Q < 1
        ρᵥ = mass_density(s.model, s.p, s.T, s.y)
        ρₗ = mass_density(s.model, s.p, s.T, s.x)
        
    elseif s.Q == 1.0

    elseif s.Q == 0.0

    end
end