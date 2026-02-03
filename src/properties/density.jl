"""
Calculate mass density [kg/m^3] of stream
"""
function mstream_D_mass(ms::MaterialStream)
    if 0 < ms.Q < 1
        ρᵥ = mass_density(ms.model, ms.p, ms.T, ms.y; phase = :vapor)
        ρₗ = mass_density(ms.model, ms.p, ms.T, ms.x; phase = :liquid)
        Mᵥ = ms_get_vapor_molar_mass(ms)
        Mₗ = ms_get_liquid_molar_mass(ms)
        return (Mᵥ * ms.Q + Mₗ * (1-ms.Q)) /(Mᵥ * ms.Q / ρᵥ + Mₗ * (1-ms.Q) / ρₗ)  #ρᵥ , ρₗ, Mᵥ , Mₗ
    elseif ms.Q == 1.0
        return mass_density(ms.model, ms.p, ms.T, ms.y; phase=:vapor)
    elseif ms.Q == 0.0
        return mass_density(ms.model, ms.p, ms.T, ms.x; phase=:liquid)
    end
end

"""
Calculate molar density [mol/m^3] of stream
"""
function mstream_D_mole(ms::MaterialStream)
    D = mstream_D_mass(ms)
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x
    Mₛ = sum(ms.model.params.Mw .* A) / 1000.0      #молярная масса смеси кг / моль
    #println("A ", A, " Mₛ ", Mₛ, " D ", D)
    return D / Mₛ
end

"""
calculate density at a given pressure
"""
function mstream_D_T(ms::MaterialStream, p::Float64)
    ms_in = copy(ms)
    ms_in.p = p
    return mstream_D_mass(ms_in)
end

