

function ms_get_liquid_density(ms)
    if ms.Q < 1
        return mass_density(ms.model, ms.p, ms.T, ms.x; phase = :liquid)
    else
        println("нет жидкой фазы")
    end
end

function ms_get_vapor_density(ms)
    if ms.Q > 0
        return mass_density(ms.model, ms.p, ms.T, ms.y; phase = :vapor)
    else
        println("нет паровой фазы")
    end
end

function ms_get_liquid_entalpy(ms)
    if ms.Q < 1
        return mass_enthalpy(ms.model, ms.p, ms.T, ms.x; phase = :liquid)
    else
        println("нет жидкой фазы")
    end
end

function ms_get_vapor_entalpy(ms)
    if ms.Q > 0
        return mass_enthalpy(ms.model, ms.p, ms.T, ms.y; phase = :vapor)
    else
        println("нет паровой фазы")
    end
end