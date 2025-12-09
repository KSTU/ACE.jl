
function molar_mass_vapor(ms)
    mw = collect(ms.model.params.Mw)
    return sum(ms.y .* mw)
end

function molar_mass_liquid(ms)
    mw = collect(ms.model.params.Mw)
    return sum(ms.x .* mw)
end

