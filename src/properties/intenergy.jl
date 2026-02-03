"""
Calculate molar internal energy [J/mol] of stream
"""
function mstream_U(ms::MaterialStream)
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x         #расчет суммарной доли по фазам 
    if 0.0 < ms.Q < 1.0    #двухфазная система
        hl = internal_energy(ms.model, ms.p, ms.T, ms.x, phase=:liquid)
        hv = internal_energy(ms.model, ms.p, ms.T, ms.y, phase=:vapor)
        return ms.Q * hv + (1.0-ms.Q) * hl
    elseif ms.Q == 0.0
        return  internal_energy(ms.model, ms.p, ms.T, ms.x, phase=:liquid)
    elseif ms.Q == 1.0
        return  internal_energy(ms.model, ms.p, ms.T, ms.y, phase=:vapor)
    end
end

