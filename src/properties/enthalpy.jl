"""
Находит удельную мольную энтальпию потока

"""
function mstream_H(ms::MaterialStream)
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x         #расчет суммарной доли по фазам 
    if 0.0 < ms.Q < 1.0    #двухфазная система
        hl = enthalpy(ms.model, ms.p, ms.T, ms.x, phase=:liquid)
        hv = enthalpy(ms.model, ms.p, ms.T, ms.y, phase=:vapor)
        return ms.Q * hv + (1.0-ms.Q) * hl
    elseif ms.Q == 0.0
        return  enthalpy(ms.model, ms.p, ms.T, ms.x, phase=:liquid)
    elseif ms.Q == 1.0
        return  enthalpy(ms.model, ms.p, ms.T, ms.y, phase=:vapor)
    end
end

function mstream_H_kJkg(ms::MaterialStream)
    h = mstream_H(ms)    #Дж на моль
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x
    Mₛ = sum(ms.model.params.Mw .* A) / 1000.0      #молярная масса смеси кг /моль
    return h / Mₛ / 1000
end


"""
calculates enthalpy [J/kg] of stream
at specific temperature T
"""
function mstream_H_T(s::MaterialStream, T)
    A =  @. s.Q * s.y + (1.0-s.Q) * s.x         #расчет суммарной доли по фазам
    ins = mstream_TpA(s.N, T, s.p, A, s.model)    #температура заменена   
    return mstream_H(ins)
end

function mstream_H_Tp(ms::MaterialStream, T, p)
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x         #расчет суммарной доли по фазам
    ins = mstream_TpA(ms.N, T, p, A, ms.model)    #температура заменена
    return mstream_H(ins)
end

function mstream_H_pQ(ms::MaterialStream, p, Q)
    A =  @. ms.Q * ms.y + (1.0-ms.Q) * ms.x         #расчет суммарной доли по фазам
    ins = mstream_pQA(ms.N, p, Q, A, ms.model)  #(ms.N, T, p, A, ms.model)    #температура заменена
    return mstream_H(ins)
end
