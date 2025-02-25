"""
calculates enthalpy [J/kg] of stream
at specific temperature T
"""
function mstream_H_T(s::MaterialStream, T)
    ins = material_stream_copy(s)                                #внутренний временный поток
    A =  @. ins.Q * ins.y + (1.0-ins.Q) * ins.x         #расчет суммарной доли по фазам
    ins = mstream_TpA(ins.G, T, ins.p, A, ins.model)    #температура заменена   
    #println("Q ", ins.Q)
    if 0.0 < ins.Q < 1.0
        #двухфазная система
        M_x_sm = sum(ins.model.params.Mw .* ins.x) / 1000.0         #жидкая фаза
        hl = enthalpy(ins.model, ins.p, ins.T, ins.x) / M_x_sm
		M_y_sm = sum(ins.model.params.Mw .* ins.y) / 1000.0         #газовая фаза
        hv = enthalpy(ins.model, ins.p, ins.T, ins.y) / M_y_sm
        #println("2 ph hl", hl, " hv ", hv, " Q ", ins.Q)
        return ins.Q * hv + (1.0-ins.Q) * hl
    elseif ins.Q == 0.0
        M_sm = sum(ins.model.params.Mw .* ins.x) / 1000.0   #молярная масса смеси [кг/моль]
        #println("vap ", enthalpy(ins.model, ins.p, ins.T, ins.x), " Msm ", M_sm)
        return  enthalpy(ins.model, ins.p, ins.T, ins.x) / M_sm
    elseif ins.Q == 1.0
        M_sm = sum(ins.model.params.Mw .* ins.y) / 1000.0
        #println("liq ", enthalpy(ins.model, ins.p, ins.T, ins.y), " Msm ", M_sm)
        return  enthalpy(ins.model, ins.p, ins.T, ins.y) / M_sm
    end
end


